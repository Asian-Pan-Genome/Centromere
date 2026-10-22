"""Human satellite annotation helpers. BED coordinates are 0-based, half-open.

The scientific heuristic is retained from the source workflow; coordinate and
edge-case fixes are itemized in CHANGELOG.md. Reference boundary flanks are mapped
only to define a search window; sample satellite annotations determine boundaries.
"""
import bisect
import collections
import hashlib
import json
import re
from pathlib import Path

COLORS = {"ASat": "137,22,64", "HSat1A": "24,152,139",
          "HSat1B": "24,152,139", "HSat2": "50,51,102",
          "HSat3": "119,160,186", "BSat": "220,174,208", "GSat": "117,134,96"}


def validate_config(cfg):
    samples = cfg.get("samples", {})
    if not isinstance(samples, dict) or not samples:
        raise ValueError("samples must be a nonempty mapping of sample IDs to FASTA paths")
    for name, fasta in samples.items():
        if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", str(name)) or not isinstance(fasta, str):
            raise ValueError("Use safe sample IDs and string FASTA paths")
    modes = cfg.get("alpha_databases")
    if not isinstance(modes, list) or not modes or len(set(modes)) != len(modes) or not set(modes) <= {"T2T-CHM13_database", "Global_database"}:
        raise ValueError("alpha_databases must be [T2T-CHM13_database], [Global_database], or both")
    for mode in modes:
        if not cfg.get("databases", {}).get(mode):
            raise ValueError(f"Missing FASTA for alpha database {mode}")
    root = cfg.get("output_dir", "")
    if not isinstance(root, str) or not root or "{" in root or "}" in root or root.rstrip("/\\") in {"", ".", ".."}:
        raise ValueError("output_dir must be a distinct named output directory")
    rm = cfg.get("repeatmasker", {})
    if not rm.get("library_label") or not rm.get("expected_version"):
        raise ValueError("Set RepeatMasker expected_version and library_label")
    if not isinstance(rm.get("library_files"), list) or not rm["library_files"]:
        raise ValueError("Set repeatmasker.library_files to the actual installed repeat library/partition files; see README")
    if rm.get("species") != "human":
        raise ValueError("This workflow is for human assemblies: species must be human")
    if int(rm.get("threads", 0)) < 4:
        raise ValueError("RepeatMasker/RMBlast requires at least 4 allocated threads")
    hsat = cfg.get("hsat23", {})
    for key in ("perl_script", "kmer2", "kmer3"):
        if not isinstance(hsat.get(key), str) or not hsat[key].strip():
            raise ValueError(f"hsat23.{key} must be a nonempty file path; see README")
    ref = cfg.get("reference", {})
    if ref.get("enabled"):
        for key in ("genome", "centromere_bed"):
            if not isinstance(ref.get(key), str) or not ref[key]:
                raise ValueError(f"reference.{key} is required when reference.enabled is true")
        if ref.get("mapping_method") != "minimap2":
            raise ValueError("reference.mapping_method must be minimap2")
        for key in ("minimap2_executable", "bedtools_executable"):
            if not isinstance(ref.get(key), str) or not ref[key]:
                raise ValueError(f"reference.{key} must be a nonempty command or path")
        for key in ("flank", "edge_warning", "outside_warning", "min_mapq",
                    "max_unaligned_bp", "threads", "mem_mb"):
            if not isinstance(ref.get(key), int) or ref[key] < 0:
                raise ValueError(f"reference.{key} must be a nonnegative integer")
        if ref["flank"] == 0 or ref["threads"] == 0 or ref["mem_mb"] == 0:
            raise ValueError("reference.flank, threads and mem_mb must be positive")
        coverage = ref.get("min_query_coverage")
        if not isinstance(coverage, (int, float)) or not 0 < coverage <= 1:
            raise ValueError("reference.min_query_coverage must be in (0, 1]")
        if not isinstance(ref.get("preset"), str) or not ref["preset"]:
            raise ValueError("reference.preset must be a nonempty minimap2 preset")
    elif any(ref.get(key) for key in ("genome", "centromere_bed")):
        raise ValueError("Clear reference genome/BED when reference.enabled is false")
    for key in ("initial_flank", "extension_distance", "long_satellite", "step_padding", "stop_padding"):
        value = cfg["boundaries"][key]
        if not isinstance(value, int) or value < 0:
            raise ValueError(f"boundaries.{key} must be a nonnegative integer")
    if cfg["boundaries"]["step_padding"] < 1:
        raise ValueError("step_padding must be positive to ensure extension terminates")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def chrom_key(name):
    return tuple((0, int(x)) if x.isdigit() else (1, x) for x in re.split(r"(\d+)", name))


def canonical_chrom(name, aliases):
    if name in aliases:
        return aliases[name]
    # Exact suffix only: chr22_random or chr220 must not trigger chr22 behavior.
    match = re.search(r"(?:^|[#_.-])(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y))$", name)
    return match.group(1) if match else None


def read_fasta_lengths(path):
    """Read sequence lengths without creating an index beside a user-owned FASTA."""
    lengths, name, count = {}, None, 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    if count == 0:
                        raise ValueError(f"Empty FASTA record {name}")
                    lengths[name] = count
                fields = line[1:].split()
                if not fields or fields[0] in lengths:
                    raise ValueError("Empty or duplicate FASTA identifier")
                name, count = fields[0], 0
            elif line.strip():
                if name is None:
                    raise ValueError("Sequence encountered before FASTA header")
                count += len(line.strip())
    if name is None or count == 0:
        raise ValueError("Empty FASTA / final record")
    lengths[name] = count
    return lengths


def read_fai(path):
    lengths = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 5 or fields[0] in lengths or int(fields[1]) <= 0:
            raise ValueError(f"Invalid/duplicate FAI entry: {line}")
        lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError("Empty FASTA index")
    return lengths


def bed9(chrom, start, end, kind, strand="+", score="0"):
    return [chrom, int(start), int(end), str(kind), str(score), strand,
            int(start), int(end), COLORS.get(str(kind), COLORS["ASat"])]


def read_bed(path, lengths=None):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) not in (3, 6, 9):
                raise ValueError(f"{path}:{number}: expected BED3/6/9")
            row[1:3] = [int(row[1]), int(row[2])]
            if row[1] < 0 or row[1] >= row[2]:
                raise ValueError(f"{path}:{number}: invalid interval")
            if lengths is not None and (row[0] not in lengths or row[2] > lengths[row[0]]):
                raise ValueError(f"{path}:{number}: interval outside input FASTA")
            if len(row) >= 6 and row[5] not in ("+", "-", "."):
                raise ValueError(f"{path}:{number}: invalid strand")
            rows.append(row)
    return sorted(rows, key=lambda r: (chrom_key(r[0]), r[1], r[2], str(r[3:])) )


def write_bed(path, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


def reference_anchors(centromeres, lengths, flank):
    """Create BED6 flanks around reference centromeres for independent mapping."""
    anchors = []
    seen = set()
    for row in centromeres:
        chrom, start, end = row[:3]
        if chrom in seen:
            raise ValueError(f"Reference centromere BED has multiple rows for {chrom}")
        seen.add(chrom)
        if chrom not in lengths or end > lengths[chrom]:
            raise ValueError(f"Reference centromere outside reference FASTA: {chrom}:{start}-{end}")
        left = max(0, start - flank)
        right = min(lengths[chrom], end + flank)
        if left < start:
            anchors.append([chrom, left, start, f"{chrom}|left", "0", "+"])
        if end < right:
            anchors.append([chrom, end, right, f"{chrom}|right", "0", "+"])
    return anchors


def anchor_paf(path, target_lengths, min_mapq, min_query_coverage, max_unaligned_bp):
    """Convert qualifying minimap2 PAF anchor mappings to projected BED6 intervals.

    All passing mappings are retained so downstream code can reject ambiguity. Small
    unaligned query ends are projected 1:1; larger overhangs fail the filter.
    """
    mapped, qc = [], []
    queries_seen = set()
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"{path}:{number}: expected at least 12 PAF columns")
            query, target = fields[0], fields[5]
            qlen, qstart, qend = map(int, fields[1:4])
            strand = fields[4]
            tlen, tstart, tend = map(int, fields[6:9])
            mapq = int(fields[11])
            queries_seen.add(query)
            if target not in target_lengths or tlen != target_lengths[target]:
                continue
            coverage = (qend - qstart) / qlen if qlen else 0
            left_overhang, right_overhang = qstart, qlen - qend
            reasons = []
            if mapq < min_mapq:
                reasons.append("low_mapq")
            if coverage < min_query_coverage:
                reasons.append("low_query_coverage")
            if max(left_overhang, right_overhang) > max_unaligned_bp:
                reasons.append("excess_unaligned_query_end")
            if strand != "+":
                reasons.append("reverse_orientation")
            if reasons:
                qc.append([query, target, mapq, f"{coverage:.6f}", qstart, qend,
                           tstart, tend, ";".join(reasons), "rejected"])
                continue
            projected_start = max(0, tstart - left_overhang)
            projected_end = min(tlen, tend + right_overhang)
            if projected_start >= projected_end:
                qc.append([query, target, mapq, f"{coverage:.6f}", qstart, qend,
                           tstart, tend, "empty_projection", "rejected"])
                continue
            mapped.append([target, projected_start, projected_end, query, str(mapq), strand])
            qc.append([query, target, mapq, f"{coverage:.6f}", qstart, qend,
                       tstart, tend, "pass", "retained"])
    return sorted(mapped, key=lambda r: (chrom_key(r[0]), r[1], r[3])), qc, queries_seen


def projected_windows(lifted, target_lengths, boundary_settings):
    """Infer target search windows from unique lifted reference flanks.

    The target assembly must be chromosome-forward. Non-acrocentric chromosomes
    require both flanks; acrocentric chromosomes require the right flank only.
    """
    aliases = boundary_settings.get("chromosome_aliases", {})
    acrocentrics = set(boundary_settings["acrocentric_chromosomes"])
    grouped = collections.defaultdict(lambda: collections.defaultdict(list))
    for row in lifted:
        if len(row) < 6 or "|" not in row[3]:
            raise ValueError("Lifted anchors must be BED6 with <reference chromosome>|<left/right> names")
        ref_chrom, side = row[3].rsplit("|", 1)
        if side not in ("left", "right"):
            raise ValueError(f"Unexpected anchor side: {side}")
        if row[0] in target_lengths:
            grouped[(row[0], ref_chrom)][side].append(row)

    windows, cores, qc = [], [], []
    for target_chrom, length in target_lengths.items():
        canonical = canonical_chrom(target_chrom, aliases)
        candidates = [(ref, sides) for (target, ref), sides in grouped.items()
                      if target == target_chrom and (canonical is None or canonical_chrom(ref, {}) == canonical)]
        if len(candidates) != 1:
            qc.append([target_chrom, canonical or "unknown", "ambiguous_or_missing_reference_chromosome",
                       ".", ".", ".", ".", "yes"])
            continue
        ref_chrom, sides = candidates[0]
        required = ("right",) if canonical in acrocentrics else ("left", "right")
        if any(len(sides[side]) != 1 for side in required):
            qc.append([target_chrom, canonical or "unknown", "ambiguous_or_missing_lifted_anchor",
                       ".", ".", ".", ".", "yes"])
            continue
        used = [sides[side][0] for side in required]
        if any(row[5] != "+" for row in used):
            qc.append([target_chrom, canonical or "unknown", "reverse_orientation_not_supported",
                       ".", ".", ".", ".", "yes"])
            continue
        if canonical in acrocentrics:
            core_start, core_end = 0, used[0][1]
            start, end = 0, used[0][2]
        else:
            left, right = sides["left"][0], sides["right"][0]
            if left[1] >= right[2]:
                qc.append([target_chrom, canonical or "unknown", "lifted_anchors_out_of_order",
                           ".", ".", ".", ".", "yes"])
                continue
            core_start, core_end = left[2], right[1]
            start, end = left[1], right[2]
        start, end = max(0, start), min(length, end)
        core_start, core_end = max(0, core_start), min(length, core_end)
        if start >= end or core_start >= core_end:
            qc.append([target_chrom, canonical or "unknown", "empty_or_inverted_projected_interval",
                       ".", ".", ".", ".", "yes"])
            continue
        windows.append([target_chrom, start, end, ref_chrom, "0", "+"])
        cores.append([target_chrom, core_start, core_end, ref_chrom, "0", "+"])
        qc.append([target_chrom, canonical or "unknown", "projected_window", core_start, core_end,
                   start, end, "no"])
    return windows, cores, qc


def reference_guided_boundaries(selected, windows, baseline, lengths, boundary_settings,
                                edge_warning, full=None, outside_warning=0):
    """Use sample satellites inside projected windows to define guided boundaries."""
    aliases = boundary_settings.get("chromosome_aliases", {})
    acrocentrics = set(boundary_settings["acrocentric_chromosomes"])
    by_selected = collections.defaultdict(list)
    for row in selected:
        by_selected[row[0]].append(row)
    by_full = collections.defaultdict(list)
    for row in full or []:
        by_full[row[0]].append(row)
    by_window = {row[0]: row for row in windows}
    by_baseline = {row[0]: row for row in baseline}
    output, qc = [], []
    for chrom in lengths:
        canonical = canonical_chrom(chrom, aliases)
        base = by_baseline.get(chrom, [chrom, ".", "."])
        if chrom not in by_window:
            qc.append([chrom, canonical or "unknown", "no_projected_window", ".", ".",
                       base[1], base[2], ".", ".", "yes"])
            continue
        window = by_window[chrom]
        rows = by_selected.get(chrom, [])
        if not rows:
            qc.append([chrom, canonical or "unknown", "no_satellites_in_projected_window",
                       window[1], window[2], base[1], base[2], ".", ".", "yes"])
            continue
        start = 0 if canonical in acrocentrics else min(row[1] for row in rows)
        end = max(row[2] for row in rows)
        warnings = []
        if canonical not in acrocentrics and start - window[1] <= edge_warning:
            warnings.append("left_satellite_near_window_edge")
        if window[2] - end <= edge_warning:
            warnings.append("right_satellite_near_window_edge")
        if canonical not in acrocentrics and any(
                row[2] <= window[1] and window[1] - row[2] <= outside_warning
                for row in by_full[chrom]):
            warnings.append("nearby_satellite_left_of_window")
        if any(row[1] >= window[2] and row[1] - window[2] <= outside_warning
               for row in by_full[chrom]):
            warnings.append("nearby_satellite_right_of_window")
        status = ";".join(warnings) if warnings else "reference_guided_satellite_envelope"
        manual = "yes" if warnings or canonical is None else "no"
        output.append([chrom, start, end])
        qc.append([chrom, canonical or "unknown", status, window[1], window[2],
                   base[1], base[2], start, end, manual])
    return output, qc


def alpha_hits(path, mode, min_identity=90, min_length=100):
    """Keep best hit per query interval, then best direct overlapping hit.

    This is the legacy local-neighborhood selection, NOT a new greedy global
    nonoverlap algorithm. Ties are deterministic by input order/coordinate.
    """
    best = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            f = line.split()
            if len(f) != 12:
                raise ValueError("BLAST output must have 12 columns (outfmt 6)")
            q1, q2, s1, s2 = map(int, f[6:10])
            if min(q1, q2, s1, s2) < 1:
                raise ValueError("BLAST positions must be 1-based positive integers")
            start, end = min(q1, q2) - 1, max(q1, q2)
            key = f[0], start, end
            score = float(f[11])
            if key not in best or score > best[key][0]:
                strand = "+" if (q2 - q1) * (s2 - s1) >= 0 else "-"
                name = re.sub(r"^([0-9]+)_.*", r"\1", f[1]) if mode == "Global_database" else f[1]
                best[key] = score, float(f[2]), int(f[3]), [f[0], start, end, name, f[11], strand]
    by_chrom = collections.defaultdict(list)
    for score, identity, length, row in best.values():
        if identity >= min_identity and length >= min_length:
            by_chrom[row[0]].append(row)
    chosen = set()
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda r: (r[1], r[2], r[3]))
        starts = [r[1] for r in rows]
        # Prefix maximum ends bound the neighborhood even for nested alignments.
        maxends, largest = [], 0
        for row in rows:
            largest = max(largest, row[2])
            maxends.append(largest)
        for row in rows:
            lo = bisect.bisect_right(maxends, row[1])
            hi = bisect.bisect_left(starts, row[2])
            neighbors = (r for r in rows[lo:hi] if r[2] > row[1])
            winner = max(neighbors, key=lambda r: float(r[4]))
            chosen.add(tuple(winner))
    return sorted(chosen, key=lambda r: (chrom_key(r[0]), r[1], r[2], r[3]))


def merge_hsat1(rows, gap=5000):
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r[0], r[3], r[5]].append(r)
    output = []
    for (chrom, kind, strand), group in groups.items():
        current = None
        for row in sorted(group, key=lambda r: r[1]):
            if current is not None and row[1] <= current[2] + gap:
                current[2] = current[7] = max(current[2], row[2])
            else:
                if current is not None:
                    output.append(current)
                current = row.copy()
        if current is not None:
            output.append(current)
    return sorted(output, key=lambda r: (chrom_key(r[0]), r[1], r[2]))


def parse_repeatmasker_gff(path, lengths):
    outputs = {"hsat1": [], "bsat": [], "gsat": []}
    count = 0
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                raise ValueError("Unrecognized RepeatMasker GFF; expected 9 tab-separated fields")
            match = re.search(r'Motif:([^";\s]+)', f[8])
            if not match:
                raise ValueError("RepeatMasker GFF lacks Motif: identifiers; review parser for this version")
            count += 1
            motif = match.group(1)
            start, end = int(f[3]) - 1, int(f[4])
            if f[0] not in lengths or not 0 <= start < end <= lengths[f[0]]:
                raise ValueError("RepeatMasker GFF coordinates disagree with FASTA")
            strand = "-" if f[6] in ("-", "C") else f[6]
            if strand not in ("+", "-", "."):
                raise ValueError("Unrecognized RepeatMasker strand")
            if motif in {"SAR", "HSATI"}:
                kind = "HSat1A" if motif == "SAR" else "HSat1B"
                outputs["hsat1"].append(bed9(f[0], start, end, kind, strand))
            elif re.search(r"BSAT|Beta|LSAU", motif):
                outputs["bsat"].append(bed9(f[0], start, end, "BSat", strand))
            elif re.search(r"GSAT|TAR1", motif):
                outputs["gsat"].append(bed9(f[0], start, end, "GSat", strand))
    outputs["hsat1"] = merge_hsat1(outputs["hsat1"])
    print(f"Parsed {count} RepeatMasker GFF features; satellite counts: " +
          str({k: len(v) for k, v in outputs.items()}))
    return outputs


def search_window(alpha, annotations, chrom_length, acrocentric, settings):
    """Legacy extension-distance OR long-satellite heuristic, with bounded output.

    Distances are to interval starts (left) or ends (right), not inter-repeat gaps.
    Independent end sorting fixes bisect on nonmonotonic ends in nested records.
    """
    flank = settings["initial_flank"]
    left = max(0, min(r[1] for r in alpha) - flank)
    right = min(chrom_length, max(r[2] for r in alpha) + flank)
    distance, long = settings["extension_distance"], settings["long_satellite"]
    step, stop = settings["step_padding"], settings["stop_padding"]
    if acrocentric:
        left = 0
    else:
        ordered = sorted(annotations, key=lambda r: r[1])
        starts = [r[1] for r in ordered]
        while True:
            i = bisect.bisect_left(starts, left) - 1
            if i < 0:
                break
            r = ordered[i]
            if left - r[1] <= distance or r[2] - r[1] > long:
                left = max(0, r[1] - step)
                if left == 0:
                    break
            else:
                left = max(0, left - stop)
                break
    ordered = sorted(annotations, key=lambda r: r[2])
    ends = [r[2] for r in ordered]
    while True:
        i = bisect.bisect_right(ends, right)
        if i >= len(ordered):
            break
        r = ordered[i]
        if r[2] - right <= distance or r[2] - r[1] > long:
            right = min(chrom_length, r[2] + step)
            if right == chrom_length:
                break
        else:
            right = min(chrom_length, right + stop)
            break
    return left, right


def merge_for_display(rows):
    """Preserve legacy start-sorted adjacency/type/strand merging (not groupwise)."""
    ordered = []
    for row in rows:
        r = row.copy()
        if r[3] not in COLORS:
            r[3] = "ASat"
        r[8] = COLORS[r[3]]
        r[4] = "0"
        ordered.append(r)
    ordered.sort(key=lambda r: (chrom_key(r[0]), r[1], r[3]))
    merged = []
    gaps = {"ASat": 10000, "BSat": 5000, "GSat": 5000}
    for r in ordered:
        previous = merged[-1] if merged else None
        if previous and (previous[0], previous[3], previous[5]) == (r[0], r[3], r[5]) and r[1] <= previous[2] + gaps.get(r[3], 0):
            previous[2] = previous[7] = max(previous[2], r[2])
        else:
            merged.append(r)
    return merged


def final_boundaries(alpha, all_rows, lengths, settings):
    by_alpha, by_all = collections.defaultdict(list), collections.defaultdict(list)
    for r in alpha:
        by_alpha[r[0]].append(r)
    for r in all_rows:
        by_all[r[0]].append(r)
    windows, selected, merged, cent, qc = [], [], [], [], []
    for chrom, length in lengths.items():
        canonical = canonical_chrom(chrom, settings.get("chromosome_aliases", {}))
        acro = canonical in settings["acrocentric_chromosomes"]
        if not by_alpha[chrom]:
            qc.append([chrom, canonical or "unknown", "no_alpha_hits", ".", ".", ".", ".", "yes"])
            continue
        left, right = search_window(by_alpha[chrom], by_all[chrom], length, acro, settings)
        windows.append([chrom, left, right])
        rows = [r for r in by_all[chrom] if r[1] < right and r[2] > left]
        local_merged = merge_for_display(rows)
        selected.extend(rows)
        merged.extend(local_merged)
        contained = [r for r in local_merged if r[1] >= left and r[2] <= right]
        if not contained:
            qc.append([chrom, canonical or "unknown", "no_contained_satellites", left, right, ".", ".", "yes"])
            continue
        start = 0 if acro else min(r[1] for r in contained)
        end = max(r[2] for r in contained)
        cent.append([chrom, start, end])
        status = "acrocentric_start_fixed_0" if acro else "automatic_heuristic"
        if canonical is None:
            status += ";unknown_chromosome_alias"
        qc.append([chrom, canonical or "unknown", status, left, right, start, end, "yes"])
    return windows, selected, merged, cent, qc


def plot_annotation(bed, cent, fai, output, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Patch
    lengths = read_fai(fai)
    annotations = read_bed(bed, lengths)
    bounds = {r[0]: r[1:3] for r in read_bed(cent, lengths)}
    by_chrom = collections.defaultdict(list)
    for r in annotations:
        by_chrom[r[0]].append(r)
    with PdfPages(output) as pdf:
        for chrom in sorted(lengths, key=chrom_key):
            fig, ax = plt.subplots(figsize=(12, 2.6))
            rows = by_chrom[chrom]
            if chrom in bounds:
                start, end = bounds[chrom]
                ax.axvspan(start / 1e6, end / 1e6, color="0.93")
                ax.axvline(start / 1e6, color="black", linestyle="--", lw=0.6)
                ax.axvline(end / 1e6, color="black", linestyle="--", lw=0.6)
                margin = max((end - start) * 0.03, 10000)
                ax.set_xlim(max(0, start - margin) / 1e6, min(lengths[chrom], end + margin) / 1e6)
            else:
                ax.set_xlim(0, lengths[chrom] / 1e6)
                ax.text(0.5, 0.5, "No automatic centromere boundary; inspect boundary_qc.tsv", transform=ax.transAxes, ha="center")
            polygons, colors = [], []
            for r in rows:
                a, b = r[1] / 1e6, r[2] / 1e6
                polygons.append([(a, 0.3), (b, 0.3), (b, 0.7), (a, 0.7)])
                colors.append(tuple(int(c) / 255 for c in r[8].split(",")))
            if polygons:
                ax.add_collection(PolyCollection(polygons, facecolors=colors, edgecolors="none", rasterized=True))
            kinds = sorted({r[3] for r in rows})
            if kinds:
                ax.legend(handles=[Patch(color=tuple(int(c) / 255 for c in COLORS[k].split(",")), label=k) for k in kinds], loc="upper center", ncol=min(len(kinds), 7), fontsize=8)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.set_xlabel("Chromosome position (Mb; BED coordinates)")
            ax.set_title(f"{title} | {chrom} | automatic boundary: manual review required")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

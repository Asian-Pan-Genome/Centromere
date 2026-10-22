"""Snakemake rule entry point. External failures propagate; no shell pipelines."""
import contextlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(snakemake.scriptdir)))
import satellite as sat


def portable_config(config):
    """Retain settings while removing host-specific absolute paths from metadata."""
    data = json.loads(json.dumps(config))
    data["samples"] = {key: Path(value).name for key, value in data.get("samples", {}).items()}
    repeatmasker = data.get("repeatmasker", {})
    if repeatmasker.get("executable"):
        repeatmasker["executable"] = Path(repeatmasker["executable"]).name
    repeatmasker["library_files"] = [Path(value).name for value in repeatmasker.get("library_files", [])]
    hsat = data.get("hsat23", {})
    for key in ("perl_script", "kmer2", "kmer3"):
        if hsat.get(key):
            hsat[key] = Path(hsat[key]).name
    reference = data.get("reference", {})
    for key in ("genome", "centromere_bed", "minimap2_executable", "bedtools_executable"):
        if reference.get(key):
            reference[key] = Path(reference[key]).name
    return data


def run(argv, capture=False, cwd=None):
    argv = [str(x) for x in argv]
    print("COMMAND " + json.dumps(argv), flush=True)
    result = subprocess.run(argv, check=True, cwd=cwd,
                            stdout=subprocess.PIPE if capture else sys.stdout,
                            stderr=subprocess.STDOUT, text=True)
    if capture:
        print(result.stdout, flush=True)
    return result.stdout if capture else None


def prepare(s):
    source = Path(s.input.fasta)
    ids, name, count, lengths = set(), None, 0, {}
    with open(source) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    if not count:
                        raise ValueError(f"Empty FASTA record {name}")
                    lengths[name] = count
                tokens = line[1:].split()
                if not tokens or tokens[0] in ids:
                    raise ValueError("Empty or duplicate FASTA identifier")
                name, count = tokens[0], 0
                ids.add(name)
            elif line.strip():
                seq = line.strip()
                if name is None or re.search(r"[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]", seq):
                    raise ValueError("Expected uncompressed nucleotide FASTA with IUPAC bases")
                count += len(seq)
    if name is None or count == 0:
        raise ValueError("Empty FASTA / final record")
    lengths[name] = count
    shutil.copyfile(source, s.output.fasta)
    run(["samtools", "faidx", s.output.fasta])
    if sat.read_fai(s.output.fai) != lengths:
        raise ValueError("FASTA index disagrees with validated sequence lengths")
    sat.write_json(s.output.metadata, {"source_name": source.name,
                                     "sha256": sat.sha256(source), "lengths": lengths})


def provenance(s):
    rm = s.config["repeatmasker"]
    versions = {}
    commands = {"RepeatMasker": [rm["executable"], "-v"],
                "blastn": ["blastn", "-version"], "makeblastdb": ["makeblastdb", "-version"],
                "rmblastn": ["rmblastn", "-version"], "samtools": ["samtools", "--version"],
                "perl": ["perl", "-e", 'printf "%vd", $^V;']}
    for name, command in commands.items():
        versions[name] = run(command, capture=True).strip()
    ref = s.config.get("reference", {})
    if ref.get("enabled"):
        versions["bedtools"] = run([ref["bedtools_executable"], "--version"], capture=True).strip()
        versions["minimap2"] = run([ref["minimap2_executable"], "--version"], capture=True).strip()
    expected = str(rm["expected_version"])
    if not re.search(r"(?<![0-9.])" + re.escape(expected) + r"(?![0-9.])", versions["RepeatMasker"]):
        raise ValueError(f"RepeatMasker version does not match expected_version={expected}; update configuration deliberately")
    for module in ("snakemake", "matplotlib", "PyYAML"):
        versions[module] = importlib.metadata.version(module)
    versions["python"] = sys.version
    hashes = [{"name": Path(p).name, "sha256": sat.sha256(p)} for p in s.input]
    conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
    packages = None
    if conda:
        packages = json.loads(run([conda, "list", "--json", "--prefix", sys.prefix], capture=True))
    sat.write_json(s.output[0], {"config": portable_config(s.config), "versions": versions,
                               "file_sha256": hashes, "conda_packages": packages,
                               "library_label_status": "user-declared; verify against installed library metadata and RepeatMasker log",
                               "reference_projection": "enabled" if s.config.get("reference", {}).get("enabled") else "disabled"})


def reference_anchors(s):
    source = Path(s.input.genome).resolve()
    genome = Path(s.output.genome)
    if genome.exists() or genome.is_symlink():
        genome.unlink()
    genome.symlink_to(source)
    run(["samtools", "faidx", genome])
    lengths = sat.read_fai(s.output.fai)
    centromeres = sat.read_bed(s.input.centromeres, lengths)
    anchors = sat.reference_anchors(centromeres, lengths, s.params.flank)
    sat.write_bed(s.output.anchors, anchors)
    run([s.params.bedtools, "getfasta", "-nameOnly", "-fi", genome,
         "-bed", s.output.anchors, "-fo", s.output.fasta])
    sat.write_json(s.output.metadata, {"coordinate_system": "0-based, half-open",
                                      "flank_bp": s.params.flank,
                                      "reference_genome": source.name,
                                      "reference_centromere_bed": Path(s.input.centromeres).name,
                                      "anchor_count": len(anchors)})


def map_reference_anchors(s):
    argv = [s.params.executable, "-x", s.params.preset, "-c", "--cs=long",
            "-t", str(s.threads), s.input.target, s.input.anchors]
    print("COMMAND " + json.dumps([str(x) for x in argv]), flush=True)
    with open(s.output.paf, "w", encoding="utf-8", newline="\n") as handle:
        subprocess.run([str(x) for x in argv], check=True, stdout=handle,
                       stderr=sys.stdout, text=True)
    mapped, qc, queries = sat.anchor_paf(
        s.output.paf, sat.read_fai(s.input.fai), s.params.min_mapq,
        s.params.min_query_coverage, s.params.max_unaligned_bp)
    sat.write_bed(s.output.mapped, mapped)
    with open(s.output.qc, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("anchor\ttarget\tmapq\tquery_coverage\tquery_start\tquery_end\ttarget_start\ttarget_end\treason\tstatus\n")
        for row in qc:
            handle.write("\t".join(map(str, row)) + "\n")
    print(f"Observed {len(queries)} anchor queries; retained {len(mapped)} mappings")


def projected_windows(s):
    lengths = sat.read_fai(s.input.fai)
    lifted = sat.read_bed(s.input.mapped)
    windows, cores, qc = sat.projected_windows(lifted, lengths, s.params.boundaries)
    sat.write_bed(s.output.windows, windows)
    sat.write_bed(s.output.cores, cores)
    with open(s.output.qc, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("target_chromosome\tcanonical_chromosome\tstatus\tprojected_core_start\tprojected_core_end\twindow_start\twindow_end\tmanual_review_required\n")
        for row in qc:
            handle.write("\t".join(map(str, row)) + "\n")


def reference_intersect(s):
    # -wa/-u retains each complete target satellite record once, even if it
    # crosses a projected window edge. Sorting is restored for deterministic output.
    with tempfile.TemporaryDirectory(prefix="reference_intersect_", dir=Path(s.output[0]).parent) as tmp:
        raw = Path(tmp) / "intersect.bed"
        with open(raw, "w", encoding="utf-8") as handle:
            argv = [s.params.executable, "intersect", "-wa", "-u", "-a", s.input.satellites, "-b", s.input.windows]
            print("COMMAND " + json.dumps([str(x) for x in argv]), flush=True)
            subprocess.run([str(x) for x in argv], check=True, stdout=handle, stderr=sys.stdout, text=True)
        rows = sat.read_bed(raw, sat.read_fai(s.input.fai))
        sat.write_bed(s.output[0], rows)


def reference_boundaries(s):
    lengths = sat.read_fai(s.input.fai)
    selected = sat.read_bed(s.input.selected, lengths)
    windows = sat.read_bed(s.input.windows, lengths)
    baseline = sat.read_bed(s.input.baseline, lengths)
    full = sat.read_bed(s.input.full, lengths)
    cent, qc = sat.reference_guided_boundaries(selected, windows, baseline, lengths,
                                               s.params.boundaries, s.params.edge_warning,
                                               full, s.params.outside_warning)
    sat.write_bed(s.output.cent, cent)
    with open(s.output.qc, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("chromosome\tcanonical_chromosome\tstatus\tprojected_window_start\tprojected_window_end\tbaseline_start\tbaseline_end\treference_guided_start\treference_guided_end\tmanual_review_required\n")
        for row in qc:
            handle.write("\t".join(map(str, row)) + "\n")


def reference_merge(s):
    lengths = sat.read_fai(s.input.fai)
    rows = sat.read_bed(s.input.satellites, lengths)
    merged = sat.merge_for_display(rows)
    sat.write_bed(s.output[0], merged)


def reference_plot(s):
    sat.plot_annotation(s.input.bed, s.input.cent, s.input.fai, s.output[0],
                        f"{s.wildcards.sample} / {s.wildcards.mode} / reference-guided")


def build_database(s):
    prefix = str(Path(s.output.nin).with_suffix(""))
    # Version 4 indexes give a small, explicitly tracked .nin/.nhr/.nsq set.
    run(["makeblastdb", "-in", s.input.fasta, "-dbtype", "nucl", "-blastdb_version", "4", "-out", prefix])


def alpha_blast(s):
    run(["blastn", "-db", str(Path(s.input.nin).with_suffix("")),
         "-query", s.input.fasta, "-out", s.output[0], "-evalue", s.params.evalue,
         "-outfmt", "6", "-num_threads", s.threads, "-task", "megablast"])


def alpha_filter(s):
    settings = s.params.settings
    rows = sat.alpha_hits(s.input.blast, s.wildcards.mode,
                         settings["min_identity"], settings["min_alignment_length"])
    sat.write_bed(s.output[0], rows)
    sat.read_bed(s.output[0], sat.read_fai(s.input.fai))
    print(f"Retained {len(rows)} monomer intervals")


def hsat23(s):
    lengths = sat.read_fai(s.input.fai)
    with tempfile.TemporaryDirectory(prefix="hsat23_", dir=Path(s.output[0]).parent) as tmp:
        work = Path(tmp).resolve()
        for source, name in ((s.input.kmer2, "HSat2_kmers.txt"), (s.input.kmer3, "HSat3_kmers.txt")):
            shutil.copyfile(source, work / name)
        command = ["perl", str(Path(s.input.perl).resolve()), str(Path(s.input.fasta).resolve()), str(work)]
        print("COMMAND " + json.dumps(command), flush=True)
        result = subprocess.run(command, cwd=work, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
        print(result.stdout, flush=True)
        if result.returncode:
            bins = work / "tempfileA.bed"
            # The unmodified upstream v3 script errors when no 1-kb bin is called.
            if (bins.is_file() and bins.stat().st_size == 0 and
                    "Use of uninitialized value $line1" in result.stdout):
                sat.write_bed(s.output[0], [])
                print("No HSat2/3 bins detected; wrote empty BED.", flush=True)
                return
            result.check_returncode()
        raw = work / "input.HSat2and3_Regions.bed"
        rows = sat.read_bed(raw, lengths)
        for r in rows:
            if len(r) != 9 or r[3] not in ("HSat2", "HSat3"):
                raise ValueError("Unexpected HSat2/3 output schema")
            r[8] = sat.COLORS[r[3]]
        sat.write_bed(s.output[0], rows)


def repeatmasker(s):
    if s.threads < 4:
        raise ValueError("Use --cores >= 4 (prefer 16+) for RepeatMasker/RMBlast")
    settings = s.params.settings
    with tempfile.TemporaryDirectory(prefix="repeatmasker_", dir=Path(s.output.gff).parent) as tmp:
        work = Path(tmp).resolve()
        # Isolate all side effects and partial outputs from both source and other jobs.
        fasta = work / "input.fasta"
        shutil.copyfile(s.input.fasta, fasta)
        executable = settings["executable"]
        if os.path.sep in executable:
            executable = str(Path(executable).resolve())
        command = [executable, "-species", settings["species"], "-e", "rmblast",
                   "-s", "-pa", str(s.threads // 4), str(fasta), "-gff", "-dir", str(work)]
        run(command, cwd=work)
        out = work / "input.fasta.out"
        gff = work / "input.fasta.out.gff"
        if not out.is_file():
            raise ValueError("RepeatMasker succeeded but its .out output is missing")
        no_repeats = "no repetitive sequences" in out.read_text().lower()
        if not gff.is_file():
            if not no_repeats:
                raise ValueError("RepeatMasker GFF missing without explicit no-repeats result")
            gff.write_text("##gff-version 2\n# RepeatMasker explicitly reported no repetitive sequences\n")
        shutil.copyfile(out, s.output.out)
        shutil.copyfile(gff, s.output.gff)
        sat.write_json(s.output.metadata, {"executable": Path(command[0]).name,
                                          "species": settings["species"],
                                          "engine": "rmblast", "sensitive": True,
                                          "threads_allocated": s.threads,
                                          "parallel_batches": s.threads // 4,
                                          "explicit_no_repeats": no_repeats,
                                          "library_label_user_declared": settings["library_label"]})


def repeatmasker_satellites(s):
    rows = sat.parse_repeatmasker_gff(s.input.gff, sat.read_fai(s.input.fai))
    for key in ("hsat1", "bsat", "gsat"):
        sat.write_bed(getattr(s.output, key), sorted(rows[key], key=lambda r: (sat.chrom_key(r[0]), r[1], r[2])))


def merge_and_boundaries(s):
    lengths = sat.read_fai(s.input.fai)
    alpha = sat.read_bed(s.input.alpha, lengths)
    all_rows = [sat.bed9(r[0], r[1], r[2], r[3], r[5], r[4]) for r in alpha]
    for source in (s.input.hsat23, s.input.hsat1, s.input.bsat, s.input.gsat):
        all_rows.extend(sat.read_bed(source, lengths))
    all_rows.sort(key=lambda r: (sat.chrom_key(r[0]), r[1], r[2]))
    sat.write_bed(s.output.full, all_rows)
    windows, selected, merged, cent, qc = sat.final_boundaries(alpha, all_rows, lengths, s.params.settings)
    for key, rows in (("window", windows), ("selected", selected), ("merged", merged), ("cent", cent)):
        sat.write_bed(getattr(s.output, key), rows)
    with open(s.output.qc, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("chromosome\tcanonical_chromosome\tstatus\twindow_start\twindow_end\tcentromere_start\tcentromere_end\tmanual_review_required\n")
        for row in qc:
            handle.write("\t".join(map(str, row)) + "\n")


def plot(s):
    sat.plot_annotation(s.input.bed, s.input.cent, s.input.fai, s.output[0],
                        f"{s.wildcards.sample} / {s.wildcards.mode}")


def main(s):
    for path in list(s.output) + list(s.log):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    actions = {f.__name__: f for f in (provenance, prepare, build_database, alpha_blast,
               alpha_filter, hsat23, repeatmasker, repeatmasker_satellites, merge_and_boundaries, plot,
               reference_anchors, map_reference_anchors, projected_windows, reference_intersect,
               reference_boundaries, reference_merge, reference_plot)}
    actions["prepare_fasta"] = prepare
    actions["build_alpha_database"] = build_database
    with open(s.log[0], "w", encoding="utf-8") as log:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            try:
                actions[s.rule](s)
                for output in s.output:
                    if not Path(output).is_file():
                        raise ValueError(f"Required output missing after success: {output}")
                print("SUCCESS", flush=True)
            except Exception:
                import traceback
                traceback.print_exc()
                raise


main(snakemake)

import copy
from pathlib import Path
import sys
import tempfile
import unittest
import runpy
import subprocess
from types import SimpleNamespace
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
import satellite as sat

SETTINGS = {"initial_flank": 5000000, "extension_distance": 2000000,
            "long_satellite": 10000, "step_padding": 100000, "stop_padding": 1900000,
            "acrocentric_chromosomes": ["chr13", "chr14", "chr15", "chr21", "chr22"],
            "chromosome_aliases": {}}


class SatelliteTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def file(self, name, contents):
        p = self.root / name
        p.write_text(contents)
        return p

    def blast(self, subject, start, end, score, sstart=1, send=171):
        return f"chr8\t{subject}\t99\t171\t0\t0\t{start}\t{end}\t{sstart}\t{send}\t1e-30\t{score}\n"

    def test_blast_conversion_and_modes(self):
        path = self.file("a.tsv", self.blast("123_name", 1, 171, 300))
        self.assertEqual(sat.alpha_hits(path, "Global_database")[0][1:4], (0, 171, "123"))
        self.assertEqual(sat.alpha_hits(path, "T2T-CHM13_database")[0][3], "123_name")

    def test_reverse_subject(self):
        path = self.file("a.tsv", self.blast("x", 11, 181, 300, 171, 1))
        self.assertEqual(sat.alpha_hits(path, "T2T-CHM13_database")[0][5], "-")

    def test_query_reverse_coordinates(self):
        path = self.file("a.tsv", self.blast("x", 181, 11, 300))
        row = sat.alpha_hits(path, "T2T-CHM13_database")[0]
        self.assertEqual(row[1:3], (10, 181))
        self.assertEqual(row[5], "-")

    def test_overlap_best_and_bookended(self):
        text = self.blast("a", 1, 171, 200) + self.blast("b", 10, 180, 300) + self.blast("c", 181, 351, 250)
        rows = sat.alpha_hits(self.file("a.tsv", text), "T2T-CHM13_database")
        self.assertEqual([r[3] for r in rows], ["b", "c"])

    def test_low_identity(self):
        path = self.file("a.tsv", self.blast("a", 1, 171, 100).replace("\t99\t", "\t89\t"))
        self.assertEqual(sat.alpha_hits(path, "T2T-CHM13_database"), [])

    def test_empty_inputs(self):
        path = self.file("empty", "")
        self.assertEqual(sat.alpha_hits(path, "T2T-CHM13_database"), [])
        self.assertEqual(sat.read_bed(path), [])
        self.assertEqual(sat.merge_hsat1([]), [])
        self.assertEqual(sat.merge_for_display([]), [])

    def test_gff_and_hsat1_merge(self):
        rows = ["chr8\tRM\tsimilarity\t1\t100\t10\t+\t.\tTarget \"Motif:SAR\" 1 100",
                "chr8\tRM\tsimilarity\t5101\t5200\t10\t+\t.\tTarget \"Motif:SAR\" 1 100",
                "chr8\tRM\tsimilarity\t7001\t7100\t10\t-\t.\tTarget \"Motif:BSAT\" 1 100"]
        out = sat.parse_repeatmasker_gff(self.file("a.gff", "\n".join(rows)), {"chr8": 10000})
        self.assertEqual(out["hsat1"][0][1:3], [0, 5200])
        self.assertEqual(out["bsat"][0][1:3], [7000, 7100])
        self.assertEqual(out["gsat"], [])

    def test_unrecognized_gff_rejected(self):
        path = self.file("a.gff", "chr8\tRM\tx\t1\t100\t0\t+\t.\tunknown=foo\n")
        with self.assertRaises(ValueError):
            sat.parse_repeatmasker_gff(path, {"chr8": 1000})

    def test_invalid_bed_rejected(self):
        with self.assertRaises(ValueError):
            sat.read_bed(self.file("a.bed", "chr8\t-1\t100\n"), {"chr8": 1000})

    def test_acrocentric_final_output_starts_at_zero(self):
        chrom = "sample#mat#chr22"
        alpha = [sat.bed9(chrom, 10000000, 11000000, "123")]
        out = sat.final_boundaries(alpha, alpha, {chrom: 50000000}, SETTINGS)
        self.assertEqual(out[3], [[chrom, 0, 11000000]])
        self.assertEqual(out[4][0][2], "acrocentric_start_fixed_0")

    def test_normal_boundary_not_forced_to_zero(self):
        rows = [sat.bed9("chr8", 10000000, 11000000, "123")]
        self.assertEqual(sat.final_boundaries(rows, rows, {"chr8": 50000000}, SETTINGS)[3],
                         [["chr8", 10000000, 11000000]])

    def test_long_satellite_or_condition(self):
        alpha = [sat.bed9("chr8", 10000000, 11000000, "123")]
        distant = sat.bed9("chr8", 22000000, 22020000, "HSat3")
        left, right = sat.search_window(alpha, alpha + [distant], 50000000, False, SETTINGS)
        self.assertEqual(right, 22120000)

    def test_short_distant_satellite_does_not_extend_to_it(self):
        alpha = [sat.bed9("chr8", 10000000, 11000000, "123")]
        distant = sat.bed9("chr8", 22000000, 22001000, "HSat3")
        self.assertEqual(sat.search_window(alpha, alpha + [distant], 50000000, False, SETTINGS)[1], 17900000)

    def test_no_alpha_does_not_invent_boundary(self):
        out = sat.final_boundaries([], [sat.bed9("chr8", 1, 100, "HSat3")], {"chr8": 1000}, SETTINGS)
        self.assertEqual(out[3], [])
        self.assertEqual(out[4][0][2], "no_alpha_hits")

    def test_aliases_and_no_substring_false_positive(self):
        self.assertEqual(sat.canonical_chrom("NC_000022.11", {"NC_000022.11": "chr22"}), "chr22")
        self.assertIsNone(sat.canonical_chrom("chr220", {}))
        self.assertIsNone(sat.canonical_chrom("chr22_random", {}))

    def test_bounds_clipped_to_chromosome(self):
        alpha = [sat.bed9("chr8", 500, 1000, "123")]
        self.assertEqual(sat.search_window(alpha, alpha, 2000, False, SETTINGS), (0, 2000))

    def test_config_database_selection(self):
        cfg = {"samples": {"a": "input.fa"}, "alpha_databases": ["T2T-CHM13_database", "Global_database"],
               "databases": {"T2T-CHM13_database": "a.fa", "Global_database": "b.fa"}, "output_dir": "results/test",
               "repeatmasker": {"library_label": "Dfam 3.8", "expected_version": "4.1.6",
                                "library_files": ["library.h5"], "threads": 4, "species": "human"},
               "hsat23": {"perl_script": "HSat.pl", "kmer2": "HSat2.txt", "kmer3": "HSat3.txt"},
               "boundaries": SETTINGS, "reference": {"enabled": False}}
        for modes in (["T2T-CHM13_database"], ["Global_database"], ["T2T-CHM13_database", "Global_database"]):
            cfg["alpha_databases"] = modes
            sat.validate_config(cfg)
        cfg["alpha_databases"] = ["T2T-CHM13_database", "T2T-CHM13_database"]
        with self.assertRaises(ValueError):
            sat.validate_config(cfg)
        cfg["alpha_databases"] = ["T2T-CHM13_database"]
        missing_hsat = copy.deepcopy(cfg)
        missing_hsat["hsat23"]["kmer2"] = ""
        with self.assertRaisesRegex(ValueError, "hsat23.kmer2"):
            sat.validate_config(missing_hsat)
        cfg["reference"] = {"enabled": True}
        with self.assertRaises(ValueError):
            sat.validate_config(cfg)
        cfg["reference"] = {
            "enabled": True, "mapping_method": "minimap2", "genome": "ref.fa",
            "centromere_bed": "ref.bed", "minimap2_executable": "minimap2",
            "bedtools_executable": "bedtools", "preset": "asm5", "threads": 4,
            "mem_mb": 8000, "flank": 500000, "min_mapq": 20,
            "min_query_coverage": 0.8, "max_unaligned_bp": 10000,
            "edge_warning": 100000, "outside_warning": 500000}
        sat.validate_config(cfg)

    def test_dispatch_merge_two_modes_share_nonalpha(self):
        class NamedPaths(list):
            def __init__(self, **paths):
                super().__init__(str(p) for p in paths.values())
                self.__dict__.update({k: str(v) for k, v in paths.items()})

        fai = self.file("input.fasta.fai", "chr8\t50000000\t6\t60\t61\nchr22\t50000000\t7\t60\t61\n")
        empty = self.file("empty.bed", "")
        shared = self.file("hsat.bed", "chr22\t22000000\t22020000\tHSat3\t0\t+\t22000000\t22020000\t119,160,186\n")
        scriptdir = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
        outputs = {}
        for mode, end in (("T2T-CHM13_database", 11000000), ("Global_database", 12000000)):
            alpha = self.file(f"{mode}.bed", f"chr8\t10000000\t{end}\t123\t300\t+\nchr22\t10000000\t{end}\t123\t300\t+\n")
            out = NamedPaths(**{k: self.root / mode / f"{k}.bed" for k in ("full", "window", "selected", "merged", "cent", "qc")})
            mock = SimpleNamespace(scriptdir=str(scriptdir), rule="merge_and_boundaries",
                                   input=NamedPaths(alpha=alpha, fai=fai, hsat23=shared,
                                                    hsat1=empty, bsat=empty, gsat=empty),
                                   output=out, params=SimpleNamespace(settings=SETTINGS),
                                   log=[str(self.root / mode / "run.log")])
            runpy.run_path(str(scriptdir / "dispatch.py"), init_globals={"snakemake": mock})
            outputs[mode] = sat.read_bed(out.cent)
            self.assertEqual(outputs[mode][1], ["chr22", 0, 22020000])
            self.assertIn("SUCCESS", Path(mock.log[0]).read_text())
        self.assertEqual(outputs["T2T-CHM13_database"][0], ["chr8", 10000000, 11000000])
        self.assertEqual(outputs["Global_database"][0], ["chr8", 10000000, 12000000])

    def test_upstream_hsat_no_bins_is_empty_bed(self):
        class NamedPaths(list):
            def __init__(self, **paths):
                super().__init__(str(p) for p in paths.values())
                self.__dict__.update({k: str(v) for k, v in paths.items()})

        scriptdir = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
        input_paths = NamedPaths(
            fasta=self.file("input.fasta", ">chr8\nACGT\n"),
            fai=self.file("input.fasta.fai", "chr8\t4\t6\t4\t5\n"),
            perl=self.file("Assembly_HSat2and3_v3.pl", ""),
            kmer2=self.file("HSat2_kmers.txt", ""),
            kmer3=self.file("HSat3_kmers.txt", ""))
        output = self.root / "shared" / "HSat23.bed"
        job = SimpleNamespace(scriptdir=str(scriptdir), rule="hsat23", input=input_paths,
                              output=[str(output)], log=[str(self.root / "hsat23.log")])

        def no_bins(command, cwd, **kwargs):
            (Path(cwd) / "tempfileA.bed").touch()
            return subprocess.CompletedProcess(command, 255,
                "Use of uninitialized value $line1 in scalar chomp at Assembly_HSat2and3_v3.pl line 138.\n")

        with mock.patch("subprocess.run", side_effect=no_bins):
            runpy.run_path(str(scriptdir / "dispatch.py"), init_globals={"snakemake": job})
        self.assertEqual(output.read_text(), "")
        self.assertIn("SUCCESS", Path(job.log[0]).read_text())

    def test_reference_merge_produces_merged_bed(self):
        class NamedPaths(list):
            def __init__(self, **paths):
                super().__init__(str(p) for p in paths.values())
                self.__dict__.update({k: str(v) for k, v in paths.items()})

        fai = self.file("input.fasta.fai", "chr2\t50000000\t6\t60\t61\n")
        satellites = self.file("ref_window.bed",
            "chr2\t1000\t2000\tASat\t0\t+\t1000\t2000\t137,22,64\n"
            "chr2\t2000\t3000\tASat\t0\t+\t2000\t3000\t137,22,64\n"
            "chr2\t5000\t6000\tBSat\t0\t+\t5000\t6000\t220,174,208\n")
        out = self.root / "satellites.reference_guided.merged.bed"
        scriptdir = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
        mock = SimpleNamespace(scriptdir=str(scriptdir), rule="reference_merge",
                               input=NamedPaths(fai=fai, satellites=satellites),
                               output=[str(out)],
                               log=[str(self.root / "reference_merge.log")])
        runpy.run_path(str(scriptdir / "dispatch.py"), init_globals={"snakemake": mock})
        rows = sat.read_bed(out, {"chr2": 50000000})
        self.assertEqual([r[1:3] for r in rows], [[1000, 3000], [5000, 6000]])
        self.assertIn("SUCCESS", Path(mock.log[0]).read_text())

    def test_reference_anchors_are_separate_flanks(self):
        rows = [["chr2", 1000, 2000]]
        self.assertEqual(sat.reference_anchors(rows, {"chr2": 5000}, 500), [
            ["chr2", 500, 1000, "chr2|left", "0", "+"],
            ["chr2", 2000, 2500, "chr2|right", "0", "+"]])

    def test_anchor_paf_filters_and_projects_small_overhangs(self):
        paf = self.file("anchors.paf", "\n".join([
            "chr2|left\t500\t5\t495\t+\tsample_hap_chr2\t5000\t100\t590\t485\t490\t60\ttp:A:P",
            "chr2|right\t500\t0\t500\t+\tsample_hap_chr2\t5000\t2400\t2900\t490\t500\t10\ttp:A:P",
            "chr2|right\t500\t0\t500\t+\tsample_hap_chr2\t5000\t2400\t2900\t490\t500\t60\ttp:A:P"]) + "\n")
        mapped, qc, queries = sat.anchor_paf(
            paf, {"sample_hap_chr2": 5000}, min_mapq=20,
            min_query_coverage=0.8, max_unaligned_bp=10)
        self.assertEqual(queries, {"chr2|left", "chr2|right"})
        self.assertEqual(mapped[0][:4], ["sample_hap_chr2", 95, 595, "chr2|left"])
        self.assertEqual(mapped[1][:4], ["sample_hap_chr2", 2400, 2900, "chr2|right"])
        self.assertTrue(any(row[-2] == "low_mapq" for row in qc))

    def test_projected_windows_nonacrocentric_and_sample_name(self):
        lifted = [
            ["sample_hap_chr2", 100, 600, "chr2|left", "0", "+"],
            ["sample_hap_chr2", 2400, 2900, "chr2|right", "0", "+"]]
        windows, cores, qc = sat.projected_windows(lifted, {"sample_hap_chr2": 4000}, SETTINGS)
        self.assertEqual(windows, [["sample_hap_chr2", 100, 2900, "chr2", "0", "+"]])
        self.assertEqual(cores, [["sample_hap_chr2", 600, 2400, "chr2", "0", "+"]])
        self.assertEqual(qc[0][-1], "no")

    def test_projected_windows_acrocentric_uses_right_anchor(self):
        lifted = [["sample#hap#chr22", 2400, 2900, "chr22|right", "0", "+"]]
        windows, cores, qc = sat.projected_windows(lifted, {"sample#hap#chr22": 4000}, SETTINGS)
        self.assertEqual(windows[0][:3], ["sample#hap#chr22", 0, 2900])
        self.assertEqual(cores[0][:3], ["sample#hap#chr22", 0, 2400])
        self.assertEqual(qc[0][2], "projected_window")

    def test_reference_guided_boundary_and_edge_qc(self):
        lengths = {"sample_hap_chr2": 5000}
        windows = [["sample_hap_chr2", 500, 3500, "chr2", "0", "+"]]
        selected = [sat.bed9("sample_hap_chr2", 550, 3400, "ASat")]
        full = selected + [sat.bed9("sample_hap_chr2", 3600, 3700, "BSat")]
        cent, qc = sat.reference_guided_boundaries(
            selected, windows, [["sample_hap_chr2", 700, 3300]], lengths,
            SETTINGS, edge_warning=100, full=full, outside_warning=500)
        self.assertEqual(cent, [["sample_hap_chr2", 550, 3400]])
        self.assertIn("left_satellite_near_window_edge", qc[0][2])
        self.assertIn("right_satellite_near_window_edge", qc[0][2])
        self.assertIn("nearby_satellite_right_of_window", qc[0][2])
        self.assertEqual(qc[0][-1], "yes")


if __name__ == "__main__":
    unittest.main()

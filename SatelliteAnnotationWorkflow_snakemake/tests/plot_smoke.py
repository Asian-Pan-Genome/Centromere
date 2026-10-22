"""Optional local plotting smoke test; this is NOT a sequence annotation test."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
import satellite as sat

parser = argparse.ArgumentParser()
parser.add_argument("--output-dir", required=True)
args = parser.parse_args()
root = Path(args.output_dir)
root.mkdir(parents=True, exist_ok=True)
(root / "input.fai").write_text("chr8\t50000000\t6\t60\t61\nchr22\t50000000\t7\t60\t61\nchr3\t50000000\t7\t60\t61\n")
sat.write_bed(root / "merged.bed", [sat.bed9("chr8", 10000000, 11000000, "ASat"),
                                    sat.bed9("chr22", 10000000, 12000000, "ASat"),
                                    sat.bed9("chr22", 22000000, 22020000, "HSat3")])
sat.write_bed(root / "cent.bed", [["chr8", 10000000, 11000000], ["chr22", 0, 22020000]])
sat.plot_annotation(root / "merged.bed", root / "cent.bed", root / "input.fai",
                    root / "smoke.pdf", "Synthetic plotting check (not real sequence data)")
print(root / "smoke.pdf")

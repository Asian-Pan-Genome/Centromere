# APG Snakemake Bundle (submission copy)

This folder packages the core Snakemake workflow for the paper. Everything is self-contained except for the actual input data, HORSCAN binary, and conda environments.

## Layout
- `workflow/Snakefile` – entry point; points to `../config/config.yaml`.
- `workflow/rules/*.smk` – formatting, HORSCAN, distance matrix, mutation steps.
- `workflow/envs/formatting.yaml` – placeholder (empty in source project).
- `config/config.yaml` – template with relative paths; edit for real runs.
- `data/sample_sheet.example.tsv` – column schema for the sample sheet.
- `scripts/` – Python scripts called by rules (formatting, horscan, distance, plotting, mutation).

## How to run (after filling paths)
1) Put real inputs in a data directory and update `config/config.yaml` paths (sample_sheet, QC, blast outputs, HORSCAN binary, fasta dirs, etc.).
2) From this `summary/` directory:
```
snakemake -s workflow/Snakefile -j 8 --use-conda
```
Adjust `-j` for available cores; add `-n` to dry-run. DAG: `snakemake -s workflow/Snakefile --dag | dot -Tpdf > dag.pdf`.

## Required inputs
- Sample sheet columns: `sample, name, hap, mon_bed, merged_bed, graph_dir, HiCAT_dir, assembly_fa` (see `data/sample_sheet.example.tsv`).
- HORSCAN v2 binary path (`tools.horscan_v2`) and mode string.
- blast outputs directory (`blast_output_dir`) and optionally `blast_update_dir`.
- QC file path (`sample_QC`) and `hor_summary` template.

## Notes & gaps
- Two scripts referenced in `workflow/rules/mutation.smk` are missing from the source repo: `scripts/04_mutation/plot_synteny_sample_all.py` and `scripts/04_mutation/phylogenetic_monomer_validation.py`. Add them or disable the corresponding rules before execution.
- `workflow/envs/formatting.yaml` is empty; supply a conda env if you plan to use `--use-conda`.

## Packing for submission
From repo root:
```
tar -czf APG_snakemake_summary.tar.gz summary/
```
This creates a lightweight archive without large data files.

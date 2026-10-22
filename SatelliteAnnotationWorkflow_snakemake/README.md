# SatelliteAnnotationWorkflow

A Snakemake workflow for annotating satellite sequences and estimating candidate
centromere boundaries in human chromosome assemblies.

The workflow annotates:

- α-satellite using the T2T-CHM13 and/or Global monomer database;

- HSat2 and HSat3 using diagnostic k-mers;

- HSat1, β-satellite and γ-satellite using RepeatMasker;

- optional reference-guided centromere boundaries using minimap2.

The two α-satellite databases are processed independently and their outputs are not
combined. Automatically inferred boundaries should be manually reviewed.

## Requirements

The workflow runs on Linux. Create the core environment with:

```bash
mamba env create -f environment.yml
conda activate SatelliteAnnotationEnv
```

RepeatMasker and its repeat library must be installed and configured separately.
The workflow has been tested with RepeatMasker 4.1.6 and Dfam 3.8.
The HSat2/3 Perl script and its two k-mer files must also be obtained separately
from <https://github.com/altemose/chm13_hsat> and configured below; they are not
included in this repository.

Verify the required commands:

```bash
command -v RepeatMasker
command -v rmblastn
command -v minimap2
command -v bedtools
```

## Configuration

Copy the example configuration and edit the local copy:

```bash
cp config/config.yaml config/local.yaml
```

### Input sequences

Provide one or more uncompressed chromosome FASTA files:

```yaml
samples:
  sample_Mat_chr2: '/path/to/sample#Mat#chr2.fasta'
```

Sample mapping keys may contain letters, numbers, `_`, `.` and `-`.

### α-satellite databases

Run either database or both:

```yaml
alpha_databases: [T2T-CHM13_database, Global_database]
databases:
  T2T-CHM13_database: resources/MonomerT2TDB/MonomerT2T.consensus.fa
  Global_database: resources/MonomerGlobalDB/MonomerGlobalDB.fasta
```

### RepeatMasker

Specify the RepeatMasker executable and the files from the configured repeat library:

```yaml
repeatmasker:
  executable: RepeatMasker
  expected_version: '4.1.6'
  library_label: 'Dfam 3.8'
  library_files:
    - /path/to/RepeatMasker/Libraries/famdb/dfam38_full.0.h5
    - /path/to/RepeatMasker/Libraries/famdb/rmlib.config
  species: human
  threads: 16
  mem_mb: 32000
```

The listed library files are used for provenance tracking. RepeatMasker must already
be configured to use the intended database.

### HSat2/3 resources

Set paths to your local copies of all three upstream files:

```yaml
hsat23:
  perl_script: /path/to/chm13_hsat/Assembly_HSat2and3_v3.pl
  kmer2: /path/to/chm13_hsat/HSat2_kmers.txt
  kmer3: /path/to/chm13_hsat/HSat3_kmers.txt
```

Absolute paths or paths relative to the workflow root are accepted.&#x20;

### Reference-guided boundaries

Reference guidance is optional. Supply a reference genome and a BED3 file containing
one centromere interval per reference chromosome:

```yaml
reference:
  enabled: true
  mapping_method: minimap2
  genome: /path/to/reference.fasta
  centromere_bed: /path/to/reference.centromere.bed
  minimap2_executable: minimap2
  bedtools_executable: bedtools
  preset: asm5
  threads: 4
  mem_mb: 8000
  flank: 500000
  min_mapq: 20
  min_query_coverage: 0.80
  max_unaligned_bp: 10000
  edge_warning: 100000
  outside_warning: 500000
```

The reference FASTA identifier and BED chromosome field must match exactly. Reference
coordinates define a search window; final guided boundaries are determined from the
satellite annotations in the target assembly.

## Run

Inspect the planned jobs first:

```bash
snakemake --configfile config/local.yaml --cores 16 --dry-run
```

Then run on a compute node:

```bash
snakemake --configfile config/local.yaml \
  --cores 16 --resources mem_mb=64000 --rerun-incomplete
```

A generic single-allocation SLURM example is provided in `run_slurm.sh`.

## Tests

Run the unit tests with:

```bash
python -m unittest discover -s tests -v
```

The C001 chr2 test inputs are included in `tests/data/` (the two FASTA files are
gzip-compressed). Edit the RepeatMasker executable and library paths and all three
HSat2/3 paths in
`tests/data/test.config.yaml` for your installation, or copy that file to
`config/test.local.yaml` and edit the copy. From the workflow root, run:

```bash
bash tests/run_test.sh
```

The script checks `tests/data/md5.txt`, decompresses the FASTA files locally, and uses
`--configfile tests/data/test.config.yaml` by default. To use a local copy, run
`bash tests/run_test.sh config/test.local.yaml`. Direct Snakemake execution with
`--configfile tests/data/test.config.yaml` is possible only after decompression.
See `tests/output/` for illustrative C001 outputs.

## Outputs

For each sample and α-satellite database, the main outputs are:

```text
alpha.monomers.bed
satellites.full.bed
satellites.centromeric.bed
satellites.merged.bed
centromere.bed
boundary_qc.tsv
annotation.pdf
```

When reference guidance is enabled, additional outputs include:

```text
reference_anchors.paf
reference_anchor_mapping_qc.tsv
reference_projected_core.bed
reference_search_window.bed
satellites.reference_window.bed
centromere.reference_guided.bed
reference_boundary_qc.tsv
annotation.reference_guided.pdf
```

The satellite-only and reference-guided boundaries are retained separately.

## Citation and third-party resources

The workflow is derived from the centromere analysis repository:

<https://github.com/Asian-Pan-Genome/Centromere>

See `THIRD_PARTY_NOTICES.md` for external software and data resources.

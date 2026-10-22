# C001 chr2 test inputs

This directory contains the three inputs for the end-to-end test:

- `C001-CHA-E01#Mat#chr2.fa.gz`: compressed target chromosome assembly;
- `CHM13v2.chr2.fa.gz`: compressed chr2 reference sequence (FASTA header `chr2`);
- `CHM13v2.chr2.cent.bed`: one CHM13v2 chr2 centromere interval.

The chr2 reference sequence was extracted from the T2T-CHM13v2.0 analysis-set
FASTA, `chm13v2.0.fa.gz`:

https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

The T2T Consortium describes the analysis set and its CC0 data-release terms at
https://github.com/marbl/CHM13 . The C001 sequence is an internally generated
test assembly. The one-row centromere BED was selected from the supplied CHM13v2
centromere annotation; it is not part of the linked FASTA.

The BED file uses 0-based, half-open coordinates and LF line endings. Each gzip
file is under GitHub's 100 MiB regular-file limit and can be committed directly;
Git LFS and a separate download script are unnecessary. The uncompressed FASTA
files are generated locally by `tests/run_test.sh` and are ignored by Git.

On Linux, verify the three inputs from this directory with:

```bash
md5sum --check md5.txt
```

The test configuration is `test.config.yaml`. Its relative paths assume the
workflow root as the working directory and refer to the locally decompressed
FASTA files, because the workflow does not read gzip FASTA directly. Replace the
RepeatMasker executable and library paths with those installed on your system,
and set the three HSat2/3 paths to files obtained from
https://github.com/altemose/chm13_hsat ,
then run `bash tests/run_test.sh` from the workflow root. Do not invoke Snakemake
directly on a fresh clone before the FASTA files have been decompressed.

Before public distribution, confirm that the C001 target sequence may be shared
under the study's consent and data-use terms.

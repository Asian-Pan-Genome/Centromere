# Third-party software and data

This workflow invokes or references material from third parties. Before a public
release, verify and document the applicable license and citation for each item.

- RepeatMasker and RMBlast: installed separately by the user. RepeatMasker and Dfam
  database files are not redistributed by this workflow.

- Dfam: the validated analysis used Dfam 3.8. Users must obtain and configure the
  database separately under its applicable terms.

- Snakemake, BLAST+, SAMtools, BEDTools, minimap2, Perl, matplotlib and PyYAML:
  installed from their respective Conda channels as specified in `environment.yml`.

- HSat2/3 Perl script and two diagnostic k-mer files: users obtain and configure
  these separately from Nicolas Altemose's upstream repository,
  <https://github.com/altemose/chm13_hsat> .&#x20;

- `resources/MonomerT2TDB` and `resources/MonomerGlobalDB`: record the exact source,
  version, citation and redistribution terms before release.

- T2T-CHM13v2.0 chr2 test reference: extracted from the analysis-set FASTA at
  <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz> .
  The T2T Consortium states that its data are released under CC0 at
  <https://github.com/marbl/CHM13#data-reuse-and-license> .

- C001 test assembly: generated within this study; confirm consent and data-use
  permission before public redistribution.

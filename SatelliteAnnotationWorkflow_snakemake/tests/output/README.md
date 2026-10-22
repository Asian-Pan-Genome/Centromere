# C001 chr2 illustrative outputs

These files are selected outputs from the completed C001 chr2 validation run. The
`Global_database/` and `T2T-CHM13_database/` directories contain separate
alpha-satellite annotations and centromere-boundary results. In particular,
`satellites.full.bed` contains the alpha-monomer calls along with the other
satellite classes; a separate `alpha.monomers.bed` is not needed for this compact
example set.

Reference projections are independent of alpha-satellite database choice and are
produced once per sample under `shared/` in a new run:

```text
reference_projected_core.bed
reference_search_window.bed
reference_projection_qc.tsv
```

The projected core is the reference centromere interval mapped to the sample;
the search window includes the projected flanks. These are intermediate guide
coordinates, not the final personalized centromere BED. Final reference-guided
boundaries are recorded separately for each alpha database in
`centromere.reference_guided.bed`.

These outputs are illustrative files, not proof that every chromosome
boundary can be accepted without manual review. Results can vary with the
RepeatMasker/Dfam and other software versions used.

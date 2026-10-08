# Example Cell Painting screen

Two tab-separated tables that share the `WellId` key; locate them with
`cellpacman::exampleDataPath("features")` and `exampleDataPath("metadata")`.

| File | Rows | Columns | Content |
| --- | --- | --- | --- |
| `example_features.tsv` | 1,920 wells | `WellId` + 664 numeric features | Well-level morphological features ("NonBorder Cells - ..." columns of a Harmony export). |
| `example_metadata.tsv` | 1,920 wells | `WellId`, `Row`, `Column`, ..., `Compound`, `Concentration`, ..., `Plate` | Plate annotation. |

Design: five 384-well plates (`Plate1`–`Plate5`). Each plate holds three
annotated compounds — `Axit` (axitinib), `Cabo` (cabozantinib), `Dacti`
(dactinomycin) — at eight concentrations in duplicate, 16 wells of the `DMSO`
vehicle control, and 320 unannotated wells (`Compound` = `NA`) that the
workflow treats as unknown compounds.

Both tables were derived from one combined instrument export by splitting on
the column names (see the "Preparing your own input" section of
`vignette("getting-started")`).

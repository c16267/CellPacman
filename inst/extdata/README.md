# Synthetic Cell Painting tutorial screen

These two TSV tables contain entirely synthetic data for learning CellPacman.
They are not experimental measurements and do not represent real compound effects.

- `example_features_mock.tsv`: 2,000 wells, `WellId` + 600 numeric feature columns.
- `example_metadata_mock.tsv`: 2,000 wells, 20 annotation columns.
- The shared `WellId` is unique in each table; column names retain the input schema.
- Eight fictitious plates with 250 occupied wells each contain three fictitious references
  (`RefA`, `RefB`, `RefC`) at eight concentrations in duplicate, 16 synthetic
  `DMSO` controls, and 186 unknown wells (`Compound` blank; read as `NA`).
- Reference concentrations (arbitrary tutorial units): 0.012, 0.036, 0.108, 0.324, 0.972, 2.916, 8.748, 26.244.
- No original experimental measurements, experiment IDs, or internal server
  addresses are included. Instrument-style feature names are schema labels only.

After including these files in `inst/extdata/` and reinstalling the package, use
`cellpacman::exampleDataPath("features")` and
`cellpacman::exampleDataPath("metadata")`.

Use `RefA`, `RefB`, and `RefC` in examples requiring named references.

Regenerate these exact tables with `python data-raw/generate_toy_data_mock.py`
from the source tree (Python 3 and NumPy required). No original measurements
are needed. The companion Shiny app loads these same files by default.

# Current revision code

This directory contains the analysis modules assembled for the current
six-main-figure manuscript candidate. It is the primary code tier for author
review. The top-level `scripts/` directory is retained separately as a
historical snapshot.

`nqo1_prioritization/` and `protein_support/` now support supplementary
material rather than a seventh main figure.

Start with:

- `../docs/FIGURE_CODE_MAP.md` for the figure-to-script crosswalk;
- `../docs/DATA_MANIFEST.md` for external inputs and configuration variables;
- `../workflow/run_pipeline.py --list` for the public stage entry points; and
- each module's README for analysis-specific assumptions and limitations.

Generated figures, large matrices, serialized objects, whole-slide images, and
controlled clinical rows are not included. Selected small gene definitions,
locked formula resources, and audit tables are retained only where they are
needed to interpret the analysis. Their redistribution must be confirmed by
the authors before release.

The retained modules do not imply that all historical findings were confirmed.
See `../docs/REPRODUCIBILITY_STATUS.md` before reusing result language.

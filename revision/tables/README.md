# Manuscript support tables

`Table2_frozen_OSARS_128_genes.csv` is tracked because it is a small,
interpretation-critical definition without patient-level rows.

After all prerequisite module outputs have been regenerated, run:

```bash
python revision/scripts/build_summary_tables.py
```

This also assembles aggregate cohort and retrospective-performance tables.
Review the generated files for source-data terms and disclosure risk before
adding them to a public release; they are not automatically approved merely
because the builder completed.

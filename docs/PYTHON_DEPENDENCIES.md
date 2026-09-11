# Python dependencies

Use module-specific environments because the frozen single-cell and model
validation runs used different NumPy/SciPy generations.

## Frozen model-validation environment

Install `environment/requirements-model-validation.txt`. The key recorded
versions are Python 3.12-compatible NumPy 2.5.2, pandas 2.3.3, SciPy 1.18.1,
lifelines 0.30.3, scikit-learn 1.9.0, and scikit-survival 0.28.0. The full
transitive snapshot remains in
`revision/model_validation/references/python_requirements_frozen.txt`.

## Independent single-cell environment

Install `environment/requirements-single-cell.txt`. The analysis manifest
records scanpy 1.11.5, anndata 0.12.19, NumPy 2.4.6, SciPy 1.17.1, and
harmonypy 2.0.0.

## Plotting and supporting analyses

Install `environment/requirements-figures.txt`. It covers the pure plotting,
table, PDF-QC, and standard statistical utilities used across immune, protein,
spatial, pathology, and reviewer panels. GPU cell2location dependencies should
be installed in their own environment following the package documentation;
an exact lock for the historical GPU run was not recovered.

All version files describe the archived reference runs. Test installation in a
clean environment before public release; package versions dated after the
original analysis may change defaults or serialized-model compatibility.

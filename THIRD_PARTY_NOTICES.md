# Third-party notices

This repository contains two source-level third-party components. Their
licenses apply to the identified files independently of the repository's
eventual project-level license.

## Mime1

`revision/model_selection/references/Mime_ML.Dev.Prog.Sig_upstream_20250923.R`
is an unmodified snapshot of `R/ML.Dev.Prog.Sig.R` from
[l-magnificence/Mime](https://github.com/l-magnificence/Mime), commit
`9a9f6ac89851bf631f9df3868b2fa624bed49df2` (snapshot SHA-256
`362b817e83360813c31716aab52c9eead4175fca5d439994ab70560ee2f38737`).
Mime1 declares `MIT + file LICENSE`; copyright 2023 Hongwei Liu, Wei Zhang,
and Yihao Zhang. The license text is in `licenses/Mime1-MIT.txt`.

The HCC-OS wrapper loads this snapshot and applies documented in-memory edits
before evaluation: removal of unused plotting imports, capture of first-stage
features, and a bounded SuperPC retry loop. The upstream snapshot itself is
not edited.

## GeneNMF

`revision/single_cell_review/scripts/07_refit_genenmf_and_select_nmp.R`
contains adapted core functions from [GeneNMF v0.6.2](https://github.com/carmonalab/GeneNMF/tree/v0.6.2),
commit `d74e0254002d6aa03ac94b66c374615e71ee422a`. GeneNMF is licensed under
GPL-3.0-only and is authored by Massimo Andreatta and Santiago Carmona.
The corresponding license is in `licenses/GeneNMF-GPL-3.0.txt`; the adapted
script is marked with its SPDX identifier.

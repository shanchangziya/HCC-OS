# Pre-release checklist

## Completed in this local candidate

- [x] Clone the current public HCC-OS repository into an isolated local candidate.
- [x] Preserve the existing history and remote URL without committing or pushing.
- [x] Add the current revision code while retaining historical scripts for provenance.
- [x] Replace known personal workstation/server paths with portable configuration.
- [x] Exclude virtual environments, caches, raw/large data, models, logs, and figures.
- [x] Add a code/figure crosswalk, data manifest, dependency notes, and stage entry point.
- [x] Identify the Mime1 MIT snapshot and GeneNMF GPL-3 adapted code.
- [x] Compare the filtered revision source bundle with the candidate and document every move, replacement, de-duplication, and review-only exclusion.
- [x] Confirm that tracked CSV/JSON/TXT resources contain definitions, aggregate results, formulas, or hashes rather than patient-level rows.
- [x] Pass the repository hygiene audit, personal-path/secret scan, and final >5 MB file scan.
- [x] Parse 85 R files, compile 58 Python files, validate 2 shell scripts, parse the citation/requirements files, and pass `git diff --check`.
- [x] Generate `docs/CODE_INVENTORY.csv` and `checksums/code_sha256.txt` for all 145 code files.

## Author decisions required before upload

- [x] Update the crosswalk to six main figures and move the former Figure 7 NQO1 analysis to the supplementary section.
- [ ] Confirm that the six-main-figure map matches the final rewritten manuscript.
- [ ] Confirm final supplementary numbering and which optional sensitivity modules remain public.
- [ ] Decide whether the historical `scripts/08_functional_analysis/` archive should remain visible, move to a tagged legacy release, or be removed from the release branch.
- [ ] Approve every tracked derived CSV/JSON/TXT file for redistribution under its source-data terms.
- [ ] Choose a project-level license compatible with the GPL-3.0-only adapted GeneNMF script, or replace that script with a dependency call and re-audit licensing.
- [ ] Verify author order, spelling, affiliations/ORCIDs, manuscript title, journal citation, and DOI in `CITATION.cff`.
- [ ] Decide whether to publish frozen model objects separately (for example, a DOI-bearing archive) and document their checksums/access terms.

## Verification required before upload

- [ ] Run all retained primary stages in clean R/Python environments from approved inputs.
- [ ] Compare regenerated numerical outputs and figure-source tables against frozen hashes/tolerances.
- [ ] Review `git diff` and the full untracked-file list manually.
- [ ] Replace the pre-release `LICENSE` notice only after the license decision.
- [ ] Commit and push only after explicit author approval.

## Upload boundary

Preparation of this directory does not authorize a commit, force-push, release,
or GitHub upload. Those actions are intentionally left for a separate confirmed step.

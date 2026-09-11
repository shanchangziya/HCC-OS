# NQO1 supplementary analysis

This module produces the NQO1 supplementary figure. It provides a transparent
post-model reason for prioritizing NQO1 and does not refit OSARS or claim a
causal mechanism.

## Prioritization rule

1. Intersect the frozen 128 OSARS features with the original
   HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY set. The overlap is `NQO1`,
   `TXNRD1`, and `PRDX1`.
2. Compare their descriptive marker effects and detection-rate gaps in the
   historical discovery OS-high versus OS-low screen.
3. Assess NQO1 protein abundance in the independent PDC000198 HCC proteome.

NQO1 had the largest marker effect and detection-rate gap among the three
eligible genes and ranked first among 2,185 positive discovery markers by
those measures. It was not selected by GBM importance (`56/128`) and is not a
frozen discovery MP4 gene.

In PDC000198, NQO1 abundance was higher in 159 verified tumor–adjacent pairs
(Hodges–Lehmann shift `0.815`, 95% CI `0.648–0.998`; paired Wilcoxon
`P=1.91 × 10⁻²⁰`). Continuous tumor NQO1 was not clearly associated with OS or
RFS. The evidence therefore supports a tissue abundance difference, not
prognosis, macrophage regulation, or causality.

## Outputs

- Preferred figure: `figures/nature/Supplementary_Figure_NQO1_prioritization_and_protein_Nature.*`.
- Compact alternative: `figures/Supplementary_NQO1_prioritization_and_protein.*`.
- `data/processed/OSARS_ROS_overlap_candidate_evidence.csv`: plotted values for the three eligible genes.
- `data/processed/discovery_marker_rank_context.csv`: descriptive full marker screen.
- `data/processed/PDC000198_NQO1_verified_pairs.csv`: verified paired protein observations.
- `tables/NQO1_priority_evidence_summary.csv`: compact numerical audit.

## Concise legend

**Supplementary Fig. X | NQO1 prioritization and paired protein support.**
**a**, the three genes shared by frozen OSARS and the original Hallmark ROS
collection. **b**, the discovery marker screen with NQO1 highlighted.
**c**, NQO1 abundance in 159 verified tumor–adjacent PDC000198 pairs; lines
link specimens from the same case. The marker screen is descriptive, and the
protein panel supports an abundance difference only.

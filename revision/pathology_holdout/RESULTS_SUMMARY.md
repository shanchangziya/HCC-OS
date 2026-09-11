# Pathology strict hold-out summary

1. Final patients: 330.
2. Training: 231 patients/84 events; test: 99 patients/36 events.
3. Features: 2048 initial → 136 after training-only screening → 13 final LASSO features.
4. Apparent training C-index: 0.754 (95% bootstrap CI 0.692–0.809).
5. Strict test C-index: 0.598 (95% bootstrap CI 0.489–0.706; 2,000 patient resamples).
6. Strict test continuous Cox HR per training-set SD: 1.330 (95% CI 0.933–1.896; P=0.115).
7. With the frozen training median cutoff 2.55339, test Low/High n=52/47 and events=18/18; High-versus-Low HR 1.600 (95% CI 0.813–3.150), log-rank P=0.170.
8. Comparison with old analysis: The strict test C-index was modestly lower than the original leakage-prone saved estimate (0.598 vs 0.618), but slightly higher than the prior corrected historical-split estimate (0.590). Memberships differ, so these are descriptive, not paired, comparisons.
9. Interpretation: The strict hold-out results do not provide robust support for a pathology prognostic model; they are more appropriate as an exploratory supplementary analysis.

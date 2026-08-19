# Provenance and superseded analyses

This repository keeps the full analysis trail, including superseded scripts, for transparency. The paper's
results come from the scripts listed in `README.md` and the manuscript's reproducibility section (§8).

## Superseded / exploratory (kept, not used for headline results)

- `ctx_interaction.py` — the initial pilot on cached scGPT representations. It reported a spurious positive
  (z≈+10) driven by heteroscedastic noise and a 65-gene "stopword" panel; documented and retracted in its own
  header. The headline pipeline (`ctx_extract_maxtoki.py` + `ctx_polysemy.py`) fixes those failures.
- `ctx_feasibility.py` — tokenisation-only feasibility check used to choose the extraction settings.
- `ctx_layer_curve.py`, `ctx_devel_trajectory.py`, `ctx_switcher_test.py` — exploratory analyses; the
  switcher test in particular was underpowered (n = 1 usable gene) and is **not** counted among the paper's
  Level-2 nulls.
- `ctx_coexpr_null.py` — **DEPRECATED.** Its null gene-set modules were a fixed 400 genes while the functional
  poles are ~1,880/990 genes; that size mismatch made the "beyond co-expression" z non-reproducible. Replaced by
  `ctx_coexpr_null_v2.py`, which size-matches the null to the functional poles and reports a skew-robust
  empirical p. The corrected result (headline axis empirical p = 0.05; the organisation does not robustly exceed
  co-expression) is what the manuscript reports.

## Determinism

Every script fixes `seed = 0` for both the model/data subsampling and the null draws, so results are
bit-reproducible on the same inputs.

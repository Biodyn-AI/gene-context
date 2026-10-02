# Provenance and superseded analyses

This repository keeps the full analysis trail, including superseded scripts, for transparency. The paper's
results come from the scripts listed in `README.md` and the manuscript's reproducibility section (§8).

## Superseded / exploratory (kept, not used for headline results)

- `ctx_interaction.py` — the initial pilot on cached scGPT representations. It reported a spurious positive
  (z≈+10) driven by heteroscedastic noise and a 65-gene "stopword" panel; documented and retracted in its own
  header. The headline pipeline (`ctx_extract_maxtoki.py` + `ctx_polysemy.py`) fixes those failures.
- `ctx_feasibility.py` — tokenisation-only feasibility check used to choose the extraction settings.
- `ctx_devel_trajectory.py`, `ctx_switcher_test.py` — exploratory analyses; the switcher test in particular was
  underpowered (n = 1 usable gene) and is reported as uninformative, **not** counted among the paper's Level-2
  negatives.
- `ctx_layer_curve.py` — used for Figure 1 (full 12-layer profile, on the `ctxscan_L00..L11` extraction).
  Its stored result had been overwritten by a run on the scGPT tensors (`PREFIX=ctx_scgpt`); that file is kept as
  `results/ctx_layer_curve_ctx_scgpt.json`, and the script now writes non-default prefixes to their own file.
- `ctx_coexpr_null.py` — **DEPRECATED.** Its null gene-set modules were a fixed 400 genes while the functional
  poles are ~1,880/990 genes; that size mismatch made the "beyond co-expression" z non-reproducible.
- `ctx_coexpr_null_v2.py` — size-matched nulls. Its co-expression null (b) uses the genes MOST co-expressed with a
  random seed gene. These modules are 3–10 times more co-expressed than the functional poles, so (b) is a strong
  null, not a matched one (an earlier manuscript version described it as matched; corrected 1 Oct 2026).
- `ctx_coexpr_null_v3.py` — the null matched on size AND average co-expression (coherence within 10% of each
  functional pole), plus a fixed-size power-vs-co-expression curve. All three functional axes exceed it, while
  sitting at the upper edge of v2's strong null. The manuscript reports both.

## Determinism

Every script fixes `seed = 0` for both the model/data subsampling and the null draws.

**Fixed 1 Oct 2026:** `ctx_position_confound.mean_ranks` (used by every rank-controlled analysis) shuffled cells
while iterating over a Python `set` of cell-type names. A set's order depends on `PYTHONHASHSEED`, which is random
per process, so the 1000-cell subset behind each cell type's mean ranks, and therefore every rank-controlled
number, varied by about 0.5% between runs. The loop now iterates in sorted order (checked: identical mean ranks
under two different hash seeds), and all rank-controlled analyses were re-run (`rerun_rankfix.sh`). The result
files from before the fix are kept in `results/pre_rankfix/` for comparison.

## Input-encoding fix and full re-run (1–2 Oct 2026)

While preparing the BMC Genomics submission, an agent review found that every MaxToki input had been encoded wrongly.
The Tabula Sapiens files store log1p(counts per 10,000) in `X`; the extraction code read `X`, normalised it to 10,000
per cell, took log1p a second time and divided by the gene median before ranking. MaxToki's encoding is counts /
gene median with no log. The double transform squashed each cell's expression range, so the gene order followed the
global gene median (Spearman ~0.9) much more than the gene's level in the cell (~0.05); against the correct order,
Spearman 0.6–0.8 and about half of the top-200 genes in common. The Setty extraction also added a log to counts.

Fix: `src/ctx_tokenise.py` (counts from `raw/X`, an integer check, rank by count / median). It gives the same token
sequences as the MaxToki tokenizer on real cells, up to the order of exactly tied genes. All MaxToki extractions and
analyses were re-run by `src/rerun_tokfix.sh` (logs: `results/logs_tokfix/`). The result files from before the fix
are in `results/pre_tokfix/` for comparison. Main changes: functional organisation is now inside both co-expression
nulls; the prediction link is positive (it was slightly negative); steering works for the nuclear/surface axis but
not consistently for mitochondrion/cytoskeleton, and transcription/transport is reversed; MaxToki-1B is not
consistently more contextual than MaxToki-217M once relative depth is matched.

Other changes made in the same round, reviewed before the re-run reached them:
- Null models (`ctx_coexpr_null_v3.py`, `ctx_tightness_null.py`, `ctx_directional_probe.py`): poles matched on size and
  on average co-expression or tightness, disjoint null poles, a fallback for poles below the random level, and Monte
  Carlo p values (k + 1)/(N + 1) (also in `ctx_coexpr_null_v2.py`). The earlier tightness and directional nulls used
  400-gene poles on a fitted curve (not size-matched).
- EXCESS confidence intervals now resample genes (`src/ctx_stats.py`).
- `ctx_position_confound.mean_ranks` uses the same cells as the extraction, in a hash-seed-independent order.
- Random-weights control uses the same 600 cells and 5000-gene panel as the trained comparator; the scaling
  comparison is depth-matched (`ctx_cross_model.py`); the cap-20 control comes from the main extraction's pass.
- Setty: all 5780 cells (the earlier run used a 2515-cell subset cached from another project), spliced + unspliced
  counts, with a spliced-only robustness run.
- Prediction link: the chimeric donor never contains the target and is never padded with the real cell's genes.
- Steering: poles disjoint for the direction; natural movement along the axis recorded next to the push sizes.

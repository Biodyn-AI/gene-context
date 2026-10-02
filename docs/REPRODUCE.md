# Reproduction guide

All analyses are deterministic (`seed = 0`). Two tiers:

- **Tier A — reproduce the figures, tables and manuscript numbers from committed results.** Needs only the
  scientific-Python stack and the committed `results/ctx_*.json`. No models or GB-scale tensors required.
- **Tier B — reproduce everything from scratch.** Needs the public models/datasets (see `docs/DATA.md`), `torch`
  and `transformers`, and ~10 GB of scratch space for the extracted tensors.

## Tier A: figures and manuscript from committed results

```bash
pip install -r requirements.txt
python src/ctx_figures.py                               # figures/ctx_fig1..5.pdf from results/ctx_*.json
bash submission_bmc_genomics/build/build_bmc.sh         # every manuscript number from results/ -> DOCX + PDFs
```

`build_bmc.sh` needs pandoc, python-docx and XeLaTeX. Its first step (`build/make_numbers.py`) collects every number
quoted in the manuscript from the result files; the templates (`*.template.md`) contain only `{{key}}` placeholders,
and the build fails if any placeholder cannot be filled.

## Tier B: full pipeline from the public inputs

1. Edit the path constants at the top of each `src/ctx_extract_*.py` and analysis script to point at your local
   models/data (see `docs/DATA.md`).
2. Run the chain that produced the current results (serial, with a memory guard; about 8 hours on an Apple-silicon
   laptop):

   ```bash
   bash src/rerun_tokfix.sh        # logs in logs_tokfix/ (the committed copies are in results/logs_tokfix/)
   ```

   It runs, in order: the main MaxToki-217M extraction (12 cell types × 1000 cells, 6000 genes, layers 0/4/8/11, plus
   the cap-20 control from the same pass); every analysis on it (EXCESS, anisotropy, rank control, functional axes,
   the three null families, directional congruence, curated targets, lineage TFs, cell-type loadings, the four
   steering configurations, the prediction link); the 217M and 1B scaling extractions (600 cells, 5000 genes; layers
   2/4 and 4/7) and the random-weights extraction (same cells and panel); the cross-model comparison; the Setty
   extraction (all 5780 cells, spliced + unspliced) and its analysis; the full depth scan (all 12 hidden states,
   4000 genes, two passes of six layers); and, at the end, the spliced-only Setty robustness run.
3. The multi-model comparison on identical cells (added 2 Oct 2026; about 9 hours, serial):

   ```bash
   bash src/rerun_multimodel.sh    # one command per line, in the order they were run
   ```

   `ctx_cell_selection.py` records exactly which cells the MaxToki runs used (`results/ctx_cell_selection.npz`).
   On those 600 cells per cell type: `ctx_extract_scgpt_std.py` (scGPT with its standard binned input; also an
   untrained scGPT), `ctx_extract_state_std.py` (STATE SE-600M through its own data collator; also an untrained
   STATE), and `ctx_extract_expression.py` (the expression-only representations; the centroid version used in the paper is `MODE=centroid N_LANDMARK=50 CELL_CAP=50 OUTPREFIX=ctxexprcen50c`; `ctxexprcen50` is the same without the 50-cell cap, a check on sample size; the first, 512-component build `ctxexprcen` failed its positive control and is kept only for the record). `ctx_expression_cellcounts.py` counts the cells behind each uncapped vector. Every analysis then runs on any
   extraction through environment variables read by `src/ctx_prefix.py` (`PREFIX`, `TAPS`/`TAP`, `MIN_CTX`, `COV`,
   `RESEED_PER_TAP`); non-default runs write `results/<analysis>__<prefix>[__minctxN][__cov-expr].json` and never
   overwrite the headline files. `ctx_cross_model.py` takes `MODELS="label:prefix:layer,..."` and `COMMON=1`
   (entries count-balanced in every model). Steering (`ctx_causal.py`: `MODEL=1b`, `ALPHA_UNIT=centred`,
   `ZERO_MASSIVE=1`) and the paired prediction link (`ctx_prediction_link.py`: `PAIRED=1`) cover MaxToki-1B.
   The legacy `ctx_extract_scgpt.py` and `ctx_extract_state.py` (an old activation atlas and preprocessed files with
   non-standard inputs) are kept only for the record; the paper no longer uses them.
4. `ctx_context_composition.py` summarises the dataset, donor and assay of each cell type's cells (quoted in the Limitations).
5. Figures and manuscript as in Tier A (`build/make_additional_file.py` writes Additional file 1).

## Notes

- All MaxToki inputs go through `src/ctx_tokenise.py`: counts from `raw/X` (Tabula Sapiens) or the count layers
  (Setty), ranked by count / gene median. Never feed the log-normalised `X` of a CELLxGENE file to a rank-value
  tokenizer (this was the bug fixed on 1 Oct 2026; see `docs/PROVENANCE.md`).
- `ctx_coexpr_null.py` is **deprecated** (null-set size mismatch). Use `ctx_coexpr_null_v2.py` and `ctx_coexpr_null_v3.py`.
- Rank-controlled analyses (anything calling `ctx_position_confound.mean_ranks`) tokenise the three Tabula Sapiens
  panels and need `transformers` installed (the MaxToki tokenizer adapter imports it).
- Extraction throughput on Apple-silicon MPS (float32), MaxToki-217M: about 3 cells·s⁻¹ at sequence length 1024.
  The 1B extraction of 7200 cells took about 2 hours.

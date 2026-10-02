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
3. scGPT and STATE: `ctx_extract_scgpt.py` reformats an existing activation atlas, and `ctx_extract_state.py` reads
   preprocessed per-tissue files (see `docs/DATA.md`); their outputs (`ctx_scgpt_L*.npz`, `ctx_state_L*.npz`) are read
   by `ctx_cross_model.py`.
4. Figures and manuscript as in Tier A.

## Notes

- All MaxToki inputs go through `src/ctx_tokenise.py`: counts from `raw/X` (Tabula Sapiens) or the count layers
  (Setty), ranked by count / gene median. Never feed the log-normalised `X` of a CELLxGENE file to a rank-value
  tokenizer (this was the bug fixed on 1 Oct 2026; see `docs/PROVENANCE.md`).
- `ctx_coexpr_null.py` is **deprecated** (null-set size mismatch). Use `ctx_coexpr_null_v2.py` and `ctx_coexpr_null_v3.py`.
- Rank-controlled analyses (anything calling `ctx_position_confound.mean_ranks`) tokenise the three Tabula Sapiens
  panels and need `transformers` installed (the MaxToki tokenizer adapter imports it).
- Extraction throughput on Apple-silicon MPS (float32), MaxToki-217M: about 3 cells·s⁻¹ at sequence length 1024.
  The 1B extraction of 7200 cells took about 2 hours.

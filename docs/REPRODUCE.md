# Reproduction guide

All analyses are deterministic (`seed = 0`). Two tiers:

- **Tier A — reproduce the figures/tables from committed results.** Needs only the scientific-Python stack and
  the committed `results/ctx_*.json`. No models or GB-scale tensors required.
- **Tier B — reproduce everything from scratch.** Needs the public models/datasets (see `docs/DATA.md`), `torch`
  and `transformers`, and ~9 GB of scratch space for the extracted tensors.

## Tier A: figures from committed results

```bash
pip install -r requirements.txt
python src/ctx_figures.py            # writes figures/ctx_fig1..5.pdf from results/ctx_*.json
bash paper/make_paper_pdf.sh         # builds the manuscript PDF (needs pandoc + xelatex)
```

## Tier B: full pipeline from the public inputs

1. Edit the path constants at the top of each `src/ctx_extract_*.py` and analysis script to point at your local
   models/data (see `docs/DATA.md`).

2. **Extract** the per-gene contextual representations (writes `results/<prefix>_L*.npz`):

   ```bash
   python src/ctx_extract_maxtoki.py                 # MaxToki-217M, 12 cell types × 1000 cells (headline)
   MODEL=MaxToki-1B-HF  CELLSCTX=600 python src/ctx_extract_maxtoki.py   # 1B, matched scaling pair
   CELLSCTX=600 OUTPREFIX=ctx217m600 python src/ctx_extract_maxtoki.py   # 217M matched pair
   python src/ctx_extract_random.py                  # random-weights control
   python src/ctx_extract_scgpt.py                   # scGPT
   python src/ctx_extract_state.py                   # STATE
   python src/ctx_extract_devel.py                   # Setty developmental
   ```

3. **Analyse** (each writes a `results/ctx_*.json`):

   ```bash
   python src/ctx_polysemy.py            # §4.1  EXCESS by layer
   python src/ctx_anisotropy.py          # §4.1  anisotropy / self-similarity
   python src/ctx_position_confound.py   # §4.2  rank control
   python src/ctx_independent.py         # §4.3  independent-dataset replication
   python src/ctx_functional_axes.py     # §4.4  functional-axis modulation vs random
   python src/ctx_coexpr_null_v2.py      # §4.4  DECISIVE co-expression null (size-matched)
   python src/ctx_tightness_null.py      # §4.4  representation-tightness null
   AXIS=nuc_surf SITE=3 python src/ctx_causal.py     # §4.5  causal (repeat for mito_cyto/trans_transport, SITE=7)
   python src/ctx_directional_probe.py   # §4.6  Level-2 directional congruence
   python src/ctx_curated_targets.py     # §4.6  Level-2 curated targets
   python src/ctx_cross_model.py         # §4.7  cross-model + scaling + random-weights
   python src/ctx_prediction_link.py     # §4.8  context-benefit vs contextualisation
   python src/ctx_celltype_loadings.py   # per-cell-type axis loadings (§4.4)
   ```

4. **Figures + manuscript** as in Tier A.

## Notes

- `ctx_coexpr_null.py` is **deprecated** (a null-set size mismatch made it non-reproducible). Use
  `ctx_coexpr_null_v2.py`. See `docs/PROVENANCE.md`.
- Extraction throughput on Apple-Silicon MPS (float32): ~4.9 / 2.4 / 1.0 cells·s⁻¹ at sequence length
  512 / 1024 / 2048. The headline MaxToki-217M extraction is ~85 min.

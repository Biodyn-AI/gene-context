# gene-context

**Single-cell foundation models represent genes in context, but reveal no context-specific gene function beyond co-expression.**

This repository contains the full analysis code, intermediate result artefacts, figures, and manuscript for a
study of how single-cell foundation models (SCFMs) represent genes *contextually* — a gene's internal vector
inside a particular cell — and what that representation actually encodes.

Author: Ihor Kendiukhov · 2026 · [Biodyn-AI](https://github.com/Biodyn-AI)

---

## What the paper shows

A gene's representation inside an SCFM is reshaped by its cellular context. We measure this in MaxToki-217M in
depth, and compare four models (MaxToki-217M and -1B, scGPT, STATE). The reshaping is:

| property | evidence (MaxToki-217M unless stated) |
|---|---|
| **gene-specific** | same-gene cross-cell agreement +0.759 vs. +0.001 for a gene-shuffled null (EXCESS +0.758, layer 4) |
| **early** | highest at layers 1–2 (+0.847, +0.835), then declines to +0.651 at layer 11 |
| **not rank position** | 99% survives regression on the gene's rank |
| **replicable** | independent dataset (Setty CD34+ bone marrow, 9 clusters): EXCESS +0.713, functional-z +5.0 |
| **partly learned** | untrained model of the same architecture: EXCESS +0.441 vs. +0.740, functional-z +4.7 vs. +14.3 |
| **weakly linked to prediction** | genes that move more gain more from real context (partial ρ +0.161, 95% CI +0.083 to +0.239) |
| **read by the model (one axis)** | steering the nuclear/surface direction shifts logits at *other* genes (p = 0.0081; layer 8 30 of 30 cells); mitochondrion/cytoskeleton: no consistent effect; transcription/transport: reversed |
| **no consistent scaling** | at matched relative depth, MaxToki-1B is not more contextual than MaxToki-217M |
| **architecture-dependent** | functional-z 8.1–18.4 (MaxToki) vs. +3.2 (scGPT) vs. +1.1 (STATE) |

**The boundary.** The functional organisation of the contextualisation is no stronger than that of gene sets matched
on size and average co-expression (p = 0.08–0.38), and modules of strongly co-expressed genes reach more
(p = 0.70–0.80). Two direct tests for *context-specific* gene function are negative; the only signal beyond
co-expression is a small, static closeness of transcription factors to their curated targets (p = 0.005), which is
not stronger where the targets are active.

> **Correction (1–2 Oct 2026).** While preparing the BMC Genomics submission we found that the MaxToki inputs had been
> encoded wrongly: the Tabula Sapiens `X` matrix already holds log1p(counts per 10,000), and the code log-transformed it
> a second time before ranking. All MaxToki inputs now use the model's own encoding (raw counts / gene median,
> `src/ctx_tokenise.py`), and every MaxToki extraction and analysis was re-run (`src/rerun_tokfix.sh`; logs in
> `results/logs_tokfix/`; the pre-fix result files are kept in `results/pre_tokfix/`). Several conclusions changed
> (co-expression comparison, prediction link, steering, scaling). Other fixes in the same round: properly matched and
> size-matched null models, gene-level bootstrap CIs, depth-matched scaling, the full Setty data, the layer peak, STATE's
> cell types, and a rank-control determinism bug. See `docs/PROVENANCE.md`.

See [`paper/PAPER_context_representation.pdf`](paper/PAPER_context_representation.pdf) for the full manuscript, and
[`submission_bmc_genomics/`](submission_bmc_genomics/) for the BMC Genomics submission version.

---

## Repository structure

```
gene-context/
├── README.md                     this file
├── LICENSE                       MIT (code)
├── CITATION.cff                  citation metadata
├── requirements.txt              Python dependencies
├── environment.md                exact environment used for the runs
├── src/                          all analysis and extraction code
│   ├── ctx_tokenise.py           MaxToki input encoding (raw counts / gene median) + co-expression input
│   ├── ctx_stats.py              gene-level bootstrap for EXCESS
│   ├── rerun_tokfix.sh           the serial chain that produced every current result (1 Oct 2026)
│   ├── ctx_extract_*.py          model forward passes -> per-gene contextual representations (.npz)
│   ├── state_loader.py           STATE SE-600M loader used by ctx_extract_state.py
│   ├── ctx_polysemy.py           gene-specific context response (EXCESS) by layer
│   ├── ctx_layer_curve.py        full depth profile (Fig. 1)
│   ├── ctx_anisotropy.py         anisotropy / self-similarity
│   ├── ctx_position_confound.py  rank-position control
│   ├── ctx_independent.py        independent-dataset replication (Setty)
│   ├── ctx_functional_axes.py    functional-axis modulation vs random axes
│   ├── ctx_celltype_loadings.py  cell-type ordering along the nuclear/surface axis
│   ├── ctx_coexpr_null_v2.py     random + strongly co-expressed module nulls
│   ├── ctx_coexpr_null_v3.py     null matched on size AND average co-expression
│   ├── ctx_tightness_null.py     null matched on size and representation tightness
│   ├── ctx_causal.py             steering test (multi-axis / layer)
│   ├── ctx_directional_probe.py  Level 2: directional congruence
│   ├── ctx_curated_targets.py    Level 2: TF placement near curated targets
│   ├── ctx_switcher_test.py      Level 2: lineage TFs (uninformative; too few sampled)
│   ├── ctx_cross_model.py        cross-model, depth-matched scaling and random-weights comparison
│   ├── ctx_prediction_link.py    contextualisation vs benefit of real context for prediction
│   ├── ctx_figures.py            figure generation
│   └── ... (superseded pilots kept for transparency; see docs/PROVENANCE.md)
├── results/                      result artefacts (ctx_*.json) from which every number is computed
│   ├── logs_tokfix/              logs of the re-run chain
│   ├── pre_tokfix/               result files from before the input-encoding fix (superseded; for comparison)
│   └── ts_immune_subset_20000_cell_ids.txt   the 20,000 Tabula Sapiens immune cells used
├── figures/                      publication figures (ctx_fig1..5.pdf)
├── paper/                        general manuscript (.md/.tex/.pdf, generated from the BMC source), bibliography
├── submission_bmc_genomics/      BMC Genomics submission package (DOCX manuscript, figures, cover letter, build)
├── submission/                   earlier Elsevier (CBC) package — SUPERSEDED, predates the corrections
└── docs/                         reproduction guide, data dependencies, provenance
```

## Installation

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The extraction scripts additionally require `torch` and `transformers` (to run the models); the analysis scripts
need only the scientific-Python stack. See [`environment.md`](environment.md).

## Data and models

All inputs are public. The scripts expect local copies of the models and datasets; the exact files and their
public sources are listed in [`docs/DATA.md`](docs/DATA.md). The large per-gene contextual representation tensors
(`results/*.npz`, ~9 GB) are **not** committed — they are regenerable from the public inputs with the
`ctx_extract_*.py` scripts. The small JSON result artefacts that every figure and table is computed from **are**
committed under `results/`.

## Reproducing the analysis

Full step-by-step instructions are in [`docs/REPRODUCE.md`](docs/REPRODUCE.md). In brief: set the data/model
paths (top of each script; see [`docs/DATA.md`](docs/DATA.md)), then run `bash src/rerun_tokfix.sh`, which runs every
MaxToki extraction and analysis in order (about 8 hours on an Apple-silicon laptop) and writes `results/ctx_*.json`.
Then `python src/ctx_figures.py`, and `bash submission_bmc_genomics/build/build_bmc.sh` to regenerate every number in
the manuscript from the result files and rebuild it. All analyses are deterministic (seed 0).

## Citation

If you use this code or data, please cite the manuscript (see [`CITATION.cff`](CITATION.cff)).

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). The manuscript and figures are © the author.

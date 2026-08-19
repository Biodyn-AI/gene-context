# gene-context

**Single-cell foundation models represent genes in context, but the information is co-expression.**

This repository contains the full analysis code, intermediate result artefacts, figures, and manuscript for a
study of how single-cell foundation models (SCFMs) represent genes *contextually* — a gene's internal vector
inside a particular cell — and what that representation actually encodes.

Author: Ihor Kendiukhov · 2026 · [Biodyn-AI](https://github.com/Biodyn-AI)

---

## What the paper shows

A gene's representation inside an SCFM is reshaped by its cellular context. We measure this systematically across
three architectures (MaxToki, scGPT, STATE) and find the reshaping is:

| property | evidence |
|---|---|
| **gene-specific** | same-gene cross-cell agreement +0.759 vs. +0.001 for a gene-shuffled null (EXCESS +0.758) |
| **reproducible** | across independent halves of the cells (main-effect replication +0.996) |
| **not rank position** | 96% survives residualising on per-(gene, context) token rank |
| **replicable** | independent dataset (Setty CD34+ bone marrow), EXCESS +0.58–0.64 |
| **causally used** | steering a functional direction shifts the model's logits at *other* genes (p = 7×10⁻⁴) |
| **learned** | random-weights control: functional organisation collapses +21 → +3.8 |
| **scaling** | MaxToki-217M → 1B: EXCESS +0.740 → +0.881 (matched settings) |
| **architecture-dependent** | functional organisation z = +21 (rank-based) vs. +2.7 (value-binned) vs. +1.3 (protein-token) |

**The boundary (the load-bearing negative).** Against size-matched co-expression nulls, the functional
organisation of that contextualisation does **not robustly exceed co-expression** (headline axis at empirical
p = 0.05; two other axes non-significant), and two independent probes for *novel* context-dependent gene function
are null. The contextualisation is real, learned, and used — but the information it carries is co-expression, not
a gene's changed role in a new cell type.

See [`paper/PAPER_context_representation.pdf`](paper/PAPER_context_representation.pdf) for the full manuscript.

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
│   ├── ctx_extract_*.py          model forward passes → per-gene contextual representations (.npz)
│   ├── ctx_polysemy.py           §4.1  gene-specific context response (EXCESS)
│   ├── ctx_anisotropy.py         §4.1  Ethayarajh anisotropy / self-similarity
│   ├── ctx_position_confound.py  §4.2  rank-position control
│   ├── ctx_independent.py        §4.3  independent-dataset replication
│   ├── ctx_functional_axes.py    §4.4  functional-axis modulation vs random
│   ├── ctx_coexpr_null_v2.py     §4.4  DECISIVE co-expression null (size-matched, skew-robust)
│   ├── ctx_tightness_null.py     §4.4  representation-tightness null
│   ├── ctx_causal.py             §4.5  causal steer-here/read-there (multi-axis/layer)
│   ├── ctx_directional_probe.py  §4.6  Level-2: directional congruence
│   ├── ctx_curated_targets.py    §4.6  Level-2: TF-to-ChIP-target placement
│   ├── ctx_cross_model.py        §4.7  cross-model + scaling + random-weights metric
│   ├── ctx_prediction_link.py    §4.8  does context aid prediction?
│   ├── ctx_figures.py            figure generation
│   └── ... (superseded pilots kept for transparency; see docs/PROVENANCE.md)
├── results/                      intermediate result artefacts (ctx_*.json) — the reproducible outputs
├── figures/                      publication figures (ctx_fig1..5.pdf)
├── paper/                        manuscript (.md/.tex/.pdf), bibliography, build script
├── submission/                   Elsevier submission materials (highlights, cover letter)
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
paths (top of each `ctx_extract_*.py`), run the extractions to produce the `.npz` tensors, then run the analysis
scripts (each writes a `results/ctx_*.json`), then `python src/ctx_figures.py` and `bash paper/make_paper_pdf.sh`.
All analyses are deterministic (seed 0).

## Citation

If you use this code or data, please cite the manuscript (see [`CITATION.cff`](CITATION.cff)).

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). The manuscript and figures are © the author.

# Data and model dependencies

All inputs are **publicly available**. The scripts read them from local paths defined at the top of each file;
set those paths to your local copies. No data are redistributed in this repository, except the list of the
20,000 Tabula Sapiens immune cells used (`results/ts_immune_subset_20000_cell_ids.txt`).

## Models

| Model | Used for | Public source |
|---|---|---|
| **MaxToki-217M / 1B** (Hugging Face safetensors) | main model; extraction, steering, prediction link | Hugging Face `theodoris-lab/MaxToki`; Gómez Ortega et al., bioRxiv 2026, doi 10.64898/2026.03.30.715396 |
| **scGPT** (whole-human checkpoint) | cross-architecture comparison | Cui et al., *Nat. Methods* 2024; github.com/bowang-lab/scGPT |
| **Geneformer** dictionaries (gc104M) | gene-symbol/Ensembl map | Theodoris et al., *Nature* 2023; Hugging Face `ctheodoris/Geneformer` |
| **STATE SE-600M** | cross-architecture comparison | Adduri et al., *Cell* 2026; Hugging Face `arcinstitute/SE-600M` |

Tokenizer/annotation helper files used by the scripts:
- `token_dictionary.json` and the gene-median dictionary — from the MaxToki release (read by `maxtoki_adapter.py`).
- `gene_name_id_dict_gc104M.pkl` — Geneformer gene-symbol ↔ Ensembl map.

## Single-cell datasets

| Dataset | Used for | Public source |
|---|---|---|
| **Tabula Sapiens** immune (20,000-cell subset), kidney, lung | main context panels (12 cell types) | CZ CELLxGENE Discover dataset versions `e62b0182-368f-4875-8ac1-f49e3d5eda60` (immune), `65ca6e36-73b0-4c88-b0f3-7b23b48844ad` (kidney), `40f8b1a3-9f76-4ac4-8761-32078555ed4e` (lung): `https://datasets.cellxgene.cziscience.com/<id>.h5ad` |
| **Setty CD34+ bone marrow** (5780 cells) | independent replication (9 differentiation clusters) | scVelo `bonemarrow()` file, figshare file 27686835 (`https://ndownloader.figshare.com/files/27686835`); Setty et al., *Nat. Biotechnol.* 2019 |

**Important — which matrix to read.** The CELLxGENE Tabula Sapiens files store log1p(counts per 10,000) in `X`
(computed from ambient-corrected counts) and the original integer counts in `raw/X` (same genes). MaxToki's
encoding needs counts, so every MaxToki script reads `raw/X` through `src/ctx_tokenise.py` (which refuses
non-integer input). The Setty/scVelo file has no `raw` group; its `X` holds velocyto **spliced** counts and its
layers hold `spliced` and `unspliced`. The main Setty extraction ranks genes by spliced + unspliced counts (the closer
match to the Tabula Sapiens counts, which include intronic reads); `SETTY_COUNTS=spliced` gives the spliced-only
robustness run.

**STATE inputs.** `ctx_extract_state.py` reads preprocessed per-tissue copies of the same Tabula Sapiens datasets that
keep about 4800 highly variable genes (seurat_v3) with ambient-corrected counts (`tabula_sapiens_*_processed.h5ad`).
These were made in earlier work and are not produced by a script in this repository; they are available from the
author on request.

**scGPT inputs.** `ctx_extract_scgpt.py` reformats an existing per-gene activation atlas of 3000 Tabula Sapiens cells
(from earlier work, Kendiukhov 2026, arXiv 2603.02952), in which each cell's 1200 highest-expressed genes were given
to scGPT with log-normalised values in the value channel (unbinned). The atlas is available from the author on request.

## Gene annotations

| Resource | Used for | Source |
|---|---|---|
| **Gene Ontology** (`gene2go_all.pkl`: symbol → GO term set) | functional-axis poles (direct annotations only) | GEARS file, Harvard Dataverse file 6153417 (`https://dataverse.harvard.edu/api/access/datafile/6153417`) |
| **TRRUST v2** (all human interactions) and **DoRothEA** (confidence A/B) | Level-2 curated-target test | www.grnpedia.org/trrust; github.com/saezlab/dorothea |

## Expected local layout

The scripts as published use absolute paths from the development machine. Before running, edit the path constants
at the top of each script (search for `NAME_ID`, `G2G`, `TS`, `MSETUP`, `MDIR`, `PANELS`, `DATA`, `NET`) to point at
your local copies. A minimal set for the main MaxToki analyses:

```
<models>/MaxToki-217M-HF/                     # config.json + model.safetensors
<models>/MaxToki-1B-HF/
<models>/token_dictionary.json
<geneformer>/gene_name_id_dict_gc104M.pkl
<data>/tabula_sapiens_immune_subset_20000.h5ad
<data>/tabula_sapiens_kidney.h5ad
<data>/tabula_sapiens_lung.h5ad
<data>/setty19_cd34_bm.h5ad                    # scVelo bonemarrow.h5ad
<annot>/gene2go_all.pkl
<networks>/trrust_human.tsv  <networks>/dorothea_human.tsv
```

## Extracted tensors (not committed)

The per-gene contextual representations produced by the extraction scripts are stored as `results/*.npz`
(shape `[partition, context, gene, dim]`, ~200–500 MB each, ~9 GB total). They are excluded from git
(`.gitignore`) because of size and are regenerable from the public inputs above (MaxToki, random-weights, Setty;
STATE given the preprocessed files). They are available from the author on reasonable request.

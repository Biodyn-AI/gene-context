# Data and model dependencies

All inputs are **publicly available**. The scripts read them from local paths defined at the top of each file;
set those paths to your local copies. No data are redistributed in this repository.

## Models

| Model | Used for | Public source |
|---|---|---|
| **MaxToki-217M / 1B** (HF safetensors) | primary model; extraction + steering | Hugging Face `theodoris-lab/MaxToki`; Gómez Ortega & Theodoris, bioRxiv 2026, DOI 10.64898/2026.03.30.715396 |
| **scGPT** | cross-architecture comparison | Cui et al., *Nat. Methods* 2024; public release |
| **Geneformer** | tokenizer / gene-median dictionaries | Theodoris et al., *Nature* 2023; HF `ctheodoris/Geneformer` |
| **STATE / STATE-SE** | cross-architecture comparison | Adduri et al., bioRxiv 2025 (Arc Institute) |

Tokenizer/annotation helper files used by the scripts:
- `token_dictionary.json` — MaxToki gene→token map (from the MaxToki release).
- `gene_name_id_dict_gc104M.pkl` — Geneformer gene-symbol↔Ensembl map.

## Single-cell datasets

| Dataset | Used for | Public source |
|---|---|---|
| **Tabula Sapiens** (immune, kidney, lung) | primary context panels (12 cell types) | Tabula Sapiens Consortium, *Science* 2022; CZ CELLxGENE Discover |
| **Setty CD34+ bone marrow** | independent-dataset replication (developmental states) | Setty et al. 2019; GEO |
| **Replogle K562 / RPE1** | (exploratory only) | Replogle et al. 2022 |

## Gene annotations

| Resource | Used for | Source |
|---|---|---|
| **Gene Ontology** (`gene2go_all.pkl`: symbol → GO term set) | functional-axis pole definitions | Gene Ontology Consortium |
| **DoRothEA / TRRUST / ChIP targets** | Level-2 curated-target test | DoRothEA; TRRUST v2 |

## Expected local layout

The scripts as published use absolute paths from the development machine. Before running, edit the path constants
at the top of each script (search for `NAME_ID`, `G2G`, `TS`, `MSETUP`, `MDIR`, `PANELS`) to point at your local
copies. A minimal set for the headline MaxToki analyses:

```
<models>/MaxToki-217M-HF/                     # config.json + model.safetensors
<models>/token_dictionary.json
<geneformer>/gene_name_id_dict_gc104M.pkl
<data>/tabula_sapiens_immune_subset_20000.h5ad
<data>/tabula_sapiens_kidney.h5ad
<data>/tabula_sapiens_lung.h5ad
<data>/setty19_cd34_bm.h5ad
<annot>/gene2go_all.pkl
```

## Extracted tensors (not committed)

The per-gene contextual representations produced by the extraction scripts are stored as `results/*.npz`
(shape `[partition, context, gene, dim]`, ~200–500 MB each, ~9 GB total). They are excluded from git
(`.gitignore`) because of size and are fully regenerable from the public inputs above. They are available from
the author on reasonable request.

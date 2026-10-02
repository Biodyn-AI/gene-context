# Environment used for the reported runs

- **Hardware/OS:** Apple-silicon Mac (macOS), PyTorch MPS backend, float32.
- **Python:** 3.12.
- **Versions used for the current results (1 Oct 2026 re-run):** torch 2.12.1, transformers 5.13.0, numpy 2.4.6,
  scipy 1.18.0, scikit-learn 1.9.0.
- **Analysis stack:** numpy, scipy, scikit-learn, h5py, matplotlib (see `requirements.txt`).
- **Extraction stack:** torch, transformers (loads the Hugging Face `LlamaForCausalLM` MaxToki checkpoints; MaxToki-1B
  needs transformers >= 4.43). STATE additionally needs the Arc Institute `state` package (arc_state 0.11.1 was used)
  and anndata.
- **Manuscript build:** pandoc 3.x; DOCX via a reference document and python-docx post-processing; PDFs via XeLaTeX.

The analyses are deterministic (`seed = 0`).

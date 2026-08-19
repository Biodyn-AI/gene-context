# Environment used for the reported runs

- **OS:** macOS (Apple Silicon), MPS backend, float32.
- **Python:** 3.12.
- **Analysis stack:** numpy, scipy, scikit-learn, h5py, matplotlib (see `requirements.txt`).
- **Extraction stack:** torch ≥ 2.2, transformers ≥ 4.40 (loads the HF `LlamaForCausalLM` MaxToki checkpoints).
- **Manuscript build:** pandoc 3.x + XeLaTeX (MacTeX), `--citeproc` with `paper/paper_ctx.bib`.

Exact package versions are not pinned to a lockfile; the analyses depend only on stable APIs of the libraries
above and are deterministic (`seed = 0`).

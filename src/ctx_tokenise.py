"""Shared input handling for the context-representation analyses: MaxToki's rank-value encoding and the
expression matrix used for co-expression.

WHY THIS FILE EXISTS (bug fixed 1 Oct 2026). The Tabula Sapiens h5ad files store log1p(counts per 10,000) in `X`
(computed from ambient-corrected counts); the original integer counts are in `raw/X` (same genes, same order).
The earlier code read `X`, normalised those LOG values to 10,000 per cell, applied log1p a SECOND time, and divided
by the gene median. MaxToki's encoding (inherited from Geneformer) is counts normalised to the cell total and divided
by each gene's median, with no log. The double transform squashed each cell's expression range, so the gene order
was set mostly by the global gene median (Spearman ~0.9 with 1/median) rather than by the gene's level in the cell.
The Setty file stores integer counts in `X`, and its extraction added an unneeded log1p as well.

THE CORRECT ORDER. A per-cell scale factor (dividing by the cell total, multiplying by 10,000) does not change the
order of genes within a cell, so the order is simply   raw_count(g) / median(g)   over the cell's non-zero genes,
descending, truncated to max_len - 2 (room for <bos>/<eos>). That is what `rank_order` computes.

CO-EXPRESSION uses the standard log1p(counts per 10,000) of the raw counts (`log_cp10k`).
"""
import numpy as np


def count_matrix(f):
    """The CSR group holding integer counts: `raw/X` when present (CELLxGENE files, where `X` is log-normalised),
    otherwise `X` (e.g. the Setty/scVelo file, whose `X` is counts). Checked to be integer-valued by `check_counts`."""
    if "raw" in f and "X" in f["raw"]:
        return f["raw"]["X"]
    return f["X"]


def check_counts(X, n_rows=5):
    """Raise if the first rows of a CSR group are not non-negative integers (guards against reading log values)."""
    indptr = X["indptr"]
    for r in range(min(n_rows, int(X.attrs["shape"][0]))):
        s, e = int(indptr[r]), int(indptr[r + 1])
        v = np.asarray(X["data"][s:e], dtype=np.float64)
        if len(v) and (v.min() < 0 or not np.allclose(v, np.round(v))):
            raise ValueError("expected integer counts but found non-integer values: is this a log-normalised matrix?")


def rank_order(counts, medians):
    """Indices (into `counts`) of the non-zero genes, in MaxToki's order: descending counts / gene median.
    `counts` and `medians` are aligned 1-D arrays for the genes in the model vocabulary. Stable sort, so ties keep
    their input order (deterministic)."""
    nz = np.nonzero(counts > 0)[0]
    norm = counts[nz] / np.maximum(medians[nz], 1e-9)
    return nz[np.argsort(-norm, kind="stable")]


def log_cp10k(counts, total):
    """log1p(counts per 10,000); `total` is the cell's total count over ALL genes."""
    return np.log1p(counts / (float(total) or 1.0) * 1e4)

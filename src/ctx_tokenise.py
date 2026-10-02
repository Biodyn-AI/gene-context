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


def count_rows(f, rows):
    """Dense float32 COUNT rows (len(rows) x n_genes, in the gene order of f['var']) of a cell-by-gene h5ad.
    Uses `raw/X` when present (its genes must equal `var`'s), else `X` when it holds integers (Setty, aging files),
    else expm1(`X`) for a log1p(CP10k) matrix with no raw counts (Replogle K562): that gives counts up to a per-cell
    scale factor, which does not change MaxToki's order (added 2 Oct 2026, encoding fix for the non-ctx scripts)."""
    import h5py
    rows = np.asarray(rows)
    src, from_raw = f["X"], False
    if "raw" in f and "X" in f["raw"]:
        vk = "_index" if "_index" in f["var"] else "index"
        if not np.array_equal(f["raw"]["var"]["_index"][:], f["var"][vk][:]):
            raise ValueError("raw/var genes differ from var genes")
        src, from_raw = f["raw"]["X"], True
    if isinstance(src, h5py.Group):                       # CSR
        n_g = int(src.attrs["shape"][1]); ip = src["indptr"][:]
        E = np.zeros((len(rows), n_g), np.float32)
        for i, r in enumerate(rows):
            s, e = int(ip[r]), int(ip[r + 1]); E[i, src["indices"][s:e]] = src["data"][s:e]
    else:                                                  # dense
        E = np.stack([np.asarray(src[int(r), :], dtype=np.float32) for r in rows])
    head = E[: min(5, len(E))]
    if not np.allclose(head, np.round(head)):
        if from_raw:
            raise ValueError("raw/X is not integer-valued")
        E = np.expm1(E)                                    # log1p(CP10k) -> CP10k (order-equivalent to counts)
    return E


def maxtoki_order(row, var_idx, medians, max_len):
    """MaxToki/Geneformer gene order for one cell: indices into var_idx (the contract the old loaders used), by
    descending count / gene median, truncated to max_len - 2. `row` must hold counts (or counts up to a cell scale)."""
    return rank_order(np.asarray(row, dtype=np.float64)[var_idx], medians)[: max_len - 2]


def log_cp10k_rows(f, rows):
    """Dense log1p(counts per 10,000) rows, each cell normalised by its total over ALL genes: computed from counts
    (`raw/X`, or an integer `X`), or `X` returned as it is when it already holds log1p(CP10k) and no raw counts exist.
    Replaces the old pattern log1p(X / rowsum * 1e4), which logged an already-logged matrix a second time."""
    import h5py
    rows = np.asarray(rows)
    has_raw = "raw" in f and "X" in f["raw"]
    src = f["raw"]["X"] if has_raw else f["X"]
    if isinstance(src, h5py.Group):
        n_g = int(src.attrs["shape"][1]); ip = src["indptr"][:]
        E = np.zeros((len(rows), n_g), np.float32)
        for i, r in enumerate(rows):
            s, e = int(ip[r]), int(ip[r + 1]); E[i, src["indices"][s:e]] = src["data"][s:e]
    else:
        E = np.stack([np.asarray(src[int(r), :], dtype=np.float32) for r in rows])
    head = E[: min(5, len(E))]
    if np.allclose(head, np.round(head)):                  # counts
        tot = E.sum(1, keepdims=True); tot[tot == 0] = 1
        return np.log1p(E / tot * 1e4).astype(np.float32)
    if has_raw:
        raise ValueError("raw/X is not integer-valued")
    return E                                               # already log1p(CP10k)


def select_rows(f, n, seed, cell_line=None, gene_label=None):
    """Sorted random row indices (seeded) of an h5ad, optionally restricted to obs/cell_line == cell_line and to
    obs/gene containing gene_label (case-insensitive). Added 2 Oct 2026: the Replogle 'K562' files hold four cell
    lines (K562, Jurkat, RPE1, HepG2) and mostly knockdown cells; K562 analyses must select K562 controls."""
    import h5py

    def col(k):
        o = f["obs"][k]
        if isinstance(o, h5py.Group):
            c = np.array([x.decode() if isinstance(x, bytes) else x for x in o["categories"][:]]).astype(str)
            return c[o["codes"][:]]
        return np.array([x.decode() if isinstance(x, bytes) else x for x in o[:]]).astype(str)
    X = f["X"]
    n_rows = int(X.attrs["shape"][0]) if "shape" in X.attrs else X.shape[0]
    keep = np.ones(n_rows, bool)
    if cell_line is not None:
        keep &= np.char.lower(col("cell_line")) == cell_line.lower()
    if gene_label is not None:
        keep &= np.char.find(np.char.lower(col("gene")), gene_label.lower()) >= 0
    idx = np.where(keep)[0]
    if not len(idx):
        raise ValueError(f"no rows match cell_line={cell_line} gene={gene_label}")
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(idx, min(n, len(idx)), replace=False))

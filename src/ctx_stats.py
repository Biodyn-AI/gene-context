"""Shared statistics for the context-representation analyses."""
import numpy as np


def gene_bootstrap_ci(S, Dg, gene_ids, rng, n_boot=2000):
    """95% CI for EXCESS = mean(S) - mean(Dg), pooled over (gene, context-pair) entries, resampling GENES.

    S[i], Dg[i] are the same-gene and different-gene cosines of entry i, and gene_ids[i] is the gene of entry i.
    Each gene enters up to n_pairs entries, so resampling entries (the earlier approach) treats correlated entries as
    independent and gives too-narrow intervals. Here every resampled gene brings all of its entries with it.
    """
    S, Dg = np.asarray(S, float), np.asarray(Dg, float)
    _, inv = np.unique(np.asarray(gene_ids), return_inverse=True)
    sS, sD = np.bincount(inv, S), np.bincount(inv, Dg)
    cnt = np.bincount(inv).astype(float)
    n_genes = len(cnt); out = np.empty(n_boot)
    for b in range(n_boot):
        w = np.bincount(rng.integers(0, n_genes, n_genes), minlength=n_genes)
        out[b] = (w @ sS - w @ sD) / (w @ cnt)
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]

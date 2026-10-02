"""EXPRESSION-ONLY contextual gene representation -- the direct test of whether co-expression alone reproduces the
functional organisation found in the models.

For each cell partition p and cell type c (the SAME cells, contexts and partitions as the 600-cell MaxToki runs, from
results/ctx_cell_selection.npz), gene g is represented by its co-expression profile in those cells: the vector of
Pearson correlations, across the cells of (p, c), between g's log1p(CP10k) expression and that of K landmark genes
(the K most variable panel genes over all selected cells). No model is involved: this is a gene's co-expression
neighbourhood in a given cell type. MODE=pcs instead correlates each gene with the scores of K global principal
components of expression (computed once over all selected cells), a denoised version of the same idea.
MODE=centroid (added 2 Oct 2026 after review: both versions above use correlations WITHIN one cell type, so they cannot
see co-expression BETWEEN cell types, which the pooled co-expression nulls include): gene g in (p, c) is the
expression-weighted mean position of the cells of (p, c) in the K global principal components of all selected cells
pooled, i.e. where in pooled cell space the gene is expressed within that cell type. Genes co-expressed across the
pooled cells get nearby positions, so this construction does respond to pooled co-expression. CELL_CAP=50 (3 Oct 2026)
matches the models' sample size: each vector then averages at most 50 expressing cells, as a model vector averages 50
occurrences; without it the centroid averages all expressing cells (median ~200) and its EXCESS is not comparable.

Layout identical to the model extractions (results/ctxexpr_L00.npz: M[partition, context, gene, K], counts, genes,
contexts, cap), on the SAME gene panel and with the SAME count-balanced mask as results/ctx217m600_L04.npz, so the
identical pipeline (EXCESS, functional axes, co-expression and tightness nulls) runs on it unchanged.
"L00" is only a file-name convention (there are no layers).
"""
import os, sys, json
import numpy as np, h5py
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ctx_tokenise as TK
TS = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/raw"
K = int(os.environ.get("N_LANDMARK", 512))
MODE = os.environ.get("MODE", "landmark")       # landmark: correlation with K landmark genes; pcs: with K global PCs
REF = os.environ.get("REF_PREFIX", "ctx217m600")
OUT = os.environ.get("OUTPREFIX", {"landmark": "ctxexpr", "pcs": "ctxexprpc", "centroid": "ctxexprcen"}[MODE])
# centroid only: use at most CELL_CAP expressing cells per (gene, cell type, partition), the first ones in a seeded random
# cell order, so each vector averages as many cells as a model vector averages occurrences (cap 50). 0 = all cells.
CELL_CAP = int(os.environ.get("CELL_CAP", 0))


def main():
    sel = np.load(os.path.join(HERE, "results", "ctx_cell_selection.npz"), allow_pickle=True)
    cells = sel["cells_600"]; panels = sel["panels"].astype(str); contexts = sel["contexts"].astype(str)
    ref = np.load(os.path.join(HERE, "results", f"{REF}_L04.npz"), allow_pickle=True)
    genes = ref["genes"].astype(str); assert list(ref["contexts"].astype(str)) == list(contexts)
    nC, nP = len(contexts), 2
    # read log1p(CP10k) of raw counts for the selected cells, panel genes only
    E = np.zeros((len(cells), len(genes)), np.float32)
    for pi, p in enumerate(panels):
        rows_here = np.where(cells[:, 0] == pi)[0]
        if not len(rows_here): continue
        with h5py.File(os.path.join(TS, p), "r") as f:
            ens = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]).astype(str)
            col = {e.split(".")[0]: j for j, e in enumerate(ens)}
            take = np.array([col.get(g, -1) for g in genes]); vcol = np.full(len(ens), -1, np.int64)
            for gi, vr in enumerate(take):
                if vr >= 0: vcol[vr] = gi
            X = TK.count_matrix(f); TK.check_counts(X); indptr = X["indptr"][:]
            for k in rows_here:
                r = int(cells[k, 1]); s, e = int(indptr[r]), int(indptr[r + 1])
                ii, vv = X["indices"][s:e], X["data"][s:e].astype(np.float32)
                gj = vcol[ii]; keep = gj >= 0
                E[k, gj[keep]] = TK.log_cp10k(vv[keep], vv.sum())
        print(f"[read] {p}: {len(rows_here)} cells", flush=True)
    land = np.argsort(-E.var(0))[:K]
    if MODE in ("pcs", "centroid"):                           # global PCs of z-scored expression over all cells
        Zg = (E - E.mean(0)) / (E.std(0) + 1e-8)
        U, S, Vt = np.linalg.svd(Zg, full_matrices=False); scores = U[:, :K] * S[:K]
    M = np.full((nP, nC, len(genes), K), np.nan, np.float32)
    for c in range(nC):
        for p in range(nP):
            idx = np.where((cells[:, 2] == c) & (cells[:, 3] == p))[0]
            Z = E[idx]; Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-8)
            if MODE == "centroid":                            # expression-weighted mean cell position (pooled PCs)
                if CELL_CAP > 0:                              # first CELL_CAP expressing cells in a seeded random order
                    perm = np.random.default_rng(1000 * c + p).permutation(len(idx)); idx = idx[perm]
                    W = E[idx]; on = W > 0
                    W = W * (on & (np.cumsum(on, 0) <= CELL_CAP))
                else:
                    W = E[idx]
                tot = W.sum(0)
                M[p, c] = (W.T @ scores[idx]) / np.maximum(tot, 1e-8)[:, None]
                M[p, c][tot <= 0] = np.nan
            elif MODE == "pcs":
                R = scores[idx]; R = (R - R.mean(0)) / (R.std(0) + 1e-8)
                M[p, c] = (Z.T @ R) / len(idx)                # correlation of every gene with every global PC
            else:
                M[p, c] = (Z.T @ Z[:, land]) / len(idx)      # correlation of every gene with every landmark
        print(f"[corr] context {c} ({contexts[c][:30]}): {len(idx)} cells per partition", flush=True)
    counts = ref["counts"]; cap = int(ref["cap"])            # same count-balanced (gene, context) entries as MaxToki
    if MODE == "centroid":                                    # float16 range: put PC scores on a unit scale
        M = M / (np.nanstd(M) + 1e-8)
    np.savez_compressed(os.path.join(HERE, "results", f"{OUT}_L00.npz"), M=np.nan_to_num(M).astype(np.float16),
                        counts=counts, genes=genes, contexts=contexts, cap=cap, landmarks=genes[land], cell_cap=CELL_CAP)
    print(f"[done] results/{OUT}_L00.npz  M{M.shape}", flush=True)


if __name__ == "__main__":
    main()

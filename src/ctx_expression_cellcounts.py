"""How many cells does each expression-only vector average? For every count-balanced (gene, cell type, partition)
entry of the 600-cell runs (mask from results/ctx217m600_L04.npz), count the selected cells that express the gene
(raw count > 0). A model vector averages exactly 50 occurrences per entry (the cap); the uncapped centroid version
averages all expressing cells, so this shows how unequal the sample sizes were. Out: results/ctx_expression_cellcounts.json
"""
import os, sys, json
import numpy as np, h5py
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ctx_tokenise as TK
TS = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/raw"


def main():
    sel = np.load(os.path.join(HERE, "results", "ctx_cell_selection.npz"), allow_pickle=True)
    cells = sel["cells_600"]; panels = sel["panels"].astype(str)
    ref = np.load(os.path.join(HERE, "results", "ctx217m600_L04.npz"), allow_pickle=True)
    genes = ref["genes"].astype(str); counts, cap = ref["counts"], int(ref["cap"])
    nP, nC, nG = counts.shape
    n_expr = np.zeros((nP, nC, nG), np.int32)
    for pi, p in enumerate(panels):
        rows = cells[cells[:, 0] == pi]
        if not len(rows): continue
        with h5py.File(os.path.join(TS, p), "r") as f:
            ens = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]).astype(str)
            col = {e.split(".")[0]: j for j, e in enumerate(ens)}
            vcol = np.full(len(ens), -1, np.int64)
            for gi, g in enumerate(genes):
                if g in col: vcol[col[g]] = gi
            X = TK.count_matrix(f); TK.check_counts(X); indptr = X["indptr"][:]
            for _, r, c, part in rows:
                s, e = int(indptr[r]), int(indptr[r + 1])
                ii, vv = X["indices"][s:e], X["data"][s:e]
                gj = vcol[ii]; ok = (gj >= 0) & (vv > 0)
                n_expr[part, c, gj[ok]] += 1
    bal = (counts == cap).all(0)                              # (nC, nG) count-balanced entries
    vals = np.concatenate([n_expr[p][bal] for p in range(nP)])
    out = dict(n_entries=int(bal.sum()), model_occurrences_per_entry=cap,
               expressing_cells_per_entry=dict(median=float(np.median(vals)), p10=float(np.percentile(vals, 10)),
                                               p90=float(np.percentile(vals, 90)), min=int(vals.min()), max=int(vals.max())))
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, "results", "ctx_expression_cellcounts.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

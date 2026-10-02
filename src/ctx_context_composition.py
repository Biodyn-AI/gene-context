"""What else differs between the 12 cell types besides cell type: source dataset, tissue, donor and assay of the
selected cells (results/ctx_cell_selection.npz; 1000-cell main selection and the 600-cell selection). Both cell
partitions of a cell type share these, so differences between cell types in them count toward EXCESS.

Out: results/ctx_context_composition.json (per cell type: counts by dataset, tissue, assay and donor; the largest
dataset / donor share; the 10x 3' v3 share) and a short summary.
"""
import os, json
from collections import Counter
import numpy as np, h5py
HERE = os.path.dirname(os.path.abspath(__file__))
TS = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/raw"
COLS = ("donor_id", "assay", "tissue")


def obs_col(f, c):
    g = f["obs"][c]
    if isinstance(g, h5py.Group):                                  # anndata categorical: codes + categories
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]).astype(str)
        return cats[g["codes"][:]]
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v]).astype(str)


def main():
    sel = np.load(os.path.join(HERE, "results", "ctx_cell_selection.npz"), allow_pickle=True)
    panels = sel["panels"].astype(str); contexts = sel["contexts"].astype(str)
    meta = {}
    for pi, p in enumerate(panels):
        with h5py.File(os.path.join(TS, p), "r") as f:
            meta[pi] = {c: obs_col(f, c) for c in COLS}
    out = {"note": __doc__.strip().splitlines()[0]}
    for key in ("cells_1000", "cells_600"):
        cells = sel[key]; res = {}
        for ci, ctx in enumerate(contexts):
            rows = cells[cells[:, 2] == ci]
            ds = Counter(panels[r[0]].replace("tabula_sapiens_", "").replace(".h5ad", "") for r in rows)
            per = {c: Counter(meta[r[0]][c][r[1]] for r in rows) for c in COLS}
            n = len(rows)
            res[ctx] = dict(n=n, datasets=dict(ds), tissues=dict(per["tissue"]), assays=dict(per["assay"]),
                            n_donors=len(per["donor_id"]), donors=dict(per["donor_id"].most_common()),
                            top_dataset_share=max(ds.values()) / n, top_donor=per["donor_id"].most_common(1)[0][0],
                            top_donor_n=per["donor_id"].most_common(1)[0][1],
                            top_donor_share=per["donor_id"].most_common(1)[0][1] / n,
                            share_10x_3v3=sum(v for k, v in per["assay"].items() if k == "10x 3' v3") / n)
        one_ds = [c for c, r in res.items() if len(r["datasets"]) == 1]
        out[key] = dict(per_cell_type=res, cell_types_from_one_dataset=one_ds,
                        cell_types_one_dataset_by_dataset={c: list(res[c]["datasets"])[0] for c in one_ds},
                        min_share_10x_3v3=min(r["share_10x_3v3"] for r in res.values()),
                        max_top_donor_share=max(r["top_donor_share"] for r in res.values()),
                        cell_types_top_donor_ge_half={c: f"{r['top_donor_n']} of {r['n']}" for c, r in res.items()
                                                      if r["top_donor_share"] >= 0.5})
        print(f"[{key}] one dataset only: {out[key]['cell_types_one_dataset_by_dataset']}")
        print(f"[{key}] min 10x 3' v3 share {out[key]['min_share_10x_3v3']:.3f}; "
              f"top donor >= half: {out[key]['cell_types_top_donor_ge_half']}")
        for c, r in res.items():
            print(f"   {c[:34]:<34} datasets {r['datasets']}  donors {r['n_donors']}  top donor {r['top_donor_n']}/{r['n']}"
                  f"  10x3'v3 {r['share_10x_3v3']:.2f}  assays {r['assays']}")
    json.dump(out, open(os.path.join(HERE, "results", "ctx_context_composition.json"), "w"), indent=1)
    print("[done] -> results/ctx_context_composition.json")


if __name__ == "__main__":
    main()

"""INDEPENDENT-DATASET REPLICATION of Level 1. All context claims so far rest on Tabula Sapiens (adult cell
types). This replicates EXCESS (contextualisation) and FUNC-Z (functional organisation) on a DIFFERENT dataset
and a DIFFERENT context axis: Setty CD34+ bone-marrow hematopoiesis (ctx_devel; all 5780 cells, 9 differentiation
clusters with >= 130 cells, Mega dropped), 2 partitions, cap 25, floor 15. PREFIX=ctx_devel_spliced reads the
spliced-only robustness extraction. If the phenomenon replicates here it is not a
Tabula-Sapiens or adult-cell-type artefact.

Out: results/ctx_independent.json
"""
import os, sys, json, pickle, itertools, warnings; warnings.filterwarnings("ignore")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
PREFIX = os.environ.get("PREFIX", "ctx_devel")
NAME_ID = "/Volumes/Crucial X6/MacBook/Code/neuro-mechinterp/models/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl"
G2G = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
AXES = {"nuclear_vs_surface": (["GO:0005634", "GO:0000785", "GO:0003677"], ["GO:0005886", "GO:0005576", "GO:0005615"]),
        "mito_vs_cytoskeleton": (["GO:0005739"], ["GO:0005856"])}
SEED, N_RANDOM = 0, 200


def cos_rows(A, B):
    A = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-9); B = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-9)
    return (A * B).sum(1)


def main():
    ens2sym = {e: s.upper() for s, e in pickle.load(open(NAME_ID, "rb")).items()}
    g2g = {k.upper(): set(v) for k, v in pickle.load(open(G2G, "rb")).items() if isinstance(v, (set, list, tuple))}
    rng = np.random.default_rng(SEED)
    out = {"dataset": "Setty CD34+ bone marrow (developmental states)", "taps": {}}
    TAPS = [int(x) for x in os.environ.get("TAPS", "2,4").split(",")]   # verdict uses the LAST listed tap
    for tap in TAPS:
        z = np.load(os.path.join(RES, f"{PREFIX}_L{tap:02d}.npz"), allow_pickle=True)
        M, counts, cap = z["M"].astype(np.float32), z["counts"], int(z["cap"])
        genes = z["genes"].astype(str); ctxs = (z["clusters"] if "clusters" in z else z["contexts"]).astype(str)
        nP, nC, nG, d = M.shape
        full = (counts == cap).all(0)
        flat = M[:, full]; mu = flat.reshape(-1, d).mean(0); sd = flat.reshape(-1, d).std(0) + 1e-6
        Mz = (M - mu) / sd
        # EXCESS
        S, D, mr, GI = [], [], [], []
        for c1, c2 in itertools.combinations(range(nC), 2):
            keep = full[c1] & full[c2]
            if keep.sum() < 150:
                continue
            D0 = Mz[0, c2, keep] - Mz[0, c1, keep]; D1 = Mz[1, c2, keep] - Mz[1, c1, keep]
            b0, b1 = D0.mean(0), D1.mean(0); mr.append(float(np.dot(b0, b1) / (np.linalg.norm(b0) * np.linalg.norm(b1) + 1e-9)))
            d0, d1 = D0 - b0, D1 - b1
            S.append(cos_rows(d0, d1)); D.append(cos_rows(d0, d1[rng.permutation(len(d1))])); GI.append(np.where(keep)[0])
        S, D = np.concatenate(S), np.concatenate(D); excess = float(S.mean() - D.mean())
        from ctx_stats import gene_bootstrap_ci            # resample GENES (entries of a gene are correlated)
        ci_g = gene_bootstrap_ci(S, D, np.concatenate(GI), rng, 2000)
        # FUNC-Z
        a_space = np.full((nG, d), np.nan, np.float32)
        for gi in range(nG):
            cs = np.where(full[:, gi])[0]
            if len(cs): a_space[gi] = Mz[:, cs, gi].mean((0, 1))
        ok = np.isfinite(a_space[:, 0]); use = np.where(full.sum(0) >= max(2, nC // 2))[0]
        syms = [ens2sym.get(g) for g in genes]
        def power(vec):
            q = np.tensordot(Mz[:, :, use], vec, axes=([3], [0])); qm = np.where(full[:, use][None], q, np.nan)
            I = qm - np.nanmean(qm, 1, keepdims=True) - np.nanmean(qm, 2, keepdims=True) + np.nanmean(qm, (1, 2), keepdims=True)
            sel = np.isfinite(I[0]) & np.isfinite(I[1]); return float(np.nanmean(I[0][sel] * I[1][sel]))
        def axis(ia, ib):
            u = a_space[ia].mean(0) - a_space[ib].mean(0); return u / (np.linalg.norm(u) + 1e-9)
        pool = np.where(ok)[0]; fz = {}
        for name, (A, B) in AXES.items():
            ia = [i for i, s in enumerate(syms) if ok[i] and s in g2g and g2g[s] & set(A)]
            ib = [i for i, s in enumerate(syms) if ok[i] and s in g2g and g2g[s] & set(B)]
            both = set(ia) & set(ib); ia = [i for i in ia if i not in both]; ib = [i for i in ib if i not in both]
            if len(ia) < 15 or len(ib) < 15: fz[name] = None; continue
            pf = power(axis(ia, ib))
            null = np.array([power(axis(rng.choice(pool, len(ia), False), rng.choice(pool, len(ib), False))) for _ in range(N_RANDOM)])
            fz[name] = dict(z=float((pf - null.mean()) / (null.std() + 1e-12)), poleA=len(ia), poleB=len(ib))
        out["taps"][f"L{tap:02d}"] = dict(n_ctx=int(nC), n_genes=int(nG), contexts=[str(c) for c in ctxs],
                                          pairs_scored=len(mr), excess=excess,
                                          excess_ci=ci_g,
                                          diff_null=float(D.mean()), main_rep=float(np.mean(mr)), func_z=fz)
        fzp = fz.get("nuclear_vs_surface")
        print(f"L{tap}: EXCESS {excess:+.4f} CI[{ci_g[0]:+.3f},{ci_g[1]:+.3f}] "
              f"diff {D.mean():+.4f} mainrep {np.mean(mr):+.3f} | FUNC-z(nuc/surf) {fzp['z']:+.1f}" if fzp else "")
    d4 = out["taps"][f"L{TAPS[-1]:02d}"]; fzp = d4["func_z"].get("nuclear_vs_surface")
    rep = fzp and d4["excess"] > 0 and d4["excess_ci"][0] > 0 and fzp["z"] > 3
    out["prefix"] = PREFIX
    out["verdict"] = ((f"REPLICATES on the Setty data: EXCESS {d4['excess']:+.3f} (CI lower {d4['excess_ci'][0]:+.3f}), "
                       f"diff-null {d4['diff_null']:+.4f}, functional-z {fzp['z']:+.1f}.") if rep else
                      (f"DOES NOT fully replicate: EXCESS {d4['excess']:+.3f}, functional-z "
                       f"{(fzp or {}).get('z', float('nan')):+.1f} (needs EXCESS CI > 0 and functional-z > 3)."))
    print("VERDICT:", out["verdict"])
    json.dump(out, open(os.path.join(RES, "ctx_independent.json" if PREFIX == "ctx_devel" else f"ctx_independent_{PREFIX}.json"), "w"), indent=1)
    print(f"[done] -> results/{'ctx_independent.json' if PREFIX == 'ctx_devel' else f'ctx_independent_{PREFIX}.json'}")


if __name__ == "__main__":
    main()

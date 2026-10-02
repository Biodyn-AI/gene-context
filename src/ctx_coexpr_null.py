"""IS THE FUNCTIONAL-AXIS CONTEXT MODULATION BEYOND CO-EXPRESSION? (the decisive control)

`ctx_functional_axes.py` found that context modulates genes along functional axes (nuclear-vs-surface etc.) far
more than along RANDOM-PARTITION axes (z ~ +20, surviving rank/abundance control). But functional gene sets are
CO-EXPRESSION MODULES -- genes in a pathway are co-regulated -- and every other pipeline in this project found
SCFM attention captures co-expression, not deeper function. So the random-partition null is too weak: random
genes are not co-expressed, so of course a coherent functional set beats them. The honest question is whether
functional axes beat a CO-EXPRESSION-MATCHED null.

DESIGN -- a dose-response, not a single matched null (more robust to imperfect matching). Co-expression
COHERENCE of a gene set = its mean pairwise expression correlation, measured on the same Tabula Sapiens cells.
We generate many null axes spanning the whole coherence range:
    - loose (near-random) gene sets  -> low coherence
    - tight modules grown from a seed's top co-expressed neighbours -> high coherence
    - everything between, via a looseness knob
For each null axis we record (coherence, interaction power). This traces the CO-EXPRESSION BASELINE CURVE:
how much reproducible gene x context modulation you get purely as a function of how co-expressed the axis genes
are. Then we place each FUNCTIONAL axis on the same plot.

READING IT
    functional point ON the curve   -> modulation is fully explained by co-expression coherence; the functional
                                        labels add nothing. This is the project's recurring deflation.
    functional point ABOVE the curve-> functional structure beyond co-expression: context modulates genes along
                                        the functional direction more than equally-co-expressed non-functional
                                        directions do.
The statistic is the functional axis's residual above a power~coherence fit, in units of the null scatter (z).

All power is measured RANK-CONTROLLED (projection residualised on per-(gene,context) mean rank) to keep the
abundance control from ctx_functional_axes.

Out: results/ctx_coexpr_null.json
"""
import os, sys, json, pickle, warnings; warnings.filterwarnings("ignore")
import numpy as np, h5py

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ctx_position_confound as CP
RES = os.path.join(HERE, "results")
NAME_ID = "/Volumes/Crucial X6/MacBook/Code/neuro-mechinterp/models/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl"
G2G = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
TS = CP.TS; PANELS = CP.PANELS
TAPS = [int(x) for x in os.environ.get("TAPS","4,8").split(",")]
PREFIX = os.environ.get("PREFIX", "ctx_maxtoki")
MIN_CTX = int(os.environ.get("MIN_CTX", 9))
N_NULL = 400
COEXPR_CELLS = 6000
SEED = 0
from scipy.stats import spearmanr

AXES = {
    "nuclear_vs_surface":   (["GO:0005634", "GO:0000785", "GO:0003677"], ["GO:0005886", "GO:0005576", "GO:0005615"]),
    "mito_vs_cytoskeleton": (["GO:0005739"], ["GO:0005856"]),
    "transcription_vs_transport": (["GO:0006355", "GO:0003700"], ["GO:0006811", "GO:0038023"]),
}


def coexpr_matrix(genes_ens):
    """gene x gene expression correlation on pooled TS cells (the co-expression we are controlling FOR)."""
    G = len(genes_ens)
    rows = []
    for p in PANELS:
        path = os.path.join(TS, p)
        if not os.path.exists(path):
            continue
        with h5py.File(path, "r") as f:
            ens = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]).astype(str)
            ens = np.array([e.split(".")[0] for e in ens])
            col = {e: j for j, e in enumerate(ens)}
            take = np.array([col.get(g.split(".")[0], -1) for g in genes_ens])   # panel col -> var row
            vcol = np.full(len(ens), -1, np.int64)                                # var row -> panel col
            for pj, vr in enumerate(take):
                if vr >= 0:
                    vcol[vr] = pj
            import ctx_tokenise as TK                       # raw counts -> log1p(CP10k) (X itself is already log)
            X = TK.count_matrix(f); TK.check_counts(X); n = int(X.attrs["shape"][0]); indptr = X["indptr"][:]
            sel = np.sort(np.random.default_rng(SEED).choice(n, min(COEXPR_CELLS // 2, n), replace=False))
            E = np.zeros((len(sel), G), np.float32)
            for i, r in enumerate(sel):
                s, e = int(indptr[r]), int(indptr[r + 1])
                ii, vv = X["indices"][s:e], X["data"][s:e].astype(np.float32)
                pj = vcol[ii]; keep = pj >= 0
                if not keep.any():
                    continue
                E[i, pj[keep]] = TK.log_cp10k(vv[keep], vv.sum())   # total over ALL genes of the raw row
            rows.append(E)
    Z = np.vstack(rows)
    Z = Z - Z.mean(0); Z = Z / (Z.std(0) + 1e-8)
    return ((Z.T @ Z) / len(Z)).astype(np.float32)


def coherence(C, idx):
    if len(idx) < 2:
        return 0.0
    sub = C[np.ix_(idx, idx)]
    iu = np.triu_indices(len(idx), 1)
    return float(sub[iu].mean())


def main():
    ens2sym = {e: s.upper() for s, e in pickle.load(open(NAME_ID, "rb")).items()}
    g2g = {k.upper(): set(v) for k, v in pickle.load(open(G2G, "rb")).items() if isinstance(v, (set, list, tuple))}
    rng = np.random.default_rng(SEED)

    z0 = np.load(os.path.join(RES, f"{PREFIX}_L{TAPS[0]:02d}.npz"), allow_pickle=True)
    ctxs = z0["contexts"].astype(str); genes = z0["genes"].astype(str)
    syms = [ens2sym.get(g) for g in genes]
    tokmap = json.load(open(f"{CP.MSETUP}/token_dictionary.json")); ens2tid = {k: int(v) for k, v in tokmap.items()}
    tids = np.array([ens2tid.get(g, -1) for g in genes])

    print("[1/3] per-context mean rank (abundance control)", flush=True)
    MR = CP.mean_ranks(set(ctxs))
    rank = np.full((len(ctxs), len(genes)), np.nan)
    for ci, c in enumerate(ctxs):
        d = MR.get(c, {})
        for gi, t in enumerate(tids):
            if t in d:
                rank[ci, gi] = d[t]

    print("[2/3] co-expression matrix on Tabula Sapiens cells", flush=True)
    C = coexpr_matrix(list(genes))

    out = {"taps": {}}
    for tap in TAPS:
        z = np.load(os.path.join(RES, f"{PREFIX}_L{tap:02d}.npz"), allow_pickle=True)
        M, counts, cap = z["M"].astype(np.float32), z["counts"], int(z["cap"])
        full = (counts == cap).all(0)
        flat = M[:, full]; mu = flat.reshape(-1, M.shape[-1]).mean(0); sd = flat.reshape(-1, M.shape[-1]).std(0) + 1e-6
        Mz = (M - mu) / sd
        a_space = np.full((len(genes), M.shape[-1]), np.nan, np.float32)
        for gi in range(len(genes)):
            cs = np.where(full[:, gi])[0]
            if len(cs):
                a_space[gi] = Mz[:, cs, gi].mean((0, 1))
        gene_ok = np.isfinite(a_space[:, 0])
        use = np.where(full.sum(0) >= MIN_CTX)[0]
        Ruse = rank[:, use]; muse = full[:, use]
        Msub = Mz[:, :, use]

        def power(vec):
            q = np.tensordot(Msub, vec, axes=([3], [0]))
            qm = np.where(muse[None], q, np.nan)
            fin = np.isfinite(qm[0]) & np.isfinite(Ruse)
            r = Ruse[fin]; A = np.column_stack([np.ones_like(r), r, r ** 2])
            for p in range(2):
                yv = qm[p][fin]; qm[p][fin] = yv - A @ np.linalg.lstsq(A, yv, rcond=None)[0]
            I = qm - np.nanmean(qm, 1, keepdims=True) - np.nanmean(qm, 2, keepdims=True) + np.nanmean(qm, (1, 2), keepdims=True)
            sel = np.isfinite(I[0]) & np.isfinite(I[1])
            return float(np.nanmean(I[0][sel] * I[1][sel]))

        def axis(ia, ib):
            u = a_space[ia].mean(0) - a_space[ib].mean(0); return u / (np.linalg.norm(u) + 1e-9)

        pool = np.where(gene_ok)[0]

        def module(size, loose):
            """a gene set of `size` with tunable co-expression coherence via looseness knob."""
            seed = rng.choice(pool)
            nb = pool[np.argsort(-C[seed, pool])]           # panel genes by co-expression to the seed
            top = nb[: min(len(nb), int(size * loose))]
            return rng.choice(top, min(size, len(top)), replace=False)

        # ----- co-expression baseline curve: many null axes spanning coherence -----
        print(f"\n=== layer {tap}: building co-expression baseline curve ({N_NULL} null axes) ===", flush=True)
        # size the null poles like a representative functional axis
        iaN, ibN = 400, 400   # DEPRECATED: size mismatch vs functional poles inflated the z. Use ctx_coexpr_null_v2.py
        nx, ny = [], []
        for k in range(N_NULL):
            loose = rng.choice([1, 1, 2, 3, 5, 8, 15, 40, 120])   # spans tight module -> near-random
            ga, gb = module(iaN, loose), module(ibN, loose)
            coh = 0.5 * (coherence(C, ga) + coherence(C, gb))
            nx.append(coh); ny.append(power(axis(ga, gb)))
        nx, ny = np.array(nx), np.array(ny)
        # quadratic fit power ~ coherence, residual scatter for the z
        A = np.column_stack([np.ones_like(nx), nx, nx ** 2])
        coef, *_ = np.linalg.lstsq(A, ny, rcond=None)
        resid_sd = float((ny - A @ coef).std())
        rho_cc = float(spearmanr(nx, ny).statistic)
        print(f"   null coherence range [{nx.min():+.3f}, {nx.max():+.3f}]; power~coherence Spearman {rho_cc:+.2f}"
              f"; residual scatter {resid_sd:.3f}")

        # dump the null scatter + fit so the paper figure can draw a faithful band (rank-controlled power)
        if tap == 4:
            np.savez(os.path.join(RES, "ctx_coexpr_null_scatter.npz"),
                     null_coh=nx, null_pow=ny, coef=coef, resid_sd=resid_sd)
        out["taps"][f"L{tap:02d}"] = {"coexpr_power_coherence_rho": rho_cc, "axes": {}}
        for name, (Ag, Bg) in AXES.items():
            ia = [i for i, s in enumerate(syms) if gene_ok[i] and s in g2g and g2g[s] & set(Ag)]
            ib = [i for i, s in enumerate(syms) if gene_ok[i] and s in g2g and g2g[s] & set(Bg)]
            both = set(ia) & set(ib); ia = [i for i in ia if i not in both]; ib = [i for i in ib if i not in both]
            coh = 0.5 * (coherence(C, ia) + coherence(C, ib))
            pw = power(axis(ia, ib))
            pred = float(np.array([1, coh, coh ** 2]) @ coef)
            z = (pw - pred) / (resid_sd + 1e-9)
            # where does its coherence sit among null axes?
            coh_pctl = float((nx < coh).mean())
            print(f"   {name:<28} coherence {coh:+.3f} (pctl {coh_pctl:.2f})  power {pw:+.3f}  "
                  f"curve-predicts {pred:+.3f}  -> ABOVE-CURVE z = {z:+.1f}")
            out["taps"][f"L{tap:02d}"]["axes"][name] = dict(coherence=coh, coherence_pctl=coh_pctl,
                                                            power=pw, curve_pred=pred, z_above_coexpr=float(z))

    allz = [(t, n, d["z_above_coexpr"]) for t, tv in out["taps"].items() for n, d in tv["axes"].items()]
    best = max(allz, key=lambda x: x[2])
    out["verdict"] = (
        f"strongest: {best[0]}/{best[1]}, {best[2]:+.1f} sigma above the co-expression-coherence curve. " +
        ("BEYOND CO-EXPRESSION — functional axes carry context modulation exceeding what equally-co-expressed "
         "non-functional gene sets produce. The model's per-context gene movement has functional structure not "
         "reducible to co-expression coherence."
         if best[2] > 3 else
         "NOT BEYOND CO-EXPRESSION — the functional axes sit on the co-expression-coherence curve. The "
         "context modulation is explained by how co-expressed the axis genes are, not by their function. This "
         "is the project's recurring result: the model tracks co-expression, and 'functional' structure is "
         "co-expression structure wearing a GO label.")
    )
    print(f"\nVERDICT: {out['verdict']}")
    json.dump(out, open(os.path.join(RES, "ctx_coexpr_null.json"), "w"), indent=1)
    print("[done] -> results/ctx_coexpr_null.json")


if __name__ == "__main__":
    main()

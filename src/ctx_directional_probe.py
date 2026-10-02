"""LEVEL-2 FEASIBILITY: does a gene move toward its OWN functional pole in contexts that share that function?

The results so far (ctx_polysemy .. ctx_coexpr_null) reach Level 1: genes undergo gene-specific, reproducible,
abundance-independent context modulation organised along functional axes beyond co-expression. Level 2 -- "the
model tells you a gene's role in a new cell" -- needs the movement to be DIRECTIONALLY appropriate: a
functionally-surface gene should move toward the surface pole specifically in surface/secretory cell types, not
everywhere. This probe asks, cheaply, on the existing extraction, whether any such directional signal exists.

THE STATISTIC -- a cross-validated Tukey non-additivity test along a functional axis u.
For a gene g and context c, project the representation onto u:  q(g,c) = <v(g,c), u>. Decompose into
  f(g)  gene loading    (mean over contexts)   -- the gene's own functional identity on u
  h(c)  context loading (mean over genes)      -- the context's functional character on u
  e(g,c) interaction    (what is left)
The DIRECTIONAL / context-appropriate prediction is that the interaction follows the PRODUCT of the loadings:
e(g,c) ~ f(g) * h(c), with positive slope -- i.e. a gene moves further toward its own pole exactly in contexts
whose character shares that pole. beta = corr(e, f*h) measures this.

NO CIRCULARITY: f, h (hence the predictor f*h) are built from ONE cell partition; the interaction e is measured
on the OTHER partition; beta averages the two cross-directions. A predictor fit on one half cannot manufacture
correlation with an independent half's residual.

CONTROLS: projection rank-residualised (abundance); and the decisive one -- functional beta is compared against 300
null axes whose poles match the functional poles on SIZE and average co-expression (coherence within 10%), as in
ctx_coexpr_null_v3; empirical one-sided p. (Rewritten 1 Oct 2026: the earlier version used 400-gene null poles and a
fitted beta-vs-coherence curve, i.e. not size-matched.)

Reads the headline ctx_maxtoki_L*.npz. Out: results/ctx_directional_probe.json
"""
import os, sys, json, pickle, warnings; warnings.filterwarnings("ignore")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ctx_position_confound as CP
from ctx_coexpr_null import coexpr_matrix, coherence
RES = os.path.join(HERE, "results")
NAME_ID = "/Volumes/Crucial X6/MacBook/Code/neuro-mechinterp/models/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl"
G2G = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
TAPS = [4, 8]   # headline ctx_maxtoki set has L0/L4/L8/L11; L4 is where the null comparisons were run
MIN_CTX = 9
N_NULL = 300
TOL = 0.10
MAX_TRIES = 40
SEED = 0
from scipy.stats import spearmanr

AXES = {
    "nuclear_vs_surface":   (["GO:0005634", "GO:0000785", "GO:0003677"], ["GO:0005886", "GO:0005576", "GO:0005615"]),
    "mito_vs_cytoskeleton": (["GO:0005739"], ["GO:0005856"]),
    "transcription_vs_transport": (["GO:0006355", "GO:0003700"], ["GO:0006811", "GO:0038023"]),
}


def main():
    ens2sym = {e: s.upper() for s, e in pickle.load(open(NAME_ID, "rb")).items()}
    g2g = {k.upper(): set(v) for k, v in pickle.load(open(G2G, "rb")).items() if isinstance(v, (set, list, tuple))}
    rng = np.random.default_rng(SEED)

    z0 = np.load(os.path.join(RES, "ctx_maxtoki_L04.npz"), allow_pickle=True)
    genes = z0["genes"].astype(str); ctxs = z0["contexts"].astype(str)
    syms = [ens2sym.get(g) for g in genes]
    tokmap = json.load(open(f"{CP.MSETUP}/token_dictionary.json")); ens2tid = {k: int(v) for k, v in tokmap.items()}
    tids = np.array([ens2tid.get(g, -1) for g in genes])

    print("[1/3] per-context mean rank (abundance control)", flush=True)
    MR = CP.mean_ranks(set(ctxs))
    rank = np.full((len(ctxs), len(genes)), np.nan)
    for ci, c in enumerate(ctxs):
        for gi, t in enumerate(tids):
            if t in MR.get(c, {}):
                rank[ci, gi] = MR[c][t]
    print("[2/3] co-expression matrix", flush=True)
    C = coexpr_matrix(list(genes))

    out = {"taps": {}}
    for tap in TAPS:
        z = np.load(os.path.join(RES, f"ctx_maxtoki_L{tap:02d}.npz"), allow_pickle=True)
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
        muse = full[:, use]; Ruse = rank[:, use]

        def proj_resid(vec):
            """per-partition projection onto vec, rank-residualised; returns (2, nctx, |use|) with nan mask."""
            q = np.tensordot(Mz[:, :, use], vec, axes=([3], [0]))
            qm = np.where(muse[None], q, np.nan)
            fin = np.isfinite(qm[0]) & np.isfinite(Ruse)
            r = Ruse[fin]; A = np.column_stack([np.ones_like(r), r, r ** 2])
            for p in range(2):
                yv = qm[p][fin]; qm[p][fin] = yv - A @ np.linalg.lstsq(A, yv, rcond=None)[0]
            return qm

        def beta(vec):
            """cross-validated non-additivity: predictor f*h from one partition, interaction from the other."""
            qm = proj_resid(vec)
            F, H, E, P = {}, {}, {}, {}
            for p in range(2):
                Q = qm[p]
                grand = np.nanmean(Q)
                F[p] = np.nanmean(Q, 0)                      # gene loading  (|use|,)
                H[p] = np.nanmean(Q, 1)                      # context loading (nctx,)
                E[p] = Q - F[p][None, :] - H[p][:, None] + grand
                P[p] = np.outer(H[p], F[p])                  # predictor f*h  (nctx, |use|)
            out = []
            for a, b in [(0, 1), (1, 0)]:
                sel = np.isfinite(E[a]) & np.isfinite(P[b])
                x, y = P[b][sel], E[a][sel]
                if x.std() < 1e-9:
                    continue
                out.append(float(np.corrcoef(x, y)[0, 1]))
            return float(np.mean(out)) if out else 0.0

        def axis(ia, ib):
            u = a_space[ia].mean(0) - a_space[ib].mean(0); return u / (np.linalg.norm(u) + 1e-9)

        pool = np.where(gene_ok)[0]

        def mixed(seed, n, k, base):
            """the k genes nearest/most co-expressed to `seed` plus (n-k) random genes, all from `base`"""
            top = base[np.argsort(-C[seed, base])][:k]
            rest = np.setdiff1d(base, top)
            return np.concatenate([top, rng.choice(rest, n - k, replace=False)])

        def matched(n, target, excl=None):
            """size-n set, disjoint from `excl`, with SCORE within TOL of target; returns (genes, score, below_random)
            or None. If even a purely random set (k=0) scores above the target, a random set is returned and flagged
            (the pole is less coherent/tight than random genes; power rises with the score, so this is conservative)."""
            base = pool if excl is None else np.setdiff1d(pool, excl)
            for _ in range(MAX_TRIES):
                sd_ = rng.choice(base); lo, hi = 0, n; best = None
                for _ in range(12):
                    k = (lo + hi) // 2; g = mixed(sd_, n, k, base); t = coherence(C, list(g))
                    if best is None or abs(t - target) < abs(best[1] - target):
                        best = (g, t)
                    if t < target: lo = k + 1
                    else: hi = k
                    if lo > hi: break
                if abs(best[1] - target) <= TOL * max(abs(target), 1e-3):
                    return best[0], best[1], False
                if hi == 0:
                    g = rng.choice(base, n, replace=False); return g, coherence(C, list(g)), True
            return None

        def draw_pair(nA, tA, nB, tB):
            ga = matched(nA, tA)
            if ga is None: return None
            gb = matched(nB, tB, excl=ga[0])
            return None if gb is None else (ga, gb)

        out["taps"][f"L{tap:02d}"] = {"axes": {}}
        print(f"\n=== layer {tap}: beta vs {N_NULL} size- and co-expression-matched null axes per functional axis ===",
              flush=True)
        for name, (Ag, Bg) in AXES.items():
            ia = [i for i, s in enumerate(syms) if gene_ok[i] and s in g2g and g2g[s] & set(Ag)]
            ib = [i for i, s in enumerate(syms) if gene_ok[i] and s in g2g and g2g[s] & set(Bg)]
            both = set(ia) & set(ib); ia = [i for i in ia if i not in both]; ib = [i for i in ib if i not in both]
            bt = beta(axis(ia, ib)); tA, tB = coherence(C, ia), coherence(C, ib)
            nb_, n_fail = [], 0
            for _ in range(N_NULL):
                pr = draw_pair(len(ia), tA, len(ib), tB)
                if pr is None: n_fail += 1; continue
                nb_.append(beta(axis(pr[0][0], pr[1][0])))
            nb_ = np.array(nb_)
            if len(nb_) < 50:
                raise SystemExit(f"{name}: only {len(nb_)} matched null axes ({n_fail} failures)")
            n_ge = int((nb_ >= bt).sum()); p_m = (n_ge + 1) / (len(nb_) + 1)
            print(f"   {name:<28} beta={bt:+.3f}  matched null {nb_.mean():+.3f} (sd {nb_.std():.3f}, 95th "
                  f"{np.percentile(nb_, 95):+.3f})  {n_ge}/{len(nb_)} >= functional  p = {p_m:.4f}", flush=True)
            out["taps"][f"L{tap:02d}"]["axes"][name] = dict(
                beta=bt, nA=len(ia), nB=len(ib), coherence_A=tA, coherence_B=tB, null_n=int(len(nb_)),
                null_mean=float(nb_.mean()), null_sd=float(nb_.std()), null_p95=float(np.percentile(nb_, 95)),
                n_ge=n_ge, p=p_m, null_failures=n_fail, null_beta=nb_.tolist())

    allc = [(t, n, d) for t, tv in out["taps"].items() for n, d in tv["axes"].items()]
    n_tests = len(allc); srt = sorted(allc, key=lambda x: x[2]["p"]); sig = []
    for i_, x in enumerate(srt):                      # Holm step-down over all axes x layers, positive beta only
        if x[2]["p"] * (n_tests - i_) < 0.05 and x[2]["beta"] > 0: sig.append(x)
        else: break
    t, n, d = (sig or srt)[0]
    out["holm_n_tests"] = n_tests; out["holm_significant"] = [f"{a}/{b}" for a, b, _ in sig]
    out["verdict"] = (
        f"smallest empirical p: {t}/{n} beta={d['beta']:+.3f} vs matched null {d['null_mean']:+.3f} "
        f"(sd {d['null_sd']:.3f}), p={d['p']:.4f}; Holm-corrected over {n_tests} tests: " +
        (f"SIGNAL in {len(sig)} test(s)." if sig else "NO SIGNAL beyond size- and co-expression-matched null axes."))
    print(f"\nVERDICT: {out['verdict']}")
    json.dump(out, open(os.path.join(RES, "ctx_directional_probe.json"), "w"), indent=1)
    print("[done] -> results/ctx_directional_probe.json")


if __name__ == "__main__":
    main()

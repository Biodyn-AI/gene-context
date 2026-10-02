"""IS THE FUNCTIONAL-AXIS MODULATION JUST REPRESENTATION-SPACE COMPACTNESS? (the reviewer's other null)

`ctx_coexpr_null.py` showed the functional-axis context modulation exceeds a CO-EXPRESSION-coherence-matched
null. A reviewer will ask for the other obvious match variable: maybe functional gene sets simply CLUSTER TIGHTLY
in the model's representation space, and ANY tightly-clustered set gives a strong axis regardless of function.

Design (rewritten 1 Oct 2026, mirroring ctx_coexpr_null_v3): the match variable is REPRESENTATION-SPACE TIGHTNESS,
the mean pairwise cosine of a gene set's members in the gene-main-effect space a(g). For each functional pole we draw
null poles of the SAME SIZE whose tightness is within 10% of the pole's own (the k nearest a(g)-neighbours of a random
seed gene plus random genes, k by bisection), and report the functional axis's power as an empirical percentile
among 300 such null axes. The earlier version compared the functional axes against 400-gene null poles placed on a
fitted power-vs-tightness curve, i.e. NOT size-matched (the same mismatch that invalidated ctx_coexpr_null.py).

This is the strict companion to the co-expression control: co-expression matches the DATA structure of the gene
set, tightness matches its REPRESENTATION structure. Surviving both is the strong claim.

Power is the same rank-controlled reproducible interaction power used throughout.
Out: results/ctx_tightness_null.json
"""
import os, sys, json, pickle, warnings; warnings.filterwarnings("ignore")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ctx_position_confound as CP
RES = os.path.join(HERE, "results")
NAME_ID = "/Volumes/Crucial X6/MacBook/Code/neuro-mechinterp/models/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl"
G2G = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
TAPS = [4, 8]
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


def tightness(S, idx):
    if len(idx) < 2:
        return 0.0
    sub = S[np.ix_(idx, idx)]; iu = np.triu_indices(len(idx), 1)
    return float(sub[iu].mean())


def main():
    ens2sym = {e: s.upper() for s, e in pickle.load(open(NAME_ID, "rb")).items()}
    g2g = {k.upper(): set(v) for k, v in pickle.load(open(G2G, "rb")).items() if isinstance(v, (set, list, tuple))}
    rng = np.random.default_rng(SEED)

    z0 = np.load(os.path.join(RES, "ctx_maxtoki_L04.npz"), allow_pickle=True)
    ctxs = z0["contexts"].astype(str); genes = z0["genes"].astype(str)
    syms = [ens2sym.get(g) for g in genes]
    tokmap = json.load(open(f"{CP.MSETUP}/token_dictionary.json")); ens2tid = {k: int(v) for k, v in tokmap.items()}
    tids = np.array([ens2tid.get(g, -1) for g in genes])

    print("[1/2] per-context mean rank (abundance control)", flush=True)
    MR = CP.mean_ranks(set(ctxs))
    rank = np.full((len(ctxs), len(genes)), np.nan)
    for ci, c in enumerate(ctxs):
        d = MR.get(c, {})
        for gi, t in enumerate(tids):
            if t in d:
                rank[ci, gi] = d[t]

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
        Ruse = rank[:, use]; muse = full[:, use]; Msub = Mz[:, :, use]

        # tightness similarity matrix S = cosine in a(g) space (representation compactness)
        A = np.nan_to_num(a_space)
        An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-9)
        pool = np.where(gene_ok)[0]

        def power(vec):
            q = np.tensordot(Msub, vec, axes=([3], [0]))
            qm = np.where(muse[None], q, np.nan)
            fin = np.isfinite(qm[0]) & np.isfinite(Ruse)
            r = Ruse[fin]; B = np.column_stack([np.ones_like(r), r, r ** 2])
            for p in range(2):
                yv = qm[p][fin]; qm[p][fin] = yv - B @ np.linalg.lstsq(B, yv, rcond=None)[0]
            I = qm - np.nanmean(qm, 1, keepdims=True) - np.nanmean(qm, 2, keepdims=True) + np.nanmean(qm, (1, 2), keepdims=True)
            sel = np.isfinite(I[0]) & np.isfinite(I[1])
            return float(np.nanmean(I[0][sel] * I[1][sel]))

        def axis(ia, ib):
            u = a_space[ia].mean(0) - a_space[ib].mean(0); return u / (np.linalg.norm(u) + 1e-9)

        # cosine Gram matrix in a(g) space, computed once per layer; tightness of a set = mean off-diagonal entry
        Sg = An @ An.T

        def tset_tightness(idx):
            idx = np.asarray(idx); n = len(idx)
            if n < 2:
                return 0.0
            sub = Sg[np.ix_(idx, idx)]
            return float((sub.sum() - np.trace(sub)) / (n * (n - 1)))

        def mixed(seed, n, k, base):
            """the k genes nearest/most co-expressed to `seed` plus (n-k) random genes, all from `base`"""
            top = base[np.argsort(-Sg[seed, base])][:k]
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
                    k = (lo + hi) // 2; g = mixed(sd_, n, k, base); t = tset_tightness(g)
                    if best is None or abs(t - target) < abs(best[1] - target):
                        best = (g, t)
                    if t < target: lo = k + 1
                    else: hi = k
                    if lo > hi: break
                if abs(best[1] - target) <= TOL * max(abs(target), 1e-3):
                    return best[0], best[1], False
                if hi == 0:
                    g = rng.choice(base, n, replace=False); return g, tset_tightness(g), True
            return None

        def draw_pair(nA, tA, nB, tB):
            ga = matched(nA, tA)
            if ga is None: return None
            gb = matched(nB, tB, excl=ga[0])
            return None if gb is None else (ga, gb)

        out["taps"][f"L{tap:02d}"] = {"axes": {}}
        print(f"\n=== layer {tap}: size- and tightness-matched null ({N_NULL} axes per functional axis) ===", flush=True)
        for name, (Ag, Bg) in AXES.items():
            ia = [i for i, s in enumerate(syms) if gene_ok[i] and s in g2g and g2g[s] & set(Ag)]
            ib = [i for i, s in enumerate(syms) if gene_ok[i] and s in g2g and g2g[s] & set(Bg)]
            both = set(ia) & set(ib); ia = [i for i in ia if i not in both]; ib = [i for i in ib if i not in both]
            tA, tB = tset_tightness(ia), tset_tightness(ib)
            pw = power(axis(ia, ib))
            nul, cA, cB, n_fail, n_below = [], [], [], 0, [0, 0]
            for _ in range(N_NULL):
                pr = draw_pair(len(ia), tA, len(ib), tB)
                if pr is None:
                    n_fail += 1; continue
                ga, gb = pr
                nul.append(power(axis(ga[0], gb[0]))); cA.append(ga[1]); cB.append(gb[1])
                n_below[0] += int(ga[2]); n_below[1] += int(gb[2])
            nul = np.array(nul)
            if len(nul) < 50:
                raise SystemExit(f"{name}: only {len(nul)} matched null axes ({n_fail} failures)")
            n_ge = int((nul >= pw).sum()); p_m = (n_ge + 1) / (len(nul) + 1)
            rand_t = float(np.mean([tset_tightness(rng.choice(pool, len(ia), replace=False)) for _ in range(20)]))
            print(f"   {name:<28} poles {len(ia)}/{len(ib)} tightness {tA:+.4f}/{tB:+.4f} (random {rand_t:+.4f})  "
                  f"power {pw:.3f}  matched null mean {nul.mean():.3f} 95th {np.percentile(nul, 95):.3f}  "
                  f"{n_ge}/{len(nul)} >= functional  p = {p_m:.4f}  (poles below random {n_below})", flush=True)
            out["taps"][f"L{tap:02d}"]["axes"][name] = dict(
                power=pw, nA=len(ia), nB=len(ib), tightness_A=tA, tightness_B=tB, random_set_tightness=rand_t,
                matched_n=int(len(nul)), matched_tightness_A=float(np.mean(cA)), matched_tightness_B=float(np.mean(cB)),
                matched_mean=float(nul.mean()), matched_p95=float(np.percentile(nul, 95)), matched_n_ge=n_ge,
                matched_p=p_m, matched_failures=n_fail, matched_poles_below_random=n_below, null_matched=nul.tolist())

    summ = {f"{t}/{n}": d["matched_p"] for t, tv in out["taps"].items() for n, d in tv["axes"].items()}
    out["verdict"] = ("Empirical p of each functional axis against size- and tightness-matched null axes: " +
                      ", ".join(f"{k} p={v:.3f}" for k, v in summ.items()) +
                      ". (Replaces the earlier 400-gene, not size-matched tightness curve, 1 Oct 2026.)")
    print(f"\nVERDICT: {out['verdict']}")
    json.dump(out, open(os.path.join(RES, "ctx_tightness_null.json"), "w"), indent=1)
    print("[done] -> results/ctx_tightness_null.json")


if __name__ == "__main__":
    main()

"""CO-EXPRESSION NULL MATCHED ON SIZE *AND* AVERAGE CO-EXPRESSION (completes ctx_coexpr_null_v2.py).

WHY THIS EXISTS. ctx_coexpr_null_v2.py fixed a set-SIZE mismatch, but its co-expression null (b) builds each null
pole from the n genes MOST co-expressed with a random seed gene. Those modules are NOT matched to the functional
poles on co-expression: they are 3-10x more co-expressed (mean pairwise expression correlation, "coherence") than
the functional poles, whose coherence is about that of random genes (nuclear pole 0.022 vs random 0.023). So
null (b) is a STRONG null (strongly co-expressed modules of the same size), not a matched one. The paper's
Methods described it as matched on both size and coherence; this script runs the test that description names.

THREE OUTPUTS, per functional axis (same axes, data, rank-controlled power and seed as v2):
  (1) MATCHED NULL. Each null pole has the functional pole's SIZE and a coherence within TOL (10%) of the
      functional pole's coherence: the top-k genes most co-expressed with a random seed gene plus (n-k) random
      genes, with k found by bisection. Empirical one-sided p = fraction of null axes with power >= functional.
  (2) FIXED-SIZE CURVE. Power for null poles with a fraction f in {0..1} of seed-co-expressed genes (25 draws
      each): how power grows with co-expression at the functional poles' sizes.
  (3) COHERENCE OF THE v2 STRONG MODULES (100 fresh draws of the v2 null-(b) construction, separate RNG so (1)
      and (2) are unchanged), to document how much more co-expressed they are than the functional poles.

Null (b) p-values themselves stay in results/ctx_coexpr_null_v2.json (not recomputed here).
Needs the transformers environment (rank control tokenises the panels): ../../.venv_state/bin/python
Out: results/ctx_coexpr_null_v3.json
"""
import os, sys, json, pickle, warnings; warnings.filterwarnings("ignore")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); RES = os.path.join(HERE, "results")
import ctx_position_confound as CP
from ctx_coexpr_null import coexpr_matrix, coherence
from ctx_coexpr_null_v2 import NAME_ID, G2G, AXES, TAP
N_NULL, SEED, TOL, MAX_TRIES = 300, 0, 0.10, 40
FRACS, N_CURVE, N_STRONG = [0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0], 25, 100


def main():
    ens2sym = {e: s.upper() for s, e in pickle.load(open(NAME_ID, "rb")).items()}
    g2g = {k.upper(): set(v) for k, v in pickle.load(open(G2G, "rb")).items() if isinstance(v, (set, list, tuple))}
    rng = np.random.default_rng(SEED); rng_strong = np.random.default_rng(SEED + 1)
    z = np.load(os.path.join(RES, f"ctx_maxtoki_L{TAP:02d}.npz"), allow_pickle=True)
    M, counts, cap = z["M"].astype(np.float32), z["counts"], int(z["cap"])
    genes = z["genes"].astype(str); ctxs = z["contexts"].astype(str); syms = [ens2sym.get(g) for g in genes]
    full = (counts == cap).all(0); d = M.shape[-1]
    flat = M[:, full]; mu = flat.reshape(-1, d).mean(0); sd = flat.reshape(-1, d).std(0) + 1e-6; Mz = (M - mu) / sd
    a = np.full((len(genes), d), np.nan, np.float32)
    for gi in range(len(genes)):
        cs = np.where(full[:, gi])[0]
        if len(cs): a[gi] = Mz[:, cs, gi].mean((0, 1))
    ok = np.isfinite(a[:, 0]); use = np.where(full.sum(0) >= 9)[0]
    tokmap = json.load(open(f"{CP.MSETUP}/token_dictionary.json")); ens2tid = {k: int(v) for k, v in tokmap.items()}
    tids = np.array([ens2tid.get(g, -1) for g in genes]); MR = CP.mean_ranks(set(ctxs))
    rank = np.full((len(ctxs), len(genes)), np.nan)
    for ci, c in enumerate(ctxs):
        dd = MR.get(c, {})
        for gi, t in enumerate(tids):
            if t in dd: rank[ci, gi] = dd[t]
    Ruse, muse, Msub = rank[:, use], full[:, use], Mz[:, :, use]
    del M, flat

    def power(vec):                                   # identical to ctx_coexpr_null_v2.power
        q = np.tensordot(Msub, vec, axes=([3], [0])); qm = np.where(muse[None], q, np.nan)
        fin = np.isfinite(qm[0]) & np.isfinite(Ruse); r = Ruse[fin]; A = np.column_stack([np.ones_like(r), r, r ** 2])
        for p in range(2):
            yv = qm[p][fin]; qm[p][fin] = yv - A @ np.linalg.lstsq(A, yv, rcond=None)[0]
        I = qm - np.nanmean(qm, 1, keepdims=True) - np.nanmean(qm, 2, keepdims=True) + np.nanmean(qm, (1, 2), keepdims=True)
        sel = np.isfinite(I[0]) & np.isfinite(I[1]); return float(np.nanmean(I[0][sel] * I[1][sel]))

    def axis(ia, ib):
        u = a[ia].mean(0) - a[ib].mean(0); return u / (np.linalg.norm(u) + 1e-9)

    C = coexpr_matrix(list(genes)); pool = np.where(ok)[0]

    def mixed(seed, n, k, base):
        """the k genes most co-expressed with `seed` plus (n-k) random genes, all drawn from `base`"""
        top = base[np.argsort(-C[seed, base])][:k]
        rest = np.setdiff1d(base, top); return np.concatenate([top, rng.choice(rest, n - k, False)])

    def matched(n, target, excl=None):
        """size-n set with coherence within TOL of target, disjoint from `excl` (the other null pole). Returns
        (genes, coherence, below_random) or None. If even a purely random set (k=0) is MORE coherent than the target
        (a pole less co-expressed than random genes), a random size-matched set is returned and flagged: power
        rises with coherence, so this errs on the conservative side."""
        base = pool if excl is None else np.setdiff1d(pool, excl)
        for _ in range(MAX_TRIES):
            sd_ = rng.choice(base); lo, hi = 0, n; best = None
            for _ in range(12):
                k = (lo + hi) // 2; g = mixed(sd_, n, k, base); c = coherence(C, list(g))
                if best is None or abs(c - target) < abs(best[1] - target): best = (g, c)
                if c < target: lo = k + 1
                else: hi = k
                if lo > hi: break
            if abs(best[1] - target) <= TOL * max(abs(target), 1e-3): return best[0], best[1], False
            if hi == 0:
                g = rng.choice(base, n, False); return g, coherence(C, list(g)), True
        return None

    out = {"tap": TAP, "n_null": N_NULL, "tol": TOL, "axes": {}}
    for name, (Ag, Bg) in AXES.items():
        ia = [i for i, s in enumerate(syms) if ok[i] and s in g2g and g2g[s] & set(Ag)]
        ib = [i for i, s in enumerate(syms) if ok[i] and s in g2g and g2g[s] & set(Bg)]
        both = set(ia) & set(ib); ia = [i for i in ia if i not in both]; ib = [i for i in ib if i not in both]
        nA, nB = len(ia), len(ib); pf = power(axis(ia, ib)); tA, tB = coherence(C, ia), coherence(C, ib)
        nul, cA, cB, n_fail, n_below = [], [], [], 0, [0, 0]
        for _ in range(N_NULL):
            ga = matched(nA, tA)
            if ga is None: n_fail += 1; continue
            gb = matched(nB, tB, excl=ga[0])                  # disjoint poles, like the functional ones
            if gb is None: n_fail += 1; continue
            nul.append(power(axis(ga[0], gb[0]))); cA.append(ga[1]); cB.append(gb[1])
            n_below[0] += int(ga[2]); n_below[1] += int(gb[2])
        nul = np.array(nul)
        if len(nul) < 50:
            raise SystemExit(f"{name}: only {len(nul)} matched null axes could be built ({n_fail} failures)")
        n_ge = int((nul >= pf).sum()); p_m = (n_ge + 1) / (len(nul) + 1)    # Monte Carlo p, never exactly 0
        curve = []
        for f in FRACS:
            for _ in range(N_CURVE):
                ga = mixed(rng.choice(pool), nA, int(f * nA), pool)
                bb = np.setdiff1d(pool, ga); gb = mixed(rng.choice(bb), nB, int(f * nB), bb)
                curve.append(dict(f=f, cohA=coherence(C, list(ga)), cohB=coherence(C, list(gb)),
                                  power=power(axis(ga, gb))))
        sA, sB, rA = [], [], []                       # v2 null-(b) construction: n genes most co-expressed w/ seed
        for _ in range(N_STRONG):
            s1, s2 = rng_strong.choice(pool), rng_strong.choice(pool)
            sA.append(coherence(C, list(pool[np.argsort(-C[s1, pool])][:nA])))
            sB.append(coherence(C, list(pool[np.argsort(-C[s2, pool])][:nB])))
            rA.append(coherence(C, list(rng_strong.choice(pool, nA, False))))
        rec = dict(power=pf, nA=nA, nB=nB, coh_A=tA, coh_B=tB,
                   matched_n=int(len(nul)), matched_cohA_mean=float(np.mean(cA)), matched_cohB_mean=float(np.mean(cB)),
                   matched_mean=float(nul.mean()), matched_p95=float(np.percentile(nul, 95)),
                   matched_n_ge=n_ge, matched_p=float(p_m), matched_failures=n_fail,
                   matched_poles_below_random=n_below,
                   strong_module_cohA_mean=float(np.mean(sA)), strong_module_cohB_mean=float(np.mean(sB)),
                   random_set_cohA_mean=float(np.mean(rA)),
                   null_matched=nul.tolist(), curve=curve)
        out["axes"][name] = rec
        print(f"\n=== {name} poles {nA}/{nB}  power {pf:.3f}  coherence {tA:.4f}/{tB:.4f}", flush=True)
        print(f"  matched null (n={len(nul)}): coherence {np.mean(cA):.4f}/{np.mean(cB):.4f}  mean power {nul.mean():.3f}  "
              f"95th {np.percentile(nul, 95):.3f}  {n_ge}/{len(nul)} >= functional  p = {p_m:.4f}  (poles below random: {n_below})")
        print(f"  v2 strong modules: coherence {np.mean(sA):.4f}/{np.mean(sB):.4f} "
              f"(x{np.mean(sA)/tA:.1f}/x{np.mean(sB)/tB:.1f} the functional poles); random sets {np.mean(rA):.4f}")
        for f in FRACS:
            rows = [c for c in curve if c["f"] == f]
            print(f"  curve f={f:<4} coherence {np.mean([r['cohA'] for r in rows]):.4f}/{np.mean([r['cohB'] for r in rows]):.4f}"
                  f"  power mean {np.mean([r['power'] for r in rows]):.3f}  max {np.max([r['power'] for r in rows]):.3f}")
    v2p = os.path.join(RES, "ctx_coexpr_null_v2.json")
    if os.path.getmtime(v2p) < os.path.getmtime(os.path.join(RES, f"ctx_maxtoki_L{TAP:02d}.npz")):
        raise SystemExit("ctx_coexpr_null_v2.json is older than the current extraction (stale); re-run ctx_coexpr_null_v2.py")
    v2 = json.load(open(v2p))["axes"]
    out["v2_source_mtime"] = os.path.getmtime(v2p)
    out["summary"] = {k: dict(power=v["power"], p_matched=v["matched_p"], p_strong_modules_v2=v2[k]["indep_coexpr_p"],
                              p_random_v2=v2[k]["random_p"]) for k, v in out["axes"].items()}
    ratios = [out["axes"][k][s_] / out["axes"][k][c_] for k in out["axes"]
              for s_, c_ in (("strong_module_cohA_mean", "coh_A"), ("strong_module_cohB_mean", "coh_B"))
              if out["axes"][k][c_] > 0]
    out["strong_over_functional_coherence_range"] = [float(min(ratios)), float(max(ratios))] if ratios else None
    n_sig = sum(v["p_matched"] < 0.05 for v in out["summary"].values()); n_ax = len(out["summary"])
    n_edge = sum(v["p_strong_modules_v2"] < 0.05 for v in out["summary"].values())
    out["verdict"] = (
        f"Matched null (size + average co-expression): {n_sig}/{n_ax} functional axes exceed it at p<0.05 ("
        + ", ".join(f"{k} p={v['p_matched']:.4f}" for k, v in out["summary"].items())
        + f"). Strong null (v2 null b, modules of the most co-expressed genes; "
        + (f"{ratios and min(ratios):.1f}-{ratios and max(ratios):.1f}x the functional poles' coherence" if ratios else "n/a")
        + f"): {n_edge}/{n_ax} axes exceed it at p<0.05 ("
        + ", ".join(f"{k} p={v['p_strong_modules_v2']:.3f}" for k, v in out["summary"].items()) + ").")
    print(f"\nVERDICT: {out['verdict']}")
    json.dump(out, open(os.path.join(RES, "ctx_coexpr_null_v3.json"), "w"), indent=1)
    print("[done] -> results/ctx_coexpr_null_v3.json")


if __name__ == "__main__":
    main()

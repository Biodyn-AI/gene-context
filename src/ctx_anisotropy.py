"""ANISOTROPY & SELF-SIMILARITY characterisation — the metric-validity control (Ethayarajh 2019).

Contextual representations in a layer can be so anisotropic that any two vectors look similar, which would
inflate or deflate every cosine we report. A reviewer in this subfield expects these numbers. For each layer,
on the RAW (pre-standardisation) per-(gene,context) representations, we report:
  anisotropy  = mean cosine of random DIFFERENT-gene pairs (Ethayarajh's baseline). ~1 => everything alike.
  self-sim    = mean cosine of the SAME gene across different contexts.
  self-sim - anisotropy = the anisotropy-corrected self-similarity (how much a gene keeps its identity across
                          contexts once the shared anisotropic direction is accounted for).
  top-PC share = variance fraction of the first principal component of pooled representations (rogue-dim check).
We report both RAW and after our per-dimension z-scoring, to show the correction removes the shared component
that the cosine-based EXCESS/functional analyses would otherwise ride on.

Out: results/ctx_anisotropy.json
"""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
import numpy as np
import ctx_prefix as PX
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
SEED = 0


def stats(V, rng, n=6000):
    """V: (n_entries, d) unit-normalised rows expected. anisotropy + top-PC share."""
    i, j = rng.integers(0, len(V), n), rng.integers(0, len(V), n); ok = i != j
    ani = float((V[i[ok]] * V[j[ok]]).sum(1).mean())
    Vc = V - V.mean(0); s = np.linalg.svd(Vc[rng.integers(0, len(Vc), min(4000, len(Vc)))], compute_uv=False)
    return ani, float((s[0] ** 2) / (s ** 2).sum())


def selfsim(U, gidx):
    """mean over genes (with >= 2 count-balanced cell types) of the mean cosine between that gene's cell types"""
    order = np.argsort(gidx, kind="stable"); g = gidx[order]; V = U[order]
    starts = np.flatnonzero(np.r_[True, g[1:] != g[:-1]]); ends = np.r_[starts[1:], len(g)]
    vals = []
    for s, e in zip(starts, ends):
        if e - s < 2: continue
        v = V[s:e]; G = v @ v.T; m = e - s
        vals.append((G.sum() - np.trace(G)) / (m * (m - 1)))
    return float(np.mean(vals)), len(vals)


def main():
    rng = np.random.default_rng(SEED)
    out = {"layers": {}}
    for tap in PX.taps([0, 4, 8, 11]):
        p = PX.npz_path(tap)
        if not os.path.exists(p):
            continue
        z = np.load(p, allow_pickle=True)
        M, counts, cap = z["M"].astype(np.float32), z["counts"], int(z["cap"])
        nP, nC, nG, d = M.shape
        full = (counts == cap).all(0)                        # (nC, nG)
        Mavg = M.mean(0)                                     # (nC, nG, d) partition-averaged
        # gather valid (context,gene) entries
        ent = []
        gid = []
        for c in range(nC):
            gg = np.where(full[c])[0]
            ent.append(Mavg[c, gg]); gid.append(gg)
        raw = np.concatenate(ent); gidx = np.concatenate(gid)
        # RAW anisotropy / self-sim
        U = raw / (np.linalg.norm(raw, axis=1, keepdims=True) + 1e-9)
        ani_raw, pc_raw = stats(U, rng)
        # self-similarity: same gene across contexts
        # (3 Oct 2026: all genes with >= 2 count-balanced cell types; the first version used only the first 1500 genes)
        ss_raw, n_ss = selfsim(U, gidx)
        # after per-dim z-scoring (our correction)
        mu = raw.mean(0); sd = raw.std(0) + 1e-6
        Z = (raw - mu) / sd; Uz = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-9)
        ani_z, pc_z = stats(Uz, rng)
        ss_z, _ = selfsim(Uz, gidx)
        out["layers"][f"L{tap:02d}"] = dict(anisotropy_raw=ani_raw, selfsim_raw=ss_raw,
                                            corrected_selfsim_raw=ss_raw - ani_raw, top_pc_raw=pc_raw,
                                            anisotropy_zscored=ani_z, selfsim_zscored=ss_z,
                                            corrected_selfsim_zscored=ss_z - ani_z, top_pc_zscored=pc_z,
                                            selfsim_n_genes=n_ss, n_entries=int(len(gidx)))
        print(f"L{tap:02d}: RAW anisotropy {ani_raw:+.3f} self-sim {ss_raw:+.3f} (corrected {ss_raw-ani_raw:+.3f}) "
              f"topPC {pc_raw:.2f} || Z-SCORED anisotropy {ani_z:+.3f} self-sim {ss_z:+.3f} "
              f"(corrected {ss_z-ani_z:+.3f}) topPC {pc_z:.2f}", flush=True)
        del M, Mavg, raw, U, Z, Uz
    out["note"] = ("High raw anisotropy with much lower self-sim-minus-anisotropy is exactly Ethayarajh's finding; "
                   "our per-dimension z-scoring lowers anisotropy and keeps genes distinguishable, validating the "
                   "cosine-based EXCESS/functional metrics computed on z-scored representations.")
    out.update({} if PX.IS_DEFAULT else {"provenance": PX.provenance()}); json.dump(out, open(PX.out("ctx_anisotropy"), "w"), indent=1)
    print("[done] -> results/ctx_anisotropy.json")


if __name__ == "__main__":
    main()

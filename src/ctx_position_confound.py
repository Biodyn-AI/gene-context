"""IS THE GENE-SPECIFIC CONTEXT RESPONSE JUST RANK-POSITION ENCODING?

`ctx_polysemy.py` found a large, reproducible gene-specific context response (EXCESS +0.758 at L4, +0.665 at
L8, +0.627 at L11; exactly 0.000 at L0, the context-free embedding layer, which validates the pipeline).

BEFORE THAT COUNTS AS BIOLOGY, THE OBVIOUS CONFOUND MUST DIE. MaxToki reads a cell as a rank-ordered gene list
and uses RoPE, and this project has already established at length that the model encodes genomic and sequence
POSITION. A gene ranks 5th in a macrophage and 500th in a T cell purely because it is expressed more there.
That rank change is gene-specific, context-specific, and perfectly reproducible across cell partitions --
i.e. it would produce exactly the signal we measured, with no functional content whatsoever.

THE TEST. Recompute each gene's mean rank per context from tokenisation (CPU only, no forward pass), then:

  1. CORRELATION. Does the per-gene shift magnitude ||delta(g)|| track |mean_rank(g,c2) - mean_rank(g,c1)|?
     A strong positive correlation means the effect is substantially positional.
  2. STRATIFICATION. Split genes by how much their rank moved between the two contexts and recompute EXCESS in
     each stratum. If EXCESS is large only where rank moved a lot, and collapses where rank is stable, the
     finding is position. If EXCESS survives in the rank-stable stratum, something else is going on.
  3. RESIDUALISATION. Regress the per-gene shift on rank change (linear, on ranks) and recompute directional
     agreement on the residual.

Stratum (3) is the honest headline number: gene-specific context response AT MATCHED RANK.

Out: results/ctx_position_confound.json
"""
import os, sys, json, collections, itertools, warnings; warnings.filterwarnings("ignore")
import numpy as np, h5py
import ctx_prefix as PX

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
MSETUP = "/Volumes/Crucial X6/MacBook/Code/biomi_automation/projects/maxtoki/setup"
sys.path.insert(0, MSETUP)
RES = os.path.join(HERE, "results")
TS = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/raw"
PANELS = ["tabula_sapiens_immune_subset_20000.h5ad", "tabula_sapiens_kidney.h5ad", "tabula_sapiens_lung.h5ad"]
MAX_LEN, CELLS_CTX, SEED = 1024, 1000, 0
TAPS = PX.taps([4, 8, 11])
from scipy.stats import spearmanr


def mean_ranks(ctx_names):
    """per (context, token) mean rank position, replicating the extractor's tokenisation AND its cell selection
    exactly (same encoding via ctx_tokenise; same per-context shuffle order and seed as ctx_extract_maxtoki, so the
    ranks come from the same CELLS_CTX cells per context that the representations were averaged over)."""
    import ctx_tokenise as TK
    from maxtoki_adapter import MaxTokiTokenizer
    tok = MaxTokiTokenizer(model_input_size=MAX_LEN)
    rng = np.random.default_rng(SEED)
    cells = collections.defaultdict(list)
    for p in PANELS:
        path = os.path.join(TS, p)
        if not os.path.exists(path):
            continue
        with h5py.File(path, "r") as f:
            ens = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]).astype(str)
            ens = np.array([e.split(".")[0] for e in ens])
            ctg = f["obs"]["cell_type"]
            cats = np.array([x.decode() if isinstance(x, bytes) else x for x in ctg["categories"][:]]).astype(str)
            ctypes = cats[ctg["codes"][:]]
            X = TK.count_matrix(f); TK.check_counts(X); n = int(X.attrs["shape"][0]); indptr = X["indptr"][:]
            var_idx, token_ids, medians = tok.make_var_mapping(list(ens))
            pos = np.full(len(ens), -1, np.int64); pos[var_idx] = np.arange(len(var_idx))
            for r in range(n):
                if ctypes[r] not in ctx_names:
                    continue
                s, e = int(indptr[r]), int(indptr[r + 1])
                idx, val = X["indices"][s:e], X["data"][s:e].astype(np.float32)
                keep = pos[idx] >= 0
                if not keep.any():
                    continue
                j = pos[idx[keep]]
                order = TK.rank_order(val[keep], medians[j])[: MAX_LEN - 2]
                if not len(order):
                    continue
                cells[ctypes[r]].append(token_ids[j[order]].astype(np.int64))
    acc = {}
    # Same order as the extractor: contexts by descending cell count (stable sort over first appearance). This makes
    # the result independent of PYTHONHASHSEED (iterating a set made it vary between runs) and reuses the extractor's
    # exact per-context shuffles, so the cells match.
    for c in sorted(cells.keys(), key=lambda c: -len(cells[c])):
        rng.shuffle(cells[c])
        s = collections.Counter(); n = collections.Counter()
        for toks in cells[c][:CELLS_CTX]:
            for rank, t in enumerate(toks):
                s[int(t)] += rank; n[int(t)] += 1
        acc[c] = {g: s[g] / n[g] for g in s}
    return acc


def cos_rows(A, B):
    A = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-9)
    B = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-9)
    return (A * B).sum(1)


def main():
    z0 = np.load(PX.npz_path(TAPS[0]), allow_pickle=True)
    ctxs = z0["contexts"].astype(str); genes = z0["genes"].astype(str)
    tokmap = json.load(open(f"{MSETUP}/token_dictionary.json"))
    ens2tid = {k: int(v) for k, v in tokmap.items()}
    tids = np.array([ens2tid.get(g, -1) for g in genes])

    print("[1/2] recomputing per-gene mean rank per context (tokenisation only)", flush=True)
    rank = PX.covariate_matrix(ctxs, genes)   # MaxToki mean rank (COV=rank) on the extraction's cells
    print(f"      rank table filled for {np.isfinite(rank).mean():.1%} of (context, gene) cells")

    out = {"taps": {}}
    rng = np.random.default_rng(SEED)
    for tap in TAPS:
        z = np.load(PX.npz_path(tap), allow_pickle=True)
        M, counts, cap = z["M"].astype(np.float32), z["counts"], int(z["cap"])
        full = (counts == cap).all(0)
        flat = M[:, full]
        mu = flat.reshape(-1, M.shape[-1]).mean(0); sd = flat.reshape(-1, M.shape[-1]).std(0) + 1e-6
        Mz = (M - mu) / sd

        mags, drank, S_lo, D_lo, S_hi, D_hi, S_res, D_res = [], [], [], [], [], [], [], []
        S_nl, D_nl, D_rm, S_all = [], [], [], []               # curved rank regression; rank-matched other gene
        for c1, c2 in itertools.combinations(range(len(ctxs)), 2):
            keep = full[c1] & full[c2] & np.isfinite(rank[c1]) & np.isfinite(rank[c2])
            if keep.sum() < 200:
                continue
            D0 = Mz[0, c2, keep] - Mz[0, c1, keep]; D1 = Mz[1, c2, keep] - Mz[1, c1, keep]
            d0, d1 = D0 - D0.mean(0), D1 - D1.mean(0)
            dr = np.abs(rank[c2, keep] - rank[c1, keep])
            mags.append(np.linalg.norm((d0 + d1) / 2, axis=1)); drank.append(dr)
            med = np.median(dr)
            lo, hi = dr <= med, dr > med                       # rank-stable vs rank-moving genes
            perm = rng.permutation(keep.sum())
            S_lo.append(cos_rows(d0[lo], d1[lo])); D_lo.append(cos_rows(d0[lo], d1[perm][lo]))
            S_hi.append(cos_rows(d0[hi], d1[hi])); D_hi.append(cos_rows(d0[hi], d1[perm][hi]))
            # residualise both halves on rank change (linear in rank), then re-test direction
            A = np.column_stack([np.ones(keep.sum()), rank[c1, keep], rank[c2, keep], dr])
            proj = lambda V: V - A @ np.linalg.lstsq(A, V, rcond=None)[0]
            r0, r1 = proj(d0), proj(d1)
            S_res.append(cos_rows(r0, r1)); D_res.append(cos_rows(r0, r1[perm]))
            # (added 2 Oct 2026 after review) a CURVED rank effect: squares, product and logs of the ranks
            q1, q2, qd = rank[c1, keep] / 1000, rank[c2, keep] / 1000, dr / 1000
            A2 = np.column_stack([np.ones(keep.sum()), q1, q2, qd, q1 ** 2, q2 ** 2, qd ** 2, q1 * q2,
                                  np.log1p(rank[c1, keep]), np.log1p(rank[c2, keep])])
            proj2 = lambda V: V - A2 @ np.linalg.lstsq(A2, V, rcond=None)[0]
            n0, n1 = proj2(d0), proj2(d1)
            S_nl.append(cos_rows(n0, n1)); D_nl.append(cos_rows(n0, n1[perm]))
            # RANK-MATCHED other gene: pair each gene with the gene nearest to it in (rank in c1, rank in c2).
            # If rank drove the shift, rank-matched genes would share it and this "different gene" term would rise.
            from scipy.spatial import cKDTree
            P2 = np.column_stack([rank[c1, keep], rank[c2, keep]])
            nn = cKDTree(P2).query(P2, k=2)[1][:, 1]
            D_rm.append(cos_rows(d0, d1[nn])); S_all.append(cos_rows(d0, d1))

        mags = np.concatenate(mags); drank = np.concatenate(drank)
        rho = float(spearmanr(mags, drank).statistic)
        ex = lambda S, D: float(np.concatenate(S).mean() - np.concatenate(D).mean())
        e_lo, e_hi, e_res = ex(S_lo, D_lo), ex(S_hi, D_hi), ex(S_res, D_res)
        e_nl = ex(S_nl, D_nl)
        d_rm = float(np.concatenate(D_rm).mean()); s_all = float(np.concatenate(S_all).mean())
        print(f"\n=== layer {tap} ===")
        print(f"  shift magnitude vs |rank change|        rho = {rho:+.3f}")
        print(f"  EXCESS, rank-STABLE genes (below median) : {e_lo:+.4f}")
        print(f"  EXCESS, rank-MOVING genes (above median) : {e_hi:+.4f}")
        print(f"  EXCESS after residualising on rank       : {e_res:+.4f}   <-- the honest number")
        print(f"  EXCESS after a CURVED rank regression    : {e_nl:+.4f}")
        print(f"  agreement with the RANK-MATCHED other gene: {d_rm:+.4f} (a random other gene gives about 0); "
              f"same gene {s_all:+.4f}", flush=True)
        out["taps"][f"L{tap:02d}"] = dict(rho_mag_vs_rank=rho, excess_rank_stable=e_lo,
                                          excess_rank_moving=e_hi, excess_residualised=e_res,
                                          excess_residualised_curved=e_nl, diff_rank_matched=d_rm,
                                          excess_vs_rank_matched=s_all - d_rm, same_all=s_all,
                                          median_abs_rank_change_stable=float(np.median(drank[drank <= np.median(drank)])),
                                          median_abs_rank_change_moving=float(np.median(drank[drank > np.median(drank)])))

    best = max(out["taps"].values(), key=lambda v: v["excess_residualised"])
    out["verdict"] = (
        f"residualised EXCESS {best['excess_residualised']:+.4f}. " +
        ("SURVIVES rank control — the gene-specific context response is not merely token-position encoding."
         if best["excess_residualised"] > 0.05 else
         "COLLAPSES under rank control — the apparent gene-specific context response is substantially "
         "token-rank position, which this model is already known to encode. Not biology."))
    print(f"\nVERDICT: {out['verdict']}")
    out.update({} if PX.IS_DEFAULT else {"provenance": PX.provenance()}); json.dump(out, open(PX.out("ctx_position_confound"), "w"), indent=1)
    print("[done] -> results/ctx_position_confound.json")


if __name__ == "__main__":
    main()


def build_covariate_cache(n_cells, path):
    """Per (context, Ensembl gene) covariates on EXACTLY the extraction's cells (results/ctx_cell_selection.npz,
    cells_{n_cells}): mean_rank (MaxToki rank, averaged over the cells expressing the gene -- identical to mean_ranks),
    n_expr (cells expressing it) and mean_logcp10k_all (mean log1p CP10k over all the context's cells)."""
    import ctx_tokenise as TK
    from maxtoki_adapter import MaxTokiTokenizer
    tok = MaxTokiTokenizer(model_input_size=MAX_LEN)
    sel = np.load(os.path.join(RES, "ctx_cell_selection.npz"), allow_pickle=True)
    cells = sel[f"cells_{n_cells}"]; panels = sel["panels"].astype(str); ctxs = list(sel["contexts"].astype(str))
    tokmap = json.load(open(f"{MSETUP}/token_dictionary.json")); tid2ens = {int(v): k for k, v in tokmap.items()}
    rank_sum, n_expr = collections.defaultdict(float), collections.Counter()
    expr_sum = collections.defaultdict(float); n_cells_ctx = collections.Counter()
    for pi, p in enumerate(panels):
        rows = cells[cells[:, 0] == pi]
        if not len(rows):
            continue
        with h5py.File(os.path.join(TS, p), "r") as f:
            ens = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]).astype(str)
            ens = np.array([e.split(".")[0] for e in ens])
            X = TK.count_matrix(f); TK.check_counts(X); indptr = X["indptr"][:]
            var_idx, token_ids, medians = tok.make_var_mapping(list(ens))
            pos = np.full(len(ens), -1, np.int64); pos[var_idx] = np.arange(len(var_idx))
            for _, r, ci, _ in rows:
                s, e = int(indptr[r]), int(indptr[r + 1])
                idx, val = X["indices"][s:e], X["data"][s:e].astype(np.float32)
                n_cells_ctx[ci] += 1
                lg = TK.log_cp10k(val, val.sum())
                for g, v in zip(ens[idx], lg):
                    expr_sum[(ci, g)] += float(v)
                keep = pos[idx] >= 0
                j = pos[idx[keep]]
                order = TK.rank_order(val[keep], medians[j])[: MAX_LEN - 2]
                for rank, t in enumerate(token_ids[j[order]]):
                    g = tid2ens.get(int(t))
                    rank_sum[(ci, g)] += rank; n_expr[(ci, g)] += 1
    genes = sorted({g for _, g in list(n_expr) + list(expr_sum)})
    gi = {g: i for i, g in enumerate(genes)}
    mr = np.full((len(ctxs), len(genes)), np.nan); ne = np.zeros((len(ctxs), len(genes)), np.int32)
    me = np.zeros((len(ctxs), len(genes)))
    for (ci, g), s in rank_sum.items():
        mr[ci, gi[g]] = s / n_expr[(ci, g)]; ne[ci, gi[g]] = n_expr[(ci, g)]
    for (ci, g), s in expr_sum.items():
        me[ci, gi[g]] = s / n_cells_ctx[ci]
    np.savez_compressed(path, contexts=np.array(ctxs), genes=np.array(genes), mean_rank=mr, n_expr=ne,
                        mean_logcp10k_all=me, cells_ctx=n_cells)
    print(f"[covariates] {path}: {len(ctxs)} contexts x {len(genes)} genes", flush=True)

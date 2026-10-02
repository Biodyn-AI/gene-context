"""DOES CONTEXTUALISATION SERVE THE MODEL'S OBJECTIVE? (the 'so what', and it closes the loop)

We have shown genes are contextualised, and that the contextualisation IS co-expression. The missing link is
WHY the model does it. MaxToki is autoregressive over rank-ordered genes: at each position it predicts the next
gene. If contextualisation is functional for that objective, then the genes the model contextualises MORE should
be the genes whose prediction benefits MORE from having the real (co-expressed) context present.

INTERVENTION (a clean causal measure of 'does context help predict gene g'):
  For a target gene g at rank-position p in a real cell, compare the model's log-probability of g given
    - REAL prefix:     [BOS] + the cell's own genes ranked 1..p-1     (the true co-expression context)
    - CHIMERIC prefix: [BOS] + a random OTHER cell's genes ranked 1..p-1  (a context lacking g's co-expression)
  context_benefit(g) = logP(g | real) - logP(g | chimeric), averaged over occurrences.
Reading logits only at the final prefix position = the model's own next-gene distribution; no probe, no
circularity. The gene set differs between real and chimeric, so this isolates the contribution of WHICH genes
co-occur, not rank order alone.

CONTEXTUALISATION(g): from the extraction, how much g's representation moves across cell types --
  mean over the contexts where g is count-balanced of || z(g,c) - mean_c z(g,c) ||  (z-scored L4 space).

THE TEST: partial Spearman of context_benefit(g) vs contextualisation(g), controlling for log gene frequency
(frequent genes are both better predicted and more sampled). A POSITIVE partial correlation means the model
contextualises exactly the genes whose prediction its context improves -- i.e. contextualisation is the model
encoding the co-expression it predicts from, unifying the representation phenomenon with the co-expression
ceiling and with the training objective.

Out: results/ctx_prediction_link.json
"""
import os, sys, json, pickle, warnings; warnings.filterwarnings("ignore")
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np, h5py, torch

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import steer_lib as SL
MSETUP = "/Volumes/Crucial X6/MacBook/Code/biomi_automation/projects/maxtoki/setup"
sys.path.insert(0, MSETUP)
RES = os.path.join(HERE, "results")
NAME_ID = "/Volumes/Crucial X6/MacBook/Code/neuro-mechinterp/models/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl"
TS = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/raw"
PANEL = "tabula_sapiens_immune_subset_20000.h5ad"
N_CELLS, POS_PER_CELL, MAX_LEN, SEED = int(os.environ.get("N_CELLS", 350)), 20, 512, 0
# MODEL=1b runs it on MaxToki-1B; PREFIX/TAP name the extraction used for per-gene contextualisation
MODEL = os.environ.get("MODEL", "217m"); XPREFIX = os.environ.get("PREFIX", "ctx_maxtoki"); XTAP = int(os.environ.get("TAP", 4))
OUTNAME = ("ctx_prediction_link.json" if (MODEL, XPREFIX, XTAP) == ("217m", "ctx_maxtoki", 4)
           else f"ctx_prediction_link__{MODEL}_{XPREFIX}_L{XTAP:02d}.json")
# PAIRED=1: the cells, target positions and donors depend only on the seed (not on the model or the extraction), so
# two models are scored on identical (cell, position, donor) triples; one forward per cell (MaxToki is causal, so row p
# of a forward on [BOS]+cell gives the next-gene distribution after the cell's first p genes); and a second donor arm
# that draws the donor from the SAME cell type (co-expression beyond cell identity).
PAIRED = os.environ.get("PAIRED") == "1"
if PAIRED:
    OUTNAME = f"ctx_prediction_link_paired__{MODEL}_{XPREFIX}_L{XTAP:02d}.json"
from scipy.stats import spearmanr, rankdata


def tokenise(tok, n, rng, types=None):
    with h5py.File(os.path.join(TS, PANEL), "r") as f:
        ctg = f["obs"]["cell_type"]
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in ctg["categories"][:]]).astype(str)
        ctypes = cats[ctg["codes"][:]]
        ens = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]).astype(str)
        ens = np.array([e.split(".")[0] for e in ens])
        var_idx, token_ids, medians = tok.make_var_mapping(list(ens))
        pos = np.full(len(ens), -1, np.int64); pos[var_idx] = np.arange(len(var_idx))
        import ctx_tokenise as TK
        X = TK.count_matrix(f); TK.check_counts(X); N = int(X.attrs["shape"][0]); indptr = X["indptr"][:]
        sel = np.sort(rng.choice(N, min(n, N), replace=False)); seqs = []
        for r in sel:
            s, e = int(indptr[r]), int(indptr[r + 1]); idx, val = X["indices"][s:e], X["data"][s:e].astype(np.float32)
            keep = pos[idx] >= 0
            if not keep.any():
                continue
            j = pos[idx[keep]]; order = TK.rank_order(val[keep], medians[j])
            if len(order) < 60:
                continue
            seqs.append(token_ids[j[order[: MAX_LEN - 2]]].astype(np.int64))
            if types is not None:
                types.append(ctypes[r])
    return seqs


def contextualisation_per_gene():
    """per-gene across-context movement in z-scored L4 space; returns {tid: value}."""
    z = np.load(os.path.join(RES, f"{XPREFIX}_L{XTAP:02d}.npz"), allow_pickle=True)
    M, counts, cap, genes = z["M"].astype(np.float32), z["counts"], int(z["cap"]), z["genes"].astype(str)
    full = (counts == cap).all(0); d = M.shape[-1]
    flat = M[:, full]; mu = flat.reshape(-1, d).mean(0); sd = flat.reshape(-1, d).std(0) + 1e-6
    Mz = (M - mu) / sd; Mavg = Mz.mean(0)                      # (nctx, ngene, d)
    tokmap = json.load(open(f"{MSETUP}/token_dictionary.json")); ens2tid = {k: int(v) for k, v in tokmap.items()}
    out = {}
    for gi in range(len(genes)):
        cs = np.where(full[:, gi])[0]
        if len(cs) < 4:
            continue
        V = Mavg[cs, gi]; ctr = V - V.mean(0)
        out[ens2tid.get(genes[gi], -1)] = float(np.linalg.norm(ctr, axis=1).mean())
    out.pop(-1, None)
    return out


def main():
    st = SL.Steerer(model_dir=SL.MODELS[MODEL])
    rng = np.random.default_rng(SEED)
    seqs = tokenise(st.tok, N_CELLS, rng)
    print(f"[setup] {len(seqs)} cells", flush=True)
    ctxal = contextualisation_per_gene()
    print(f"[setup] contextualisation for {len(ctxal)} genes", flush=True)

    BOS, EOS = st.tok.BOS, st.tok.EOS
    benefit = {}                    # tid -> list of logP(real) - logP(chimeric)
    freq = {}
    def logprob(prefix_tokens, target):
        ids = np.concatenate([[BOS], prefix_tokens]).astype(np.int64)
        lg = st.logits(ids)                                   # torch tensor (seq, vocab), float32 on CPU
        lp = torch.log_softmax(lg[-1], -1)
        return float(lp[int(target)])

    n_skipped = 0
    for si, s in enumerate(seqs):
        if len(s) < 20:
            continue
        positions = rng.choice(np.arange(5, len(s)), min(POS_PER_CELL, len(s) - 5), replace=False)
        for p in positions:
            g = int(s[p])
            if g not in ctxal:
                continue
            # chimeric donor (fixed 1 Oct 2026): another cell with at least p genes whose first p genes do NOT
            # include the target. Earlier versions padded short donors with the real cell's own genes and allowed
            # the target inside the donor prefix (a repeated gene, which the model scores very low), both of which
            # inflated the benefit.
            oth = None
            for _ in range(50):
                cand = int(rng.integers(0, len(seqs)))
                if cand == si or len(seqs[cand]) < p or g in set(seqs[cand][:p].tolist()):
                    continue
                oth = seqs[cand][:p]; break
            if oth is None:
                n_skipped += 1
                continue
            real_lp = logprob(s[:p], g)
            chim_lp = logprob(oth, g)
            benefit.setdefault(g, []).append(real_lp - chim_lp)
            freq[g] = freq.get(g, 0) + 1
        if si % 25 == 0:
            print(f"    cell {si}/{len(seqs)}", flush=True)

    genes = [g for g, v in benefit.items() if len(v) >= 3 and g in ctxal]
    ben = np.array([np.mean(benefit[g]) for g in genes])
    ctx = np.array([ctxal[g] for g in genes])
    frq = np.array([np.log(freq[g]) for g in genes])
    # global sanity: does context help at all?
    all_ben = np.concatenate([benefit[g] for g in genes])
    print(f"\n[result] {len(genes)} target genes; mean context benefit (logP real - chimeric) = "
          f"{all_ben.mean():+.3f} (positive => real context helps prediction)")

    def partial(a, b, c):
        ra, rb, rc = rankdata(a), rankdata(b), rankdata(c)
        A = np.column_stack([np.ones_like(rc), rc]); res = lambda v: v - A @ np.linalg.lstsq(A, v, rcond=None)[0]
        return float(spearmanr(res(ra), res(rb)).statistic)
    raw = float(spearmanr(ben, ctx).statistic)
    par = partial(ben, ctx, frq)
    rng2 = np.random.default_rng(1)
    bs = []
    for _ in range(2000):
        i = rng2.integers(0, len(genes), len(genes))
        try: bs.append(partial(ben[i], ctx[i], frq[i]))
        except Exception: pass
    lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    print(f"[result] benefit ~ contextualisation: raw rho {raw:+.3f}; partial (control log-freq) {par:+.3f} "
          f"95% CI [{lo:+.3f},{hi:+.3f}]")
    out = dict(n_genes=len(genes), mean_context_benefit=float(all_ben.mean()),
               median_context_benefit=float(np.median(all_ben)), n_comparisons=int(len(all_ben)),
               n_positions_skipped_no_clean_donor=int(n_skipped),
               raw_rho=raw, partial_rho=par, partial_ci=[lo, hi], n_cells=len(seqs),
               model=MODEL, extraction=f"{XPREFIX}_L{XTAP:02d}")
    out["verdict"] = (
        (f"CONTEXT HELPS PREDICTION (mean benefit {all_ben.mean():+.2f} nats) AND the model contextualises the "
         f"genes it helps most (partial rho {par:+.3f}, CI [{lo:+.3f},{hi:+.3f}] controlling frequency). "
         "Contextualisation serves the objective: the model encodes the co-expression it predicts from — "
         "unifying the representation phenomenon, the co-expression ceiling, and the training objective."
         if lo > 0 else
         f"Context helps prediction (mean benefit {all_ben.mean():+.2f}) but its per-gene strength does NOT track "
         f"contextualisation once frequency is controlled (partial rho {par:+.3f}, CI [{lo:+.3f},{hi:+.3f}]). "
         "Contextualisation and context-benefit are separate; do not claim the objective link."))
    print(f"\nVERDICT: {out['verdict']}")
    json.dump(out, open(os.path.join(RES, OUTNAME), "w"), indent=1)
    print(f"[done] -> results/{OUTNAME}")


def partial_spearman(a, b, c):
    ra, rb, rc = rankdata(a), rankdata(b), rankdata(c)
    A = np.column_stack([np.ones_like(rc), rc]); res = lambda v: v - A @ np.linalg.lstsq(A, v, rcond=None)[0]
    return float(spearmanr(res(ra), res(rb)).statistic)


def main_paired():
    import gc
    from maxtoki_adapter import MaxTokiTokenizer
    ctxal = contextualisation_per_gene(); gc.collect()            # big arrays freed before the model is loaded
    print(f"[setup] contextualisation for {len(ctxal)} genes from {XPREFIX}_L{XTAP:02d}", flush=True)
    tok = MaxTokiTokenizer(model_input_size=4096)
    types = []
    seqs = tokenise(tok, N_CELLS, np.random.default_rng(SEED), types)
    types = np.array(types)
    print(f"[setup] {len(seqs)} cells, {len(set(types))} cell types", flush=True)
    lens = np.array([len(x) for x in seqs])
    # ---- queries: depend only on the seed (identical for every model / extraction) ----
    Q = []                                           # (cell, position, target, donor_random, donor_sametype)
    def draw(si, p, g, pool, r):
        for _ in range(50):
            c = int(pool[r.integers(0, len(pool))])
            if c != si and lens[c] >= p and g not in set(seqs[c][:p].tolist()):
                return c
        return -1
    allc = np.arange(len(seqs))
    for si, s in enumerate(seqs):
        if len(s) < 20:
            continue
        r = np.random.default_rng([SEED, si])
        for p in r.choice(np.arange(5, len(s)), min(POS_PER_CELL, len(s) - 5), replace=False):
            g = int(s[p]); rp = np.random.default_rng([SEED, si, int(p)])
            Q.append((si, int(p), g, draw(si, p, g, allc, rp), draw(si, p, g, np.where(types == types[si])[0], rp)))
    Q = np.array(Q, np.int64)
    print(f"[setup] {len(Q)} target positions; no clean random donor {int((Q[:, 3] < 0).sum())}, "
          f"no clean same-type donor {int((Q[:, 4] < 0).sum())}", flush=True)
    st = SL.Steerer(model_dir=SL.MODELS[MODEL])
    BOS = st.tok.BOS
    real = np.full(len(Q), np.nan); rnd = np.full(len(Q), np.nan); same = np.full(len(Q), np.nan)
    for c in range(len(seqs)):
        qr = np.where(Q[:, 0] == c)[0]; qd = np.where(Q[:, 3] == c)[0]; qs = np.where(Q[:, 4] == c)[0]
        if not (len(qr) or len(qd) or len(qs)):
            continue
        rows = np.unique(np.concatenate([Q[qr, 1], Q[qd, 1], Q[qs, 1]]))
        lg = st.logits(np.concatenate([[BOS], seqs[c]]).astype(np.int64))
        lp = torch.log_softmax(lg[torch.as_tensor(rows)].float(), -1).cpu().numpy(); ri = {int(x): i for i, x in enumerate(rows)}
        for qi in qr: real[qi] = lp[ri[int(Q[qi, 1])], Q[qi, 2]]
        for qi in qd: rnd[qi] = lp[ri[int(Q[qi, 1])], Q[qi, 2]]
        for qi in qs: same[qi] = lp[ri[int(Q[qi, 1])], Q[qi, 2]]
        if c % 50 == 0:
            print(f"    cell {c}/{len(seqs)}", flush=True)
    # equivalence check against separate prefix forwards (the slow method) on 20 queries
    dmax = 0.0
    for qi in np.random.default_rng(1).choice(len(Q), 20, replace=False):
        si, p, g = Q[qi, :3]
        lg = st.logits(np.concatenate([[BOS], seqs[si][:p]]).astype(np.int64))
        dmax = max(dmax, abs(float(torch.log_softmax(lg[-1].float(), -1)[g]) - real[qi]))
    print(f"[check] one-forward vs prefix-forward log-prob: max |diff| over 20 queries = {dmax:.2e}", flush=True)

    out = dict(model=MODEL, extraction=f"{XPREFIX}_L{XTAP:02d}", n_cells=len(seqs), n_queries=int(len(Q)),
               equivalence_max_abs_diff=dmax, arms={})
    rng2 = np.random.default_rng(1)
    for arm, don, col in (("random_donor", rnd, 3), ("same_celltype_donor", same, 4)):
        ok = (Q[:, col] >= 0) & np.isfinite(don)
        ben = {}
        for qi in np.where(ok)[0]:
            ben.setdefault(int(Q[qi, 2]), []).append(real[qi] - don[qi])
        genes = [g for g, v in ben.items() if len(v) >= 3 and g in ctxal]
        b = np.array([np.mean(ben[g]) for g in genes]); cx = np.array([ctxal[g] for g in genes])
        fq = np.array([np.log(len(ben[g])) for g in genes])
        if len(genes) < 10:
            out["arms"][arm] = dict(n_genes=len(genes), note="too few genes with >= 3 comparisons"); continue
        allb = np.concatenate([ben[g] for g in genes])
        par = partial_spearman(b, cx, fq); bs = []
        for _ in range(2000):
            i = rng2.integers(0, len(genes), len(genes))
            try: bs.append(partial_spearman(b[i], cx[i], fq[i]))
            except Exception: pass
        lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
        out["arms"][arm] = dict(n_genes=len(genes), n_comparisons=int(len(allb)), mean_benefit=float(allb.mean()),
                                median_benefit=float(np.median(allb)), raw_rho=float(spearmanr(b, cx).statistic),
                                partial_rho=par, partial_ci=[lo, hi],
                                per_gene_benefit={str(g): float(v) for g, v in zip(genes, b)})
        print(f"[{arm}] {len(genes)} genes; mean benefit {allb.mean():+.3f} nats; partial rho {par:+.3f} "
              f"CI [{lo:+.3f},{hi:+.3f}]", flush=True)
    json.dump(out, open(os.path.join(RES, OUTNAME), "w"), indent=1)
    print(f"[done] -> results/{OUTNAME}")


if __name__ == "__main__":
    main_paired() if PAIRED else main()

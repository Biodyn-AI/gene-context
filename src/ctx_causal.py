"""IS THE FUNCTIONAL-CONTEXT DIRECTION CAUSALLY USED? (upgrades Level 1 from 'represented' to 'used')

Level 1 established that genes move along functional axes (nuclear/transcriptional vs surface/secreted) with cell
context. That is a REPRESENTATION fact. This asks whether the model's COMPUTATION uses that direction: if we push
some of a cell's genes along the functional axis, does the model raise its predicted probability of
functionally-matching genes at the OTHER, untouched positions?

DESIGN (the chromosome-steering playbook, per STEERING_TOOL.md):
  - Direction = the functional-context axis at layer 4 (mean raw hidden state of nuclear-pole genes minus
    surface-pole genes; the same axis Level 1 is about), injected after layer 3 (= the hidden_states[4] tap).
  - Push it into a RANDOM HALF of a cell's gene positions; READ the model's own next-gene logits at the OTHER
    half. Because read positions are disjoint from steer positions, any effect travels through ATTENTION, not
    local pass-through.
  - Readout = the MODEL'S OWN logits (no fitted probe -> non-circular): mean logit of nuclear-pole gene tokens
    (target) vs surface-pole gene tokens (control). specificity = d_target - d_control.
  - Dose-response over alpha; norm-matched RANDOM direction as the control (must stay flat).
NON-CIRCULARITY: the direction is built from layer-4 representations; the readout is the model's output logits
for GO-defined gene sets; read positions never overlap steer positions. Nothing is read in the basis it was
defined in.

A positive result -- target logits rise with alpha, control flat, random flat, at disjoint positions -- means
the functional-context direction is a channel the model's computation actually uses.

Out: results/ctx_causal.json
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
G2G = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/perturb/gene2go_all.pkl"
TS = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/raw"
PANEL = "tabula_sapiens_immune_subset_20000.h5ad"
# Axis is parameterised so the causal test can be shown for MORE than one functional direction and layer
# (a reviewer's cherry-pick worry). Defaults reproduce the headline nuclear-vs-surface @ layer 4 run.
AX_DEFS = {
    "nuc_surf":       (["GO:0005634", "GO:0000785", "GO:0003677"], ["GO:0005886", "GO:0005576", "GO:0005615"]),
    "mito_cyto":      (["GO:0005739"], ["GO:0005856"]),
    "trans_transport": (["GO:0006355", "GO:0003700"], ["GO:0006811", "GO:0038023"]),
}
AXIS = os.environ.get("AXIS", "nuc_surf")
NUC, SURF = AX_DEFS[AXIS]
SITE = int(os.environ.get("SITE", 3))   # inject after this layer == the hidden_states[SITE+1] tap
# MODEL=1b runs the same test on MaxToki-1B; PREFIX names the extraction the axis is built from (default ctx_maxtoki).
MODEL = os.environ.get("MODEL", "217m"); XPREFIX = os.environ.get("PREFIX", "ctx_maxtoki")
OUTNAME = "ctx_causal.json" if (AXIS == "nuc_surf" and SITE == 3) else f"ctx_causal_{AXIS}_L{SITE+1:02d}.json"
# ALPHA_UNIT=centred sets the push size from the residual norm AFTER removing the across-gene mean (averaged over 5
# cells), instead of the raw norm of one cell. Needed to compare models: MaxToki-1B has two huge, nearly constant hidden
# dimensions that dominate its raw norm, which would make its push ~3x larger in relative terms. Every non-default
# run also adds a RANDOM-GENE-SPLIT axis control (the pooled pole genes split at random into groups of the same sizes,
# fresh per cell), which is matched to the functional axis in how it is built.
ALPHA_UNIT = os.environ.get("ALPHA_UNIT", "raw")
ZERO_MASSIVE = os.environ.get("ZERO_MASSIVE") == "1"     # set the axis to 0 in the two largest (massive) dimensions
if MODEL != "217m" or XPREFIX != "ctx_maxtoki" or ALPHA_UNIT != "raw" or ZERO_MASSIVE:
    OUTNAME = (OUTNAME[:-5] + f"__{MODEL}_{XPREFIX}" + ("" if ALPHA_UNIT == "raw" else f"_{ALPHA_UNIT}")
               + ("_nomassive" if ZERO_MASSIVE else "") + ".json")
EXTRA = OUTNAME != ("ctx_causal.json" if (AXIS == "nuc_surf" and SITE == 3) else f"ctx_causal_{AXIS}_L{SITE+1:02d}.json")
# SPLITS=K (added 2 Oct 2026 after review): the matched control is K FIXED random splits of the pooled pole genes
# into groups of the pole sizes, each scored on the same cells, positions and push (strength 0.5) as the functional
# axis. The functional mean swing is compared with the K split mean swings (p = (1 + #{|split| >= |func|})/(K + 1)).
# This replaces the per-cell redrawn split (mean ~0 by construction, so not a valid null), which is then switched off.
SPLITS = int(os.environ.get("SPLITS", 0))
if SPLITS:
    EXTRA = False
    OUTNAME = OUTNAME[:-5] + "_splits.json"
N_CELLS, MAX_LEN, SEED = int(os.environ.get("N_CELLS", 30)), 512, 0


def tokenise_cells(tok, n):
    from maxtoki_adapter import MaxTokiTokenizer  # noqa
    with h5py.File(os.path.join(TS, PANEL), "r") as f:
        ens = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]).astype(str)
        ens = np.array([e.split(".")[0] for e in ens])
        var_idx, token_ids, medians = tok.make_var_mapping(list(ens))
        pos = np.full(len(ens), -1, np.int64); pos[var_idx] = np.arange(len(var_idx))
        import ctx_tokenise as TK
        X = TK.count_matrix(f); TK.check_counts(X); N = int(X.attrs["shape"][0]); indptr = X["indptr"][:]
        sel = np.sort(np.random.default_rng(SEED).choice(N, n, replace=False))
        seqs = []
        for r in sel:
            s, e = int(indptr[r]), int(indptr[r + 1]); idx, val = X["indices"][s:e], X["data"][s:e].astype(np.float32)
            keep = pos[idx] >= 0
            if not keep.any():
                continue
            j = pos[idx[keep]]
            order = TK.rank_order(val[keep], medians[j])
            if len(order) < 40:
                continue
            seqs.append(token_ids[j[order[: MAX_LEN - 2]]].astype(np.int64))
    return seqs


def main():
    ens2sym = {e: s.upper() for s, e in pickle.load(open(NAME_ID, "rb")).items()}
    g2g = {k.upper(): set(v) for k, v in pickle.load(open(G2G, "rb")).items() if isinstance(v, (set, list, tuple))}
    tokmap = json.load(open(f"{MSETUP}/token_dictionary.json"))
    ens2tid = {k: int(v) for k, v in tokmap.items()}

    # gene-token sets for the two poles (target/control readout + axis construction)
    def pole_tokens(terms):
        out = []
        for ens, tid in ens2tid.items():
            s = ens2sym.get(ens)
            if s and s in g2g and g2g[s] & set(terms):
                out.append(int(tid))
        return np.array(sorted(set(out)))
    nuc_tok, surf_tok = pole_tokens(NUC), pole_tokens(SURF)
    both = set(nuc_tok) & set(surf_tok)
    nuc_tok = np.array([t for t in nuc_tok if t not in both]); surf_tok = np.array([t for t in surf_tok if t not in both])
    print(f"[setup] {len(nuc_tok)} nuclear-pole tokens, {len(surf_tok)} surface-pole tokens", flush=True)

    # functional-context axis in RAW layer-4 hidden space (from the extraction; M is raw mean hidden state)
    z = np.load(os.path.join(RES, f"{XPREFIX}_L{SITE+1:02d}.npz"), allow_pickle=True)
    M, counts, cap, genes = z["M"].astype(np.float32), z["counts"], int(z["cap"]), z["genes"].astype(str)
    full = (counts == cap).all(0)
    araw = np.full((len(genes), M.shape[-1]), np.nan, np.float32)
    for gi in range(len(genes)):
        cs = np.where(full[:, gi])[0]
        if len(cs):
            araw[gi] = M[:, cs, gi].mean((0, 1))
    gsym = [ens2sym.get(g) for g in genes]
    ia = [i for i, s in enumerate(gsym) if np.isfinite(araw[i, 0]) and s in g2g and g2g[s] & set(NUC)]
    ib = [i for i, s in enumerate(gsym) if np.isfinite(araw[i, 0]) and s in g2g and g2g[s] & set(SURF)]
    both_g = set(ia) & set(ib); ia = [i for i in ia if i not in both_g]; ib = [i for i in ib if i not in both_g]
    u = araw[ia].mean(0) - araw[ib].mean(0)
    # how far does CONTEXT naturally move genes along this direction? (to compare with the push sizes below)
    uhat = u / (np.linalg.norm(u) + 1e-12)
    P = np.tensordot(M, uhat, axes=([3], [0])).mean(0)            # (n_ctx, n_gene) raw projection, partition-mean
    nat_sd, nat_rng = [], []
    for gi in range(len(genes)):
        cs = np.where(full[:, gi])[0]
        if len(cs) >= 2:
            nat_sd.append(float(P[cs, gi].std())); nat_rng.append(float(np.ptp(P[cs, gi])))
    natural = dict(n_genes=len(nat_sd), median_sd_across_contexts=float(np.median(nat_sd)),
                   median_range_across_contexts=float(np.median(nat_rng)),
                   p95_range_across_contexts=float(np.percentile(nat_rng, 95)),
                   pole_centroid_distance=float(np.linalg.norm(u)))
    del P
    # how much of the axis lies in the few largest (massive) hidden dimensions of the gene representations
    gm = np.nanmean(araw, 0); top2 = np.argsort(-np.abs(gm))[:2]
    natural["massive_dims"] = [int(x) for x in top2]
    natural["axis_share_in_massive_dims"] = float((u[top2] ** 2).sum() / (u ** 2).sum())
    if ZERO_MASSIVE:
        u = u.copy(); u[top2] = 0.0; natural["axis_massive_dims_zeroed"] = True
        natural["note"] = "natural-movement numbers describe the axis before the two dimensions were zeroed"
    pooled = np.array(ia + ib); n_a = len(ia)
    del M, z, counts, full
    import gc; gc.collect()
    st = SL.Steerer(model_dir=SL.MODELS[MODEL])               # model loaded only after the big arrays are freed
    d_func = SL.Direction(vec=u, name=f"{AXIS}@L{SITE+1}", basis=f"ctx_L{SITE+1}")
    d_rand = SL.random_direction(st.xt, seed=1, name="random")
    print(f"[setup] axis from {len(ia)} nuclear / {len(ib)} surface genes; ||u_raw||={np.linalg.norm(u):.2f}", flush=True)

    seqs = tokenise_cells(st.tok, N_CELLS)
    print(f"[setup] {len(seqs)} cells", flush=True)

    # calibrate alpha to the residual norm at the injection site
    h = st.hidden(np.concatenate([[st.tok.BOS], seqs[0], [st.tok.EOS]]), layer=SITE + 1)
    resnorm_raw = float(np.linalg.norm(h[1:-1], axis=1).mean())
    cn = []
    for sq in seqs[:5]:
        hh = np.asarray(st.hidden(np.concatenate([[st.tok.BOS], sq, [st.tok.EOS]]), layer=SITE + 1))[1:-1]
        cn.append(float(np.linalg.norm(hh - hh.mean(0), axis=1).mean()))
    resnorm_centred = float(np.mean(cn))
    resnorm = resnorm_raw if ALPHA_UNIT == "raw" else resnorm_centred
    print(f"[setup] residual norm raw {resnorm_raw:.1f}, after removing the across-gene mean {resnorm_centred:.1f} "
          f"(push unit: {ALPHA_UNIT})", flush=True)
    alphas = [0.0, 0.25, 0.5, 1.0, 2.0]
    alpha_units = [a * resnorm for a in alphas]              # in raw hidden-norm units
    print(f"[setup] residual norm at layer {SITE} ~ {resnorm:.1f}; alphas x that = {alpha_units}", flush=True)

    # SIGNED steering: +u (toward nuclear) vs -u (toward surface). A real causal channel FLIPS the readout with
    # the sign; the baseline nuclear/surface asymmetry and generic-perturbation effects do NOT depend on sign,
    # so they cancel in the swing = spec(+a) - spec(-a). Alphas kept in the non-destructive range.
    alphas_s = [0.25, 0.5, 1.0]
    au = [a * resnorm for a in alphas_s]
    rng = np.random.default_rng(SEED)
    swing = {"functional": {a: [] for a in alphas_s}, "random": {a: [] for a in alphas_s},
             "randsplit": {a: [] for a in alphas_s}}
    spec = {"functional_+": {a: [] for a in alphas_s}, "functional_-": {a: [] for a in alphas_s}}
    split_dirs = []
    for k in range(SPLITS):                                       # fixed splits, the same for every cell
        pr = np.random.default_rng(7000 + k).permutation(pooled)
        u_k = araw[pr[:n_a]].mean(0) - araw[pr[n_a:]].mean(0)
        if ZERO_MASSIVE:
            u_k = u_k.copy(); u_k[top2] = 0.0
        split_dirs.append(SL.Direction(vec=u_k, name=f"split{k}", basis=f"ctx_L{SITE+1}"))
    split_swing = np.zeros((SPLITS, len(seqs)))
    a_mid = au[alphas_s.index(0.5)]

    def spec_at(ids, d, a, steer_pos, read_pos):
        with st.steering(d, alpha=a, positions=steer_pos, site=SITE):
            lg = st.logits(ids)
        return float(np.mean([SL.mean_logit(lg[p], nuc_tok) for p in read_pos])
                     - np.mean([SL.mean_logit(lg[p], surf_tok) for p in read_pos]))

    for si, s in enumerate(seqs):
        ids = np.concatenate([[st.tok.BOS], s, [st.tok.EOS]]).astype(np.int64)
        gene_pos = np.arange(1, 1 + len(s)); rng.shuffle(gene_pos)
        steer_pos = list(gene_pos[: len(gene_pos) // 2]); read_pos = list(gene_pos[len(gene_pos) // 2:])
        d_rand_cell = SL.random_direction(st.xt, seed=1000 + si, name="random")   # FRESH per cell -> swing->0
        dirs = [("functional", d_func), ("random", d_rand_cell)]
        if EXTRA:                                                 # random split of the pooled pole genes
            pr = np.random.default_rng(5000 + si).permutation(pooled)
            u_rs = araw[pr[:n_a]].mean(0) - araw[pr[n_a:]].mean(0)
            dirs.append(("randsplit", SL.Direction(vec=u_rs, name="randsplit", basis=f"ctx_L{SITE+1}")))
        for label, d in dirs:
            pos_rows = SL.dose_response(st, ids, d, steer_pos, read_pos, nuc_tok, surf_tok, au, site=SITE)
            neg_rows = SL.dose_response(st, ids, d, steer_pos, read_pos, nuc_tok, surf_tok, [-a for a in au], site=SITE)
            for k, a in enumerate(alphas_s):
                swing[label][a].append(pos_rows[k]["specificity"] - neg_rows[k]["specificity"])
                if label == "functional":
                    spec["functional_+"][a].append(pos_rows[k]["specificity"])
                    spec["functional_-"][a].append(neg_rows[k]["specificity"])
        for k, d in enumerate(split_dirs):                        # baseline cancels in the signed swing
            split_swing[k, si] = (spec_at(ids, d, a_mid, steer_pos, read_pos)
                                  - spec_at(ids, d, -a_mid, steer_pos, read_pos))
        if si % 10 == 0:
            print(f"    cell {si}/{len(seqs)}", flush=True)

    from math import comb
    out = {"model": MODEL, "extraction": XPREFIX, "alpha_unit": ALPHA_UNIT, "resid_norm_raw": resnorm_raw,
           "resid_norm_centred": resnorm_centred,
           "alphas_xResidNorm": alphas_s, "site": SITE, "n_cells": len(seqs), "resid_norm": resnorm,
           "push_size_raw": [float(a) for a in au], "natural_movement_along_axis": natural,
           "n_nuclear_tokens": len(nuc_tok), "n_surface_tokens": len(surf_tok), "signed": {}}
    print(f"\n{'alpha':<8} {'spec(+u)':<12} {'spec(-u)':<12} {'FUNC swing':<20} {'RAND swing':<18} {'func>rand'}")
    for a in alphas_s:
        fp = np.array(spec["functional_+"][a]); fm = np.array(spec["functional_-"][a])
        fs = np.array(swing["functional"][a]); rs = np.array(swing["random"][a])
        k = int((fs > rs).sum()); n = len(fs)
        p = float(sum(comb(n, i) for i in range(k, n + 1)) / 2 ** n)
        out["signed"][f"alpha_{a}"] = dict(spec_plus=float(fp.mean()), spec_minus=float(fm.mean()),
                                           func_swing=float(fs.mean()), func_swing_sem=float(fs.std()/np.sqrt(n)),
                                           rand_swing=float(rs.mean()), func_gt_rand=k, n=n, sign_p=p)
        if EXTRA:
            ss_ = np.array(swing["randsplit"][a]); k2 = int((fs > ss_).sum())
            out["signed"][f"alpha_{a}"].update(
                randsplit_swing=float(ss_.mean()), randsplit_swing_sd=float(ss_.std()), func_gt_randsplit=k2,
                sign_p_vs_randsplit=float(sum(comb(n, i) for i in range(k2, n + 1)) / 2 ** n),
                z_vs_randsplit=float((fs.mean() - ss_.mean()) / (np.sqrt(fs.var() / n + ss_.var() / n) + 1e-12)))
        print(f"  {a:<6} {fp.mean():+.3f}       {fm.mean():+.3f}       {fs.mean():+.4f}±{fs.std()/np.sqrt(n):.4f}     "
              f"{rs.mean():+.4f}          {k}/{n} p={p:.1e}")

    mid = out["signed"]["alpha_0.5"]
    clean = (mid["spec_plus"] > 0 > mid["spec_minus"] and mid["func_swing"] > 2 * mid["func_swing_sem"]
             and mid["sign_p"] < 0.05 and mid["func_swing"] > 2 * abs(mid["rand_swing"]))
    out["verdict"] = (
        f"at 0.5xResid: spec(+u)={mid['spec_plus']:+.3f}, spec(-u)={mid['spec_minus']:+.3f}, "
        f"signed swing {mid['func_swing']:+.4f} (random {mid['rand_swing']:+.4f}), func>rand {mid['func_gt_rand']}/"
        f"{mid['n']} p={mid['sign_p']:.1e}. " +
        ("CAUSALLY USED — the readout FLIPS with the sign of the functional push (nuclear-gene logits rise toward "
         "+u, fall toward -u) at disjoint positions, dose-dependently, far beyond the sign-independent random "
         "control. Level 1 upgrades from 'represented' to 'represented AND causally used'."
         if clean else
         "NOT cleanly causal — the readout does not flip with the sign of the functional direction beyond the "
         "random control, so the earlier apparent effect was the baseline nuclear/surface asymmetry under "
         "generic perturbation. Level 1 stands as a representation result only.")
    )
    if SPLITS:
        fm = float(np.mean(swing["functional"][0.5])); sm = split_swing.mean(1)
        n_ge = int((np.abs(sm) >= abs(fm)).sum())
        out["split_null"] = dict(
            n_splits=SPLITS, alpha=0.5, func_mean_swing=fm, split_mean_swings=[float(x) for x in sm],
            split_abs_mean=float(np.abs(sm).mean()), split_abs_sd=float(np.abs(sm).std()),
            n_splits_abs_ge_func=n_ge, p_abs=float((1 + n_ge) / (SPLITS + 1)),
            z_abs=float((abs(fm) - np.abs(sm).mean()) / (np.abs(sm).std() + 1e-12)),
            n_splits_ge_func_signed=int((sm >= fm).sum()),
            p_signed=float((1 + int((sm >= fm).sum())) / (SPLITS + 1)))
        sn = out["split_null"]
        print(f"[split null] functional mean swing {fm:+.4f}; |split| mean {sn['split_abs_mean']:.4f} "
              f"(sd {sn['split_abs_sd']:.4f}); {n_ge}/{SPLITS} splits at least as large; p = {sn['p_abs']:.3f}; "
              f"z = {sn['z_abs']:+.1f}", flush=True)
    print(f"\nVERDICT: {out['verdict']}")
    json.dump(out, open(os.path.join(RES, OUTNAME), "w"), indent=1)
    print(f"[done] -> results/{OUTNAME}")


if __name__ == "__main__":
    main()

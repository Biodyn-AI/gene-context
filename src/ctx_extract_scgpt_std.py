"""scGPT (whole-human) per-gene contextual representations with scGPT's OWN pretraining input, on EXACTLY the cells,
cell types and partitions of the 600-cell MaxToki runs (results/ctx_cell_selection.npz, cells_600).

Replaces the legacy ctx_scgpt_L*.npz, which was reformatted from an old atlas: unbinned log values from the
ambient-corrected `X`, no <cls>, a different 3000-cell sample with 8 contexts, cap 20, and layer names one block off.

INPUT (checkpoint args.json: input_style=binned, n_bins=51, max_seq_len=1200, trunc_by_sample=true; scGPT
data_collator.py): per cell, the raw integer counts (`raw/X`) of all non-zero genes in scGPT's vocabulary (matched by
`var/feature_name`; duplicate symbols summed) are quantile-binned to 1..50 over ALL those genes (scGPT preprocess.binning,
random tie-breaking, seeded); then, if more than 1199 genes, 1199 are sampled at random (seeded); <cls> is prepended
with value pad_value -2 (scgpt/tasks/cell_emb.py, keep_first_n_tokens=1). Binning is invariant to per-cell scaling and log, so raw counts are the right input.
FORWARD: scGPT's own TransformerModel._encode (gene embedding + continuous value embedding, 12 post-norm blocks),
known genes only with full attention among them -- exactly how known genes were processed in generative pretraining.
Batch padding uses <pad> / pad_value -2 with the key-padding mask. dropout off (eval).

TAPS: tap k = residual after k blocks (output of layers[k-1]); tap 0 = input to block 0. Depth matching to
MaxToki-217M (11 blocks): 217M L4 (0.36) ~ scGPT L4 (0.33); 217M L2 (0.18) ~ scGPT L2 (0.17); 217M L8 (0.73) ~ scGPT L9.
GENES: the ctx217m600 panel (Ensembl), so all 600-cell representations share one gene list. Occurrence cap 50 per
(partition, context, gene), filled in a seeded random cell order (the cap is a random subsample).
RANDOM_INIT=1: same architecture, weights left at PyTorch's seeded default init (untrained control).

Out: results/{OUTPREFIX}_L{tap}.npz (M f16 [2, 12, n_genes, 512], counts, genes, contexts, cap, max_len, tap_convention)
     and results/{OUTPREFIX}_cells.npz (per-cell mean of the last tap, context, partition -> cell-type sanity check)
Run: ../../.venv_state/bin/python -u ctx_extract_scgpt_std.py
"""
import os, sys, json, time, types, warnings; warnings.filterwarnings("ignore")
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np, h5py, torch

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ctx_tokenise as TK
MI = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp"
SCGPT_PKG = os.path.join(MI, "external", "scGPT", "scgpt")
CKPT_DIR = os.path.join(MI, "external", "scGPT_checkpoints", "whole-human")
TS = os.path.join(MI, "data", "raw")
TAPS = [int(x) for x in os.environ.get("TAPS", "0,2,4,6,9,12").split(",")]
RANDOM_INIT = os.environ.get("RANDOM_INIT") == "1"
OUTPREFIX = os.environ.get("OUTPREFIX", "ctx_scgpt600rand" if RANDOM_INIT else "ctx_scgpt600")
REF = os.environ.get("REF_PREFIX", "ctx217m600")
CAP, N_BINS, MAX_SEQ, PAD_VALUE, BATCH, SEED = 50, 51, 1200, -2, int(os.environ.get("BATCH", 8)), 0
D_MODEL, N_LAYERS, N_HEADS = 512, 12, 8


def load_model_class():
    """import scgpt/model/model.py without running scgpt/__init__.py (which needs torchtext)"""
    for name, path in [("scgpt", SCGPT_PKG), ("scgpt.model", os.path.join(SCGPT_PKG, "model"))]:
        if name not in sys.modules:
            m = types.ModuleType(name); m.__path__ = [path]; sys.modules[name] = m
    from scgpt.model.model import TransformerModel
    return TransformerModel


def build_model(vocab):
    TransformerModel = load_model_class()
    torch.manual_seed(SEED)
    m = TransformerModel(ntoken=len(vocab), d_model=D_MODEL, nhead=N_HEADS, d_hid=D_MODEL, nlayers=N_LAYERS,
                         vocab=vocab, dropout=0.2, pad_token="<pad>", pad_value=PAD_VALUE,
                         input_emb_style="continuous", use_fast_transformer=False, do_mvc=True, do_dab=False,
                         use_batch_labels=False, cell_emb_style="cls", mvc_decoder_style="inner product", n_cls=1)
    if not RANDOM_INIT:
        ck = torch.load(os.path.join(CKPT_DIR, "best_model.pt"), map_location="cpu")
        sd = ck.get("model_state_dict", ck.get("model", ck)) if isinstance(ck, dict) else ck
        sd = {k.replace("Wqkv.", "in_proj_"): v for k, v in sd.items()}
        missing, unexpected = m.load_state_dict(sd, strict=False)
        print(f"[model] missing {missing}; unexpected {unexpected}", flush=True)
        enc_missing = [k for k in missing if k.startswith(("encoder.", "value_encoder.", "transformer_encoder."))]
        assert not enc_missing, f"encoder weights not loaded: {enc_missing[:5]}"
    m.transformer_encoder.enable_nested_tensor = False        # hooks must see plain (batch, seq, d) tensors
    m.transformer_encoder.use_nested_tensor = False
    torch.backends.mha.set_fastpath_enabled(False)            # plain math path (no nested-tensor fast path)
    return m.eval()


def scgpt_digitize(x, bins, rng):
    """scgpt/preprocess.py _digitize: random tie-breaking between left/right bin edges"""
    left = np.digitize(x, bins); right = np.digitize(x, bins, right=True)
    return np.ceil(rng.random(len(x)) * (right - left) + left).astype(np.int64)


def binning(vals, rng):
    bins = np.quantile(vals, np.linspace(0, 1, N_BINS - 1))
    return scgpt_digitize(vals, bins, rng)


def main():
    t0 = time.time()
    sel = np.load(os.path.join(HERE, "results", "ctx_cell_selection.npz"), allow_pickle=True)
    cells = sel["cells_600"]; panels = sel["panels"].astype(str); contexts = sel["contexts"].astype(str)
    ref = np.load(os.path.join(HERE, "results", f"{REF}_L04.npz"), allow_pickle=True)
    genes = np.array([g.split(".")[0] for g in ref["genes"].astype(str)])
    assert list(ref["contexts"].astype(str)) == list(contexts)
    vocab = json.load(open(os.path.join(CKPT_DIR, "vocab.json")))
    pad_id, cls_id = vocab["<pad>"], vocab["<cls>"]

    # ---- per panel: Ensembl -> symbol -> vocab id; panel gene index per vocab id ----
    gpos = {g: i for i, g in enumerate(genes)}
    seqs = [None] * len(cells)                                 # (ids, binned values) per selected cell
    rng_bin = np.random.default_rng(SEED)
    n_in_vocab, n_trunc, panel_cov = [], 0, None
    tid2panel = {}
    for pi, p in enumerate(panels):
        rows_here = np.where(cells[:, 0] == pi)[0]
        with h5py.File(os.path.join(TS, p), "r") as f:
            ens = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]).astype(str)
            ens = np.array([e.split(".")[0] for e in ens])
            fn = f["var"]["feature_name"]
            if isinstance(fn, h5py.Group):                     # categorical
                cats = np.array([x.decode() if isinstance(x, bytes) else x for x in fn["categories"][:]]).astype(str)
                sym = cats[fn["codes"][:]]
            else:
                sym = np.array([x.decode() if isinstance(x, bytes) else x for x in fn[:]]).astype(str)
            sym = np.array([s.split("_ENSG")[0] for s in sym])
            vid = np.array([vocab.get(s, vocab.get(s.upper(), -1)) for s in sym], np.int64)
            for e, v in zip(ens, vid):
                if v >= 0 and e in gpos:
                    tid2panel.setdefault(int(v), gpos[e])
            X = TK.count_matrix(f); TK.check_counts(X); indptr = X["indptr"][:]
            for k in rows_here:
                r = int(cells[k, 1]); s, e = int(indptr[r]), int(indptr[r + 1])
                idx, val = X["indices"][s:e], X["data"][s:e].astype(np.float32)
                v = vid[idx]; keep = (v >= 0) & (val > 0)
                u, inv = np.unique(v[keep], return_inverse=True)         # duplicate symbols -> one token, summed
                c = np.zeros(len(u), np.float32); np.add.at(c, inv, val[keep])
                b = binning(c, rng_bin).astype(np.float32)              # bin over ALL the cell's vocab genes first
                n_in_vocab.append(len(u))
                if len(u) > MAX_SEQ - 1:                                  # then sample (trunc_by_sample)
                    o = np.sort(rng_bin.permutation(len(u))[: MAX_SEQ - 1]); u, b = u[o], b[o]; n_trunc += 1
                seqs[k] = (np.concatenate([[cls_id], u]).astype(np.int64), np.concatenate([[PAD_VALUE], b]).astype(np.float32))
        print(f"[read] {p}: {len(rows_here)} cells ({time.time()-t0:.0f}s)", flush=True)
    covered = sum(1 for g in range(len(genes)) if g in set(tid2panel.values()))
    print(f"[input] median {int(np.median(n_in_vocab))} vocab genes per cell; {n_trunc} of {len(cells)} cells "
          f"sampled down to {MAX_SEQ-1}; {covered} of {len(genes)} panel genes in the scGPT vocabulary", flush=True)

    dev = torch.device("mps") if torch.backends.mps.is_available() else torch.device("cpu")
    model = build_model(vocab).to(dev)
    store = {}
    hooks = []
    for t in TAPS:
        if t == 0:
            hooks.append(model.transformer_encoder.layers[0].register_forward_pre_hook(
                lambda mod, args, t=t: store.__setitem__(t, args[0].detach())))
        else:
            hooks.append(model.transformer_encoder.layers[t - 1].register_forward_hook(
                lambda mod, inp, out, t=t: store.__setitem__(t, out.detach())))

    nP, nC, nG = 2, len(contexts), len(genes)
    acc = {t: np.zeros((nP, nC, nG, D_MODEL), np.float32) for t in TAPS}
    cnts = np.zeros((nP, nC, nG), np.int32)
    cell_emb = np.zeros((len(cells), D_MODEL), np.float32); cls_emb = np.zeros((len(cells), D_MODEL), np.float32)
    order = np.random.default_rng(SEED + 1).permutation(len(cells))
    for a in range(0, len(order), BATCH):
        chunk = order[a:a + BATCH]
        L = max(len(seqs[k][0]) for k in chunk)
        ids = np.full((len(chunk), L), pad_id, np.int64); vals = np.full((len(chunk), L), PAD_VALUE, np.float32)
        pm = np.ones((len(chunk), L), bool)
        for j, k in enumerate(chunk):
            n = len(seqs[k][0]); ids[j, :n] = seqs[k][0]; vals[j, :n] = seqs[k][1]; pm[j, :n] = False
        with torch.no_grad():
            store.clear()
            model._encode(torch.from_numpy(ids).to(dev), torch.from_numpy(vals).to(dev), torch.from_numpy(pm).to(dev))
            hs = {t: store[t].to("cpu", torch.float32).numpy() for t in TAPS}
        for j, k in enumerate(chunk):
            part, ci = int(cells[k, 3]), int(cells[k, 2]); n = len(seqs[k][0])
            cell_emb[k] = hs[TAPS[-1]][j, 1:n].mean(0); cls_emb[k] = hs[TAPS[-1]][j, 0]
            gi = np.array([tid2panel.get(int(t), -1) for t in seqs[k][0][1:]])   # position p+1 <-> gene p
            ok = np.where(gi >= 0)[0]
            ok = ok[cnts[part, ci, gi[ok]] < CAP]                  # each gene occurs once per cell -> exact cap
            cnts[part, ci, gi[ok]] += 1
            for t in TAPS:
                acc[t][part, ci, gi[ok]] += hs[t][j, 1 + ok]
        if (a // BATCH) % 100 == 0:
            print(f"    {a + len(chunk)}/{len(order)} cells | {float((cnts >= CAP).mean()):.1%} at cap "
                  f"({time.time()-t0:.0f}s)", flush=True)
    for h in hooks:
        h.remove()

    for t in TAPS:
        M = acc[t] / np.maximum(cnts[..., None], 1)
        out = os.path.join(HERE, "results", f"{OUTPREFIX}_L{t:02d}.npz")
        np.savez_compressed(out, M=M.astype(np.float16), counts=cnts, genes=genes, contexts=contexts, cap=CAP,
                            max_len=MAX_SEQ, cells_ctx=600, model="scGPT whole-human" + (" RANDOM-INIT" if RANDOM_INIT else ""),
                            tap_convention="tap k = residual after k of 12 blocks; 0 = block-0 input")
        print(f"  wrote {out}  M{M.shape}", flush=True)
    np.savez_compressed(os.path.join(HERE, "results", f"{OUTPREFIX}_cells.npz"), emb=cell_emb, cls=cls_emb, ctx=cells[:, 2],
                        part=cells[:, 3], tap=TAPS[-1])
    # sanity: nearest-centroid cell-type accuracy, centroids from partition 0, scored on partition 1 (chance 1/12)
    for nm, E_ in (("mean-pooled genes", cell_emb), ("<cls>", cls_emb)):
        Z = E_ / (np.linalg.norm(E_, axis=1, keepdims=True) + 1e-9)
        cen = np.stack([Z[(cells[:, 2] == c) & (cells[:, 3] == 0)].mean(0) for c in range(nC)])
        te = cells[:, 3] == 1
        acc_ct = float((np.argmax(Z[te] @ cen.T, 1) == cells[te, 2]).mean())
        print(f"[check] cell-type accuracy (nearest centroid) from tap {TAPS[-1]} {nm}: {acc_ct:.3f} (chance {1/nC:.3f})")
    print(f"[done] {float((cnts >= CAP).mean()):.1%} of (partition, context, gene) entries at cap {CAP}; "
          f"{(cnts == CAP).all(0).sum()} count-balanced (context, gene) pairs ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()

"""EXTRACT STATE SE-600M per-gene CONTEXTUAL representations with STATE's STANDARD input, on the SAME Tabula Sapiens
cells, contexts, partitions, processing order, occurrence cap and (by default) gene panel as ctx_extract_maxtoki.py.

WHY THIS REPLACES ctx_extract_state.py. The old script fed STATE the preprocessed HVG files (4,803-4,891 genes, only
1,864-2,414 of them in STATE's protein table, because half the HVGs are lncRNA/pseudogenes) with decontX counts, on a
greedy 600-cells-per-tissue subsample. So each cell sentence held only ~200-430 expressed genes (median), the rest
were random unexpressed HVGs, the gene universe differed by tissue, and cells/contexts/cap did not match MaxToki.

STATE's standard input (state/emb/data/loader.py, VCIDatasetSentenceCollator.sample_cell_sentences, l.457-629):
  * X = RAW integer counts over ANY gene set; genes are matched to the ESM2 protein table (19,790 HGNC symbols) by
    var_names, else by a symbol column (here raw/var 'feature_name': 19,549 of 60,606 genes match). Others are dropped.
  * per cell: log1p(counts) over the matched genes, ranked descending, ties broken at random (torch.randperm), top
    2047 genes -> positions 1..2047 (position 0 = CLS). If the cell has < 2047 expressed genes, the remaining slots
    hold zero-count genes. Count feature per position = 100 * log1p(c_g) / sum_valid log1p(c).
  * the model appends a learned dataset token at position 2048 (seq len 2049), has NO positional encoding and NO
    attention mask: the gene order is invisible to it; expression enters only through the soft-binned count feature.
We call the collator itself (in-process, one cell at a time, per-cell torch seed) so the input is exactly STATE's,
reproducible, and free of DataLoader worker processes (which crashed under memory pressure: route_state_geometry/RESULTS.md).

LAYER NAMING. L{k} = residual stream after k transformer blocks, k = 0..16 (same convention as HF hidden_states used
for MaxToki). k = 0 is the input to block 0 (ESM2 projection * sqrt(2048) + count embedding: it varies across contexts
ONLY through the count bin). k >= 1 is the output of transformer_encoder.layers[k-1] (post-norm). The old ctx_state_L04/
L08/L11 files were layers[4]/[8]/[11] outputs = after 5/9/12 blocks. Depth-matched to MaxToki-217M taps 0/4/8/11 (of 11
blocks): STATE 0/6/12/16.

VERIFY. Before any forward pass the script recomputes MaxToki's capped occurrence counts from its own replicated
selection and asserts they equal results/<MATCH>_Lxx.npz['counts'] exactly (and the same contexts and gene panel).
That proves the cells, partitions and processing order are identical. Mismatch -> abort.

Out: results/<OUTPREFIX>_L{k:02d}.npz  M[part, ctx, gene, 2048] f16, counts, genes (Ensembl), contexts, cap, max_len
     results/<OUTPREFIX>_cells.npz      per-cell manifest (file, row, ctx, part, assay, n_expressed) + CLS embeddings
Run: ../../.venv_state/bin/python -u ctx_extract_state_std.py            (DRYRUN=1 first: selection + verify only)
     MATCH=ctx217m600 CELLSCTX=600 MAXGENES=5000 ...                      (the 600-cell / 5000-gene comparison set)
     RANDOM_INIT=1 (random weights, real ESM2 table)  RANDOM_INIT=2 (random weights AND random gene table)
"""
import os, sys, glob, json, time, pickle, collections, warnings; warnings.filterwarnings("ignore")
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np, h5py, torch

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "route_state"))
MSETUP = "/Volumes/Crucial X6/MacBook/Code/biomi_automation/projects/maxtoki/setup"
sys.path.insert(0, MSETUP)

TS = "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/raw"
PANELS = ["tabula_sapiens_immune_subset_20000.h5ad", "tabula_sapiens_kidney.h5ad", "tabula_sapiens_lung.h5ad"]
GENE2GO = ("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/perturb/"
           "gene2go_all.pkl")

# ---- selection settings: MUST equal the MaxToki run named by MATCH (defaults = the A1 headline run) ------------
MATCH     = os.environ.get("MATCH", "ctx_maxtoki")    # ctx217m600 for the 600-cell / 5000-gene set
MAX_LEN   = int(os.environ.get("MAXLEN", 1024))       # MaxToki length: only used to replicate its selection + verify
N_CTX     = int(os.environ.get("NCTX", 12))
CELLS_CTX = int(os.environ.get("CELLSCTX", 1000))
FLOOR     = int(os.environ.get("FLOOR", 25))
MAX_GENES = int(os.environ.get("MAXGENES", 6000))
CAP       = int(os.environ.get("CAP", 50))
SEED, NPART = 0, 2
# ---- STATE settings --------------------------------------------------------------------------------------------
TAPS        = [int(x) for x in os.environ.get("TAPS", "0,6,12,16").split(",")]
BATCH       = int(os.environ.get("BATCH", 4))
DEVICE      = os.environ.get("STATE_DEVICE", "mps")
PANEL_MODE  = os.environ.get("PANEL", "maxtoki")      # maxtoki = MaxToki's panel (mapped); own = same rule on STATE
RANDOM_INIT = int(os.environ.get("RANDOM_INIT", 0))
OUTPREFIX   = os.environ.get("OUTPREFIX", {0: "ctx_statestd", 1: "ctx_staterand", 2: "ctx_staterandpe"}[RANDOM_INIT])
EXTRA_CAP   = int(os.environ.get("EXTRA_CAP", 0))
EXTRA_TAPS  = [int(x) for x in os.environ.get("EXTRA_TAPS", "6").split(",")] if EXTRA_CAP else []
EXTRA_PREFIX = os.environ.get("EXTRA_PREFIX", f"{OUTPREFIX}_cap{EXTRA_CAP}")
CKPT_EVERY  = int(os.environ.get("CKPT_EVERY", 2000))  # cells between resumable checkpoints (0 = off)
DRYRUN      = int(os.environ.get("DRYRUN", 0))
READOUT_N   = int(os.environ.get("READOUT_N", 0))      # >0: gate STATE's expression decoder on this many cells
NAME = "ts"                                            # dataset name inside STATE's collator (one shared mapping)
NBLOCKS = 16


def _dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a]).astype(str)


def _col(g, name):
    o = g[name]
    return _dec(o["categories"][:])[o["codes"][:]] if isinstance(o, h5py.Group) else _dec(o[:])


def uid_of(pi, row):
    return pi * 1_000_000 + int(row)                  # stable cell id -> per-cell seed (independent of order/batch)


# =========================== 1. replicate ctx_extract_maxtoki's cell selection ===================================
def stream_rows(path, tok):
    """IDENTICAL to ctx_extract_maxtoki.stream_panel (same skips, same truncation), but also yields the row index."""
    import ctx_tokenise as TK
    with h5py.File(path, "r") as f:
        ens = np.array([e.split(".")[0] for e in _dec(f["var"]["_index"][:])])
        ctg = f["obs"]["cell_type"]
        ctypes = _dec(ctg["categories"][:])[ctg["codes"][:]]
        X = TK.count_matrix(f); TK.check_counts(X); n = int(X.attrs["shape"][0]); indptr = X["indptr"][:]
        var_idx, token_ids, medians = tok.make_var_mapping(list(ens))
        pos = np.full(len(ens), -1, np.int64); pos[var_idx] = np.arange(len(var_idx))
        for r in range(n):
            s, e = int(indptr[r]), int(indptr[r + 1])
            idx, val = X["indices"][s:e], X["data"][s:e].astype(np.float32)
            keep = pos[idx] >= 0
            if not keep.any():
                continue
            j = pos[idx[keep]]
            order = TK.rank_order(val[keep], medians[j])[: MAX_LEN - 2]
            if not len(order):
                continue
            yield r, token_ids[j[order]].astype(np.int16), ctypes[r]     # MaxToki ids < 20,275 fit in int16


def select_cells(tok):
    """Same rng (default_rng(0)), same calls, same order as ctx_extract_maxtoki.main. Lists of tuples on purpose:
    Generator.shuffle on a list consumes the rng exactly as it did on MaxToki's list of token arrays."""
    rng = np.random.default_rng(SEED)
    cells = collections.defaultdict(list)
    for pi, p in enumerate(PANELS):
        path = os.path.join(TS, p)
        if not os.path.exists(path):
            print(f"  [skip] {p}"); continue
        k = 0
        for r, toks, ct in stream_rows(path, tok):
            cells[ct].append((pi, r, toks)); k += 1
        print(f"  {p}: {k} cells", flush=True)
    ctx_names = sorted([c for c in cells if len(cells[c]) >= 250], key=lambda c: -len(cells[c]))[:N_CTX]
    for c in ctx_names:
        rng.shuffle(cells[c])
        cells[c] = cells[c][:CELLS_CTX]
    cnt = {c: collections.Counter(int(t) for _, _, s in cells[c] for t in s) for c in ctx_names}
    reach = collections.Counter()
    for c in ctx_names:
        for g, n in cnt[c].items():
            if n >= FLOOR:
                reach[g] += 1
    panel = [g for g, k in reach.items() if k >= 2]
    panel.sort(key=lambda g: -sum(cnt[c].get(g, 0) for c in ctx_names))
    panel_tok = sorted(panel[:MAX_GENES])
    work = [(ci, si, pi, r, s) for ci, c in enumerate(ctx_names) for si, (pi, r, s) in enumerate(cells[c])]
    rng.shuffle(work)
    return ctx_names, panel_tok, work


def verify_against_maxtoki(ctx_names, panel_tok, work):
    """Recompute MaxToki's capped counts (no forward pass) and require exact equality with the MaxToki output."""
    cands = sorted(glob.glob(os.path.join(HERE, "results", f"{MATCH}_L??.npz")))
    assert cands, f"no results/{MATCH}_L??.npz to verify against"
    z = np.load(cands[0], allow_pickle=True)            # only the small members below are decompressed
    tid2ens = {int(v): k for k, v in json.load(open(f"{MSETUP}/token_dictionary.json")).items()}
    assert int(z["max_len"]) == MAX_LEN, (int(z["max_len"]), MAX_LEN)
    assert list(z["contexts"]) == list(ctx_names), "contexts differ"
    ens = np.array([tid2ens.get(g, str(g)) for g in panel_tok])
    assert np.array_equal(z["genes"].astype(str), ens), "gene panel differs (check MAXGENES / FLOOR / CELLSCTX)"
    cap = int(z["cap"]); gpos = {g: i for i, g in enumerate(panel_tok)}
    cnts = np.zeros((NPART, len(ctx_names), len(panel_tok)), np.int32)
    for ci, si, _, _, s in work:
        part = si % NPART
        for t in s:
            gi = gpos.get(int(t))
            if gi is None or cnts[part, ci, gi] >= cap:
                continue
            cnts[part, ci, gi] += 1
    assert np.array_equal(cnts, z["counts"]), "capped counts differ: cells / partitions / order are NOT identical"
    print(f"[verify] identical to {os.path.basename(cands[0])}: {len(ctx_names)} contexts, {len(panel_tok)} genes, "
          f"{len(work)} cells, cap {cap}", flush=True)
    return ens


# =========================== 2. STATE's standard input ===========================================================
def state_gene_mapping(keys):
    """STATE's own rule (loader.py FilteredGenesCounts l.216-266): var_names are Ensembl here (0 hits), so it falls
    back to the symbol column; feature_name -> row of the protein table, -1 if absent. Same for all three files."""
    gpos = {g: i for i, g in enumerate(keys)}
    ref = None
    for p in PANELS:
        with h5py.File(os.path.join(TS, p), "r") as f:
            rv = f["raw"]["var"] if "raw" in f else f["var"]
            fn, ens = _col(rv, "feature_name"), _col(rv, "_index")
        if ref is None:
            ref = (fn, ens)
        else:
            assert np.array_equal(fn, ref[0]) and np.array_equal(ens, ref[1]), f"gene order differs in {p}"
    fn, ens = ref
    mapping = np.array([gpos.get(g, -1) for g in fn], np.int64)
    valid = mapping >= 0
    c = collections.Counter(fn[valid])
    dup_g = np.array(sorted({gpos[s] for s, n in c.items() if n > 1}), np.int64)   # 16 symbols appear twice
    return fn, ens, mapping, valid, dup_g


def make_collator(cfg, mapping, valid):
    from state.emb.data.loader import VCIDatasetSentenceCollator
    col = VCIDatasetSentenceCollator(cfg, valid_gene_mask={NAME: valid}, ds_emb_mapping_inference={NAME: mapping},
                                     is_train=False)
    col.training = False
    return col


def collate(col, items):
    """items: list of (dense raw-count row over ALL var genes, uid). Calls STATE's collator once per cell under a
    per-cell seed (same sentence in every pass, independent of batch composition), then stacks the 9-tuple."""
    outs = []
    for v, uid in items:
        torch.manual_seed(SEED * 1_000_003 + uid)
        outs.append(col([(torch.from_numpy(v).reshape(1, -1), uid, NAME, 0)]))
    return tuple(None if outs[0][k] is None else torch.cat([o[k] for o in outs], 0) for k in range(len(outs[0])))


class RowReader:
    """Dense raw/X rows (all 60,606 genes, var order untouched -- the CLS count quirk depends on it)."""
    def __init__(self):
        import ctx_tokenise as TK
        self.f = [h5py.File(os.path.join(TS, p), "r") for p in PANELS]
        self.X = [TK.count_matrix(f) for f in self.f]
        for X in self.X:
            TK.check_counts(X)
        self.ip = [X["indptr"][:] for X in self.X]
        self.G = int(self.X[0].attrs["shape"][1])
        self.assay = [_col(f["obs"], "assay") for f in self.f]

    def row(self, pi, r):
        s, e = int(self.ip[pi][r]), int(self.ip[pi][r + 1])
        v = np.zeros(self.G, np.float32)
        v[self.X[pi]["indices"][s:e]] = self.X[pi]["data"][s:e]
        return v


def logbranch(v):
    """True if STATE's is_raw_integer_counts heuristic (loader.py l.431-453) would wrongly treat this row as log1p."""
    return v.max() <= 35 and np.expm1(v.astype(np.float64)).sum() <= 5_000_000


# =========================== 3. model ============================================================================
def build_state(dev, random_init=0):
    """0: trained SE-600M (state_loader.load_state_se). 1: same constructor, PyTorch default init (trainer.py adds no
    custom init), real frozen ESM2 table. 2: as 1 but the ESM2 table replaced by a seeded Gaussian of the same shape
    (mirrors ctx_extract_random.py, where MaxToki's token embeddings are random). dropout MUST be 0: the attention
    calls F.scaled_dot_product_attention(dropout_p=self.dropout) which applies dropout even in eval mode."""
    import torch.nn as nn
    from state_loader import load_state_se, load_protein_embeds, build_cfg
    if random_init == 0:
        model, cfg, genes, info = load_state_se(device=dev, dtype=torch.float32)
        assert not info["missing"], info["missing"][:5]
        return model, cfg, genes
    from state.emb.nn.model import StateEmbeddingModel
    cfg = build_cfg(); emb = cfg.embeddings[cfg.embeddings.current]
    _, genes, pe_mat = load_protein_embeds()
    torch.manual_seed(SEED)
    model = StateEmbeddingModel(token_dim=int(emb.size), d_model=int(cfg.model.emsize), nhead=int(cfg.model.nhead),
                                d_hid=int(cfg.model.d_hid), nlayers=int(cfg.model.nlayers),
                                output_dim=int(cfg.model.output_dim), dropout=0.0, emb_size=int(emb.size),
                                collater=None, cfg=cfg)
    if random_init == 2:
        pe_mat = torch.randn(pe_mat.shape, generator=torch.Generator().manual_seed(SEED + 1))
    model.pe_embedding = nn.Embedding.from_pretrained(pe_mat)
    model.eval().to(device=dev, dtype=torch.float32)
    return model, cfg, genes


def add_taps(model, taps, store):
    """L{k} = after k blocks: k=0 -> input of layers[0] (pre-hook); k>=1 -> output of layers[k-1]."""
    hs = []
    for k in taps:
        assert 0 <= k <= NBLOCKS, k
        if k == 0:
            hs.append(model.transformer_encoder.layers[0].register_forward_pre_hook(
                lambda m, a: store.__setitem__(0, a[0].detach())))
        else:
            hs.append(model.transformer_encoder.layers[k - 1].register_forward_hook(
                (lambda kk: (lambda m, a, o: store.__setitem__(kk, o.detach())))(k)))
    return hs


def decode_genes(model, emb, dataset_emb, Y, gene_gidx):
    """STATE's own per-gene expression readout, assembled exactly as in StateEmbeddingModel.shared_step (model.py
    l.372-402): binary_decoder(cat[gene query, cell embedding, mu, ds10]) -> predicted log1p(count), shape (B, G).
    gene query = gene_embedding_layer(ESM2[g]) WITHOUT L2-normalising (as in training); mu = mean of the non-zero
    task counts (batch[2]); ds10 = dataset_embedder(decoder output at the dataset token).
    Do NOT use Inference.decode_from_adata: it concatenates the 2058-d transform output + ds again (4117 != 4107)
    and fixes read depth at 4.0."""
    dev = emb.device
    Xg = model.gene_embedding_layer(model.pe_embedding(gene_gidx.to(dev)))                       # (G, 2048)
    mu = torch.nan_to_num(torch.nanmean(Y.float().masked_fill(Y == 0, float("nan")), dim=1), nan=0.0).to(dev)
    ds10 = model.dataset_embedder(dataset_emb)                                                   # (B, 10)
    B, G = emb.shape[0], Xg.shape[0]
    comb = torch.cat([Xg.unsqueeze(0).expand(B, G, -1), emb.unsqueeze(1).expand(B, G, -1),
                      mu.view(B, 1, 1).expand(B, G, 1), ds10.unsqueeze(1).expand(B, G, -1)], dim=2)
    return model.binary_decoder(comb).squeeze(-1)


# =========================== 4. main =============================================================================
def main():
    from maxtoki_adapter import MaxTokiTokenizer
    from state_loader import build_cfg
    t0 = time.time()
    print(f"[pass 1] replicating {MATCH}'s cell selection (MaxToki tokenisation, CPU)", flush=True)
    ctx_names, panel_tok, work = select_cells(MaxTokiTokenizer(model_input_size=MAX_LEN))
    panel_ens_mt = verify_against_maxtoki(ctx_names, panel_tok, work)
    for ci, c in enumerate(ctx_names):
        print(f"    {sum(1 for w in work if w[0] == ci):>5} cells  {c[:58]}")

    # ---- STATE gene mapping + panel -----------------------------------------------------------------------------
    from state_loader import MODEL_DIR
    import zipfile, io
    class _Keys(pickle.Unpickler):                    # protein-table KEYS only, without loading the 400 MB tensors
        def find_class(self, mod, name):
            if mod.startswith("torch"):
                return lambda *a, **k: None
            return super().find_class(mod, name)
        def persistent_load(self, pid):
            return None
    keys = list(_Keys(io.BytesIO(zipfile.ZipFile(MODEL_DIR / "protein_embeddings.pt").read(
        "archive/data.pkl"))).load().keys())
    fn, ens, mapping, valid, dup_g = state_gene_mapping(keys)
    print(f"[state] {valid.sum()} of {len(fn)} TS genes map to the {len(keys)}-gene ESM2 table; "
          f"{len(dup_g)} symbols occur twice (excluded from the panel)", flush=True)
    ens2g = {e: int(m) for e, m in zip(ens, mapping) if m >= 0}
    g2ens, g2col = {}, {}
    for c, (e, m) in enumerate(zip(ens, mapping)):
        if m >= 0:
            g2ens.setdefault(int(m), e); g2col.setdefault(int(m), c)
    dup_set = set(dup_g.tolist())
    cfg = build_cfg()
    col = make_collator(cfg, mapping, valid)
    rr = RowReader()

    if PANEL_MODE == "maxtoki":
        panel = [ens2g.get(e, -1) for e in panel_ens_mt]
        panel = np.array([g for g in panel if g >= 0 and g not in dup_set], np.int64)
        print(f"[panel] MaxToki panel: {len(panel)} of {len(panel_ens_mt)} genes usable in STATE "
              f"(dropped: not in table / duplicated symbol)", flush=True)
    else:                                              # same FLOOR / MAX_GENES rule on STATE's own sentences
        occ = np.zeros((len(ctx_names), len(keys)), np.int64)
        for a in range(0, len(work), 64):
            chunk = work[a:a + 64]
            b = collate(col, [(rr.row(pi, r), uid_of(pi, r)) for _, _, pi, r, _ in chunk])
            bs, cw = b[0].numpy(), b[7].numpy()
            for j, (ci, *_rest) in enumerate(chunk):
                g = np.unique(bs[j, 1:][cw[j, 1:] > 0])
                occ[ci, g] += 1
        reach = (occ >= FLOOR).sum(0)
        cand = [g for g in np.where(reach >= 2)[0] if g not in dup_set]
        cand.sort(key=lambda g: -occ[:, g].sum())
        panel = np.array(sorted(cand[:MAX_GENES]), np.int64)
        print(f"[panel] own STATE panel: {len(panel)} genes reaching {FLOOR} in >=2 contexts", flush=True)
    g2panel = np.full(len(keys), -1, np.int64); g2panel[panel] = np.arange(len(panel))
    panel_syms = np.array([keys[g] for g in panel]); panel_ens = np.array([g2ens[int(g)] for g in panel])

    # input statistics (documentation only)
    nlog = 0; nexp = []
    for _, _, pi, r, _ in work[:500]:
        v = rr.row(pi, r); nlog += logbranch(v); nexp.append(int(((v > 0) & valid).sum()))
    nexp = np.array(nexp)
    print(f"[input] first 500 cells: expressed table genes median {np.median(nexp):.0f}; "
          f"{(nexp > 2047).mean():.0%} truncated to 2047; {nlog} would hit STATE's log1p branch", flush=True)
    print(f"[pass 1] done in {(time.time() - t0) / 60:.1f} min", flush=True)
    if DRYRUN:
        return

    # ---- PASS 2: forward + capped accumulation ------------------------------------------------------------------
    dev = DEVICE if (DEVICE != "mps" or torch.backends.mps.is_available()) else "cpu"
    model, cfg_m, genes_m = build_state(dev, RANDOM_INIT)
    assert list(genes_m) == keys
    if not READOUT_N:
        model.binary_decoder = None; model.dataset_encoder = None          # ~0.7 GB not needed for extraction
    store = {}
    hooks = add_taps(model, TAPS, store)
    d = int(cfg_m.model.emsize)
    C, G = len(ctx_names), len(panel)
    ck_dir = os.path.join(HERE, "results", f"_ckpt_{OUTPREFIX}")
    acc = {k: np.zeros((NPART, C, G, d), np.float32) for k in TAPS}
    cnts = np.zeros((NPART, C, G), np.int32)
    wsum = np.zeros((NPART, C, G), np.float64)     # sum of STATE's count feature (100*log1p share): expression covariate
    acc2 = {k: np.zeros((NPART, C, G, d), np.float32) for k in EXTRA_TAPS}
    cnts2 = np.zeros((NPART, C, G), np.int32)
    cls = np.zeros((len(work), d), np.float16); nexp_all = np.zeros(len(work), np.int32)
    start = 0
    settings = dict(MATCH=MATCH, CELLS_CTX=CELLS_CTX, MAX_GENES=MAX_GENES, FLOOR=FLOOR, CAP=CAP, TAPS=TAPS,
                    RANDOM_INIT=RANDOM_INIT, PANEL_MODE=PANEL_MODE, EXTRA_CAP=EXTRA_CAP, EXTRA_TAPS=EXTRA_TAPS,
                    n_work=len(work), n_panel=int(len(panel)))
    if CKPT_EVERY and os.path.exists(os.path.join(ck_dir, "state.json")):
        st = json.load(open(os.path.join(ck_dir, "state.json")))
        if st.get("settings") != settings:                 # never mix a checkpoint from other settings or code
            raise SystemExit(f"checkpoint {ck_dir} was written with settings {st.get('settings')}, not {settings}; "
                             "move it away to start fresh")
        start = int(st["done"])
        for k in TAPS:
            acc[k] = np.load(os.path.join(ck_dir, f"acc_L{k:02d}.npy"))
        for k in EXTRA_TAPS:
            acc2[k] = np.load(os.path.join(ck_dir, f"acc2_L{k:02d}.npy"))
        cnts = np.load(os.path.join(ck_dir, "cnts.npy")); cnts2 = np.load(os.path.join(ck_dir, "cnts2.npy"))
        cls = np.load(os.path.join(ck_dir, "cls.npy")); nexp_all = np.load(os.path.join(ck_dir, "nexp.npy"))
        wsum = np.load(os.path.join(ck_dir, "wsum.npy"))
        print(f"[resume] from cell {start}", flush=True)
    print(f"[pass 2] {'trained' if not RANDOM_INIT else f'RANDOM_INIT={RANDOM_INIT}'} STATE on {dev}; taps {TAPS} "
          f"(after-k-blocks); accumulator {sum(a.nbytes for a in acc.values()) / 2**30:.2f} GB", flush=True)

    def checkpoint(done):
        os.makedirs(ck_dir, exist_ok=True)
        for k in TAPS:
            np.save(os.path.join(ck_dir, f"acc_L{k:02d}.npy"), acc[k])
        for k in EXTRA_TAPS:
            np.save(os.path.join(ck_dir, f"acc2_L{k:02d}.npy"), acc2[k])
        for nm, arr in (("cnts", cnts), ("cnts2", cnts2), ("cls", cls), ("nexp", nexp_all), ("wsum", wsum)):
            np.save(os.path.join(ck_dir, f"{nm}.npy"), arr)
        json.dump({"done": done, "settings": settings}, open(os.path.join(ck_dir, "state.json"), "w"))

    t1 = time.time(); last_ck = start
    for a in range(start, len(work), BATCH):
        chunk = work[a:a + BATCH]
        batch = collate(col, [(rr.row(pi, r), uid_of(pi, r)) for _, _, pi, r, _ in chunk])
        with torch.no_grad():
            _, _, _, emb, _ = model._compute_embedding_for_batch(batch)
            res = {k: store[k].float().cpu().numpy() for k in TAPS}          # (B, 2049, d); dataset token at 2048
        cls[a:a + len(chunk)] = emb.float().cpu().numpy().astype(np.float16)
        bs, cw = batch[0].numpy(), batch[7].numpy()
        for j, (ci, si, pi, r, _) in enumerate(chunk):
            part = si % NPART
            pos = np.nonzero(cw[j, 1:] > 0)[0] + 1                        # expressed genes only; skip CLS at 0
            nexp_all[a + j] = len(pos)
            gi = g2panel[bs[j, pos]]
            keep = gi >= 0; pos, gi = pos[keep], gi[keep]
            gi, first = np.unique(gi, return_index=True); pos = pos[first]  # one token per gene
            ok = cnts[part, ci, gi] < CAP
            pos1, gi1 = pos[ok], gi[ok]
            if len(gi1):
                for k in TAPS:
                    acc[k][part, ci, gi1] += res[k][j, pos1]                # gi1 unique -> fancy += is exact
                cnts[part, ci, gi1] += 1
                wsum[part, ci, gi1] += cw[j, pos1]
            if EXTRA_CAP:
                ok2 = ok & (cnts2[part, ci, gi] < EXTRA_CAP)               # first EXTRA_CAP of the same occurrences
                if ok2.any():
                    for k in EXTRA_TAPS:
                        acc2[k][part, ci, gi[ok2]] += res[k][j, pos[ok2]]
                    cnts2[part, ci, gi[ok2]] += 1
        done = a + len(chunk)
        if done % 400 < BATCH:
            el = time.time() - t1
            print(f"    {done}/{len(work)} cells | {el / max(done - start, 1):.2f} s/cell | eta "
                  f"{el / max(done - start, 1) * (len(work) - done) / 60:.0f} min | "
                  f"{float((cnts >= CAP).mean()):.1%} of (gene,ctx,part) at cap", flush=True)
        if CKPT_EVERY and done - last_ck >= CKPT_EVERY and done < len(work):
            checkpoint(done); last_ck = done
    for h in hooks:
        h.remove()

    # ---- save (ctx_maxtoki schema + provenance) -----------------------------------------------------------------
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    common = dict(counts=cnts, genes=panel_ens, mean_count_weight=(wsum / np.maximum(cnts, 1)).astype(np.float32), contexts=np.array(ctx_names), cap=CAP, max_len=2048, cells_ctx=CELLS_CTX,
                  symbols=panel_syms, layer_convention="after_k_blocks", n_blocks=NBLOCKS,
                  model=["SE-600M", "random_init_real_esm2", "random_init_random_table"][RANDOM_INIT],
                  input="raw/X integer counts, STATE VCIDatasetSentenceCollator (pad 2048, per-cell seed)",
                  matched_to=MATCH, panel_mode=PANEL_MODE)
    for k in TAPS:
        M = acc[k] / np.maximum(cnts[..., None], 1)
        out = os.path.join(HERE, "results", f"{OUTPREFIX}_L{k:02d}.npz")
        np.savez_compressed(out, M=M.astype(np.float16), tap=k, **common)
        print(f"  wrote {out}  M{M.shape}", flush=True)
    for k in EXTRA_TAPS:
        M2 = acc2[k] / np.maximum(cnts2[..., None], 1)
        out = os.path.join(HERE, "results", f"{EXTRA_PREFIX}_L{k:02d}.npz")
        np.savez_compressed(out, M=M2.astype(np.float16), tap=k, **{**common, "counts": cnts2, "cap": EXTRA_CAP})
        print(f"  wrote {out}  (cap {EXTRA_CAP})", flush=True)
    np.savez_compressed(os.path.join(HERE, "results", f"{OUTPREFIX}_cells.npz"),
                        h5file=np.array([PANELS[w[2]] for w in work]), row=np.array([w[3] for w in work]),
                        ctx=np.array([w[0] for w in work]), part=np.array([w[1] % NPART for w in work]),
                        assay=np.array([rr.assay[w[2]][w[3]] for w in work]), n_expressed=nexp_all, cls=cls,
                        contexts=np.array(ctx_names))
    if os.path.isdir(ck_dir):                               # finished: retire the checkpoint so it is never resumed
        os.replace(ck_dir, ck_dir + "_finished")
    print(f"[done] {G} genes x {C} contexts x {NPART} partitions x {len(TAPS)} taps; "
          f"{float((cnts >= CAP).mean()):.1%} at the {CAP} cap; {(time.time() - t1) / 60:.0f} min forward", flush=True)
    del acc, acc2; import gc; gc.collect()                                  # validate() reloads ~1.2 GB

    sanity_cls(cls, np.array([w[0] for w in work]))
    validate(OUTPREFIX, panel_syms, 6 if 6 in TAPS else TAPS[-1])
    if READOUT_N:
        readout_gate(model, col, rr, work[:READOUT_N], panel, g2col)


# =========================== 5. checks ===========================================================================
def sanity_cls(cls, y):
    """Cell-type decodability of STATE's own cell embedding (5-fold logistic regression). Chance = 1/n_contexts."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    acc = cross_val_score(LogisticRegression(max_iter=2000), cls.astype(np.float32), y, cv=5).mean()
    print(f"[sanity] cell type from CLS embedding: 5-fold accuracy {acc:.3f} (chance {1 / len(np.unique(y)):.3f})")


def _auc(labels, scores):
    from scipy.stats import rankdata
    labels = np.asarray(labels); r = rankdata(scores)                       # average ranks for ties
    n1 = labels.sum(); n0 = len(labels) - n1
    return float("nan") if n1 == 0 or n0 == 0 else float((r[labels == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def validate(prefix, panel_syms, tap):
    """Nuclear-vs-surface GO axis (as ctx_extract_state.validate): wrong gene->position mapping gives AUC ~0.5."""
    z = np.load(os.path.join(HERE, "results", f"{prefix}_L{tap:02d}.npz"), allow_pickle=True)
    M = z["M"].astype(np.float32); counts = z["counts"]
    tot = counts.sum(0)
    rep = (M[0] * counts[0][..., None] + M[1] * counts[1][..., None]) / np.maximum(tot[..., None], 1)
    valid = tot > 0
    ent = rep[valid]; repz = (rep - ent.mean(0)) / (ent.std(0) + 1e-8)
    G = rep.shape[1]; a = np.zeros((G, rep.shape[2]), np.float32); ok = np.zeros(G, bool)
    for g in range(G):
        cv = np.where(valid[:, g])[0]
        if len(cv):
            a[g] = repz[cv, g].mean(0); ok[g] = True
    g2g = pickle.load(open(GENE2GO, "rb"))
    NUC = {"GO:0005634", "GO:0000785", "GO:0003677"}; SUR = {"GO:0005886", "GO:0005576", "GO:0005615"}
    lab = np.full(G, -1)
    for g in range(G):
        gs = g2g.get(panel_syms[g]) or g2g.get(panel_syms[g].upper()) or set()
        n, u = bool(gs & NUC), bool(gs & SUR)
        if ok[g] and n != u:
            lab[g] = int(n)
    pole = np.where(lab >= 0)[0]; y = lab[pole]; X = a[pole]
    folds = np.array_split(np.random.default_rng(0).permutation(len(pole)), 5); ys, ss = [], []
    for k in range(5):
        te = folds[k]; tr = np.concatenate([folds[j] for j in range(5) if j != k])
        u = X[tr][y[tr] == 1].mean(0) - X[tr][y[tr] == 0].mean(0)
        ys.append(y[te]); ss.append(X[te] @ (u / (np.linalg.norm(u) + 1e-12)))
    auc = _auc(np.concatenate(ys), np.concatenate(ss))
    print(f"[validate] L{tap} nuclear-vs-surface AUC {auc:.3f} on {int(y.sum())}+{int((1 - y).sum())} pole genes "
          f"-> {'PASS' if auc > 0.60 else 'FAIL: check the gene->position mapping'}")


def readout_gate(model, col, rr, cells, panel, g2col):
    """Gate STATE's decoder before using it as a steering readout. The bar is CELL-SPECIFIC signal: for each panel
    gene, Pearson r ACROSS CELLS between predicted and observed log1p count (a decoder that only knows each gene's
    average level scores ~0 here), against a cell-shuffled control. Also reports the per-cell Spearman over genes
    for the model and for the per-gene-mean predictor (the easy part)."""
    from scipy.stats import spearmanr
    gidx = torch.from_numpy(panel); cols = np.array([g2col[int(g)] for g in panel])
    P_all, O_all = [], []
    for a in range(0, len(cells), BATCH):
        chunk = cells[a:a + BATCH]
        rows = [rr.row(pi, r) for _, _, pi, r, _ in chunk]
        b = collate(col, [(v, uid_of(pi, r)) for v, (_, _, pi, r, _) in zip(rows, chunk)])
        with torch.no_grad():
            _, Y, _, emb, ds = model._compute_embedding_for_batch(b)
            P = torch.cat([decode_genes(model, emb, ds, Y, gidx[s:s + 500]) for s in range(0, len(gidx), 500)], 1)
        P_all.append(P.float().cpu().numpy()); O_all.append(np.log1p(np.stack([v[cols] for v in rows])))
    P, O = np.concatenate(P_all), np.concatenate(O_all)

    def across_cells(P, O):
        Pc, Oc = P - P.mean(0), O - O.mean(0)
        r = (Pc * Oc).sum(0) / np.sqrt((Pc ** 2).sum(0) * (Oc ** 2).sum(0) + 1e-12)
        return float(np.nanmedian(r))
    sh = P[np.random.default_rng(0).permutation(len(P))]
    base = O.mean(0)
    rho_m = np.nanmedian([spearmanr(P[i], O[i]).correlation for i in range(len(P))])
    rho_b = np.nanmedian([spearmanr(base, O[i]).correlation for i in range(len(P))])
    print(f"[readout] {len(P)} cells x {len(panel)} genes | per-gene r across cells: model {across_cells(P, O):+.3f} "
          f"vs cell-shuffled {across_cells(sh, O):+.3f} | per-cell Spearman over genes: model {rho_m:+.3f}, "
          f"per-gene-mean predictor {rho_b:+.3f} (in-sample, favours the baseline)")


if __name__ == "__main__":
    main()

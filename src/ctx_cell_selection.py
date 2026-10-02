"""Record EXACTLY which cells the MaxToki extractions used, so every other model (scGPT, STATE, the expression-only
representation) can be run on identical cells, contexts and partitions.

Replicates ctx_extract_maxtoki.py pass 1 (same panels, same encoding via stream_panel, same >=250-cell filter, same
12 largest cell types in the same order, same seed-0 per-context shuffles) and writes, for CELLS_CTX = 600 and 1000,
the list of (panel file, row index, context index, partition) of the selected cells. Partition = position in the
shuffled per-context list modulo 2, exactly as the extractor assigns it.

Check: the token counts implied by the 600-cell selection must reproduce results/ctx217m600_L04.npz `counts`
(gene x context x partition occurrence counts, capped at 50) -- printed at the end.

Out: results/ctx_cell_selection.npz  (cells_600, cells_1000: int arrays [n, 4] = panel_idx, row, ctx_idx, partition;
     panels, contexts)
"""
import os, sys, json, collections
import numpy as np, h5py
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ctx_extract_maxtoki as EX
import ctx_tokenise as TK

SEED, NPART = 0, 2


def tokenised_cells(tok):
    """yield (panel_idx, row, cell_type, tokens) with the extractor's exact encoding and skip rules"""
    for pi, p in enumerate(EX.PANELS):
        path = os.path.join(EX.TS, p)
        with h5py.File(path, "r") as f:
            ens = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]).astype(str)
            ens = np.array([e.split(".")[0] for e in ens])
            ctg = f["obs"]["cell_type"]
            cats = np.array([x.decode() if isinstance(x, bytes) else x for x in ctg["categories"][:]]).astype(str)
            ctypes = cats[ctg["codes"][:]]
            X = TK.count_matrix(f); n = int(X.attrs["shape"][0]); indptr = X["indptr"][:]
            var_idx, token_ids, medians = tok.make_var_mapping(list(ens))
            pos = np.full(len(ens), -1, np.int64); pos[var_idx] = np.arange(len(var_idx))
            for r in range(n):
                s, e = int(indptr[r]), int(indptr[r + 1])
                idx, val = X["indices"][s:e], X["data"][s:e].astype(np.float32)
                keep = pos[idx] >= 0
                if not keep.any():
                    continue
                j = pos[idx[keep]]
                order = TK.rank_order(val[keep], medians[j])[: EX.MAX_LEN - 2]
                if not len(order):
                    continue
                yield pi, r, ctypes[r], token_ids[j[order]].astype(np.int64)


def main():
    from maxtoki_adapter import MaxTokiTokenizer
    tok = MaxTokiTokenizer(model_input_size=EX.MAX_LEN)
    cells = collections.defaultdict(list)                    # ct -> [(pi, row, tokens)] in panel order
    for pi, r, ct, toks in tokenised_cells(tok):
        cells[ct].append((pi, r, toks))
    ctx_names = sorted([c for c in cells if len(cells[c]) >= 250], key=lambda c: -len(cells[c]))[:EX.N_CTX]
    out = {"panels": np.array(EX.PANELS), "contexts": np.array(ctx_names)}
    for n_cells in (600, 1000):
        rng = np.random.default_rng(SEED)
        rows = []; sel_tokens = {}
        for ci, c in enumerate(ctx_names):
            lst = list(cells[c]); rng.shuffle(lst); lst = lst[:n_cells]
            for si, (pi, r, toks) in enumerate(lst):
                rows.append((pi, r, ci, si % NPART))
            sel_tokens[c] = [t for _, _, t in lst]
        out[f"cells_{n_cells}"] = np.array(rows, np.int64)
        if n_cells == 600:
            tok600 = sel_tokens
        print(f"[{n_cells}] {len(rows)} cells over {len(ctx_names)} contexts", flush=True)
    np.savez_compressed(os.path.join(HERE, "results", "ctx_cell_selection.npz"), **out)
    # consistency check against the 600-cell MaxToki extraction
    z = np.load(os.path.join(HERE, "results", "ctx217m600_L04.npz"), allow_pickle=True)
    tokmap = json.load(open(f"{EX.MSETUP}/token_dictionary.json")); ens2tid = {k: int(v) for k, v in tokmap.items()}
    panel_t = np.array([ens2tid.get(g, -1) for g in z["genes"].astype(str)]); gpos = {t: i for i, t in enumerate(panel_t)}
    cnt = np.zeros_like(z["counts"])
    for ci, c in enumerate(ctx_names):
        for si, toks in enumerate(tok600[c]):
            part = si % NPART
            for t in toks:
                gi = gpos.get(int(t))
                if gi is not None and cnt[part, ci, gi] < int(z["cap"]):
                    cnt[part, ci, gi] += 1
    same_ctx = list(z["contexts"].astype(str)) == ctx_names
    print(f"[check] contexts identical to ctx217m600: {same_ctx}; capped counts identical: {bool((cnt == z['counts']).all())} "
          f"(mismatching entries {int((cnt != z['counts']).sum())} of {cnt.size})", flush=True)


if __name__ == "__main__":
    main()

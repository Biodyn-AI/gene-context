"""Shared run-configuration helper so every analysis can run on any extraction prefix without touching the headline
outputs. With no environment variables set, every function reproduces the original behaviour exactly.

Environment:
  PREFIX     extraction prefix (default ctx_maxtoki)          TAPS / TAP   layers to analyse
  CELLS      cells per context of the extraction (default from a table / the npz)
  MIN_CTX    minimum count-balanced contexts per gene (default 9)
  COV        covariate for the rank control: rank (MaxToki mean rank on the same cells; default) | expr (mean log CP10k)
  OUTTAG     tag for the output file name (default: the prefix, plus non-default switches)
  RESEED_PER_TAP=1  restart the random generator for each tap (non-default runs only)
"""
import os, glob, json, pickle
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
NAME_ID = "/Volumes/Crucial X6/MacBook/Code/neuro-mechinterp/models/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl"
PREFIX = os.environ.get("PREFIX", "ctx_maxtoki")
IS_DEFAULT = PREFIX == "ctx_maxtoki" and not any(os.environ.get(k) for k in ("MIN_CTX", "COV", "OUTTAG", "CELLS", "TAPS",
                                                                            "TAP", "RESEED_PER_TAP"))
_CELLS_TABLE = {"ctx_maxtoki": 1000, "ctxscan": 1000, "ctx217m_cap20": 1000, "ctx217m600": 600, "ctx1b": 600,
                "ctxrand": 600, "ctxexpr": 600, "ctxexprpc": 600, "ctxexprcen": 600,
                "ctxexprcen50": 600, "ctxexprcen50c": 600}


def available_taps(prefix=None):
    p = prefix or PREFIX
    return sorted(int(os.path.basename(f)[len(p) + 2:len(p) + 4]) for f in glob.glob(os.path.join(RES, f"{p}_L[0-9][0-9].npz")))


def taps(default_list):
    if os.environ.get("TAPS"):
        t = [int(x) for x in os.environ["TAPS"].split(",")]
    elif PREFIX == "ctx_maxtoki":
        t = list(default_list)
    else:
        t = available_taps()
    missing = [x for x in t if not os.path.exists(os.path.join(RES, f"{PREFIX}_L{x:02d}.npz"))]
    if missing:
        raise SystemExit(f"{PREFIX}: requested taps {missing} not extracted; available {available_taps()}")
    return t


def tap(default):
    if os.environ.get("TAP"):
        return int(os.environ["TAP"])
    if PREFIX == "ctx_maxtoki":
        return default
    av = available_taps()
    if default in av:
        return default
    raise SystemExit(f"{PREFIX}: set TAP (available {av})")


def cells():
    if os.environ.get("CELLS"):
        return int(os.environ["CELLS"])
    for t in available_taps():
        z = np.load(os.path.join(RES, f"{PREFIX}_L{t:02d}.npz"), allow_pickle=True)
        if "cells_ctx" in z:
            return int(z["cells_ctx"])
        break
    if PREFIX in _CELLS_TABLE:
        return _CELLS_TABLE[PREFIX]
    raise SystemExit(f"{PREFIX}: unknown number of cells per context; set CELLS")


def min_ctx(default=9):
    return int(os.environ.get("MIN_CTX", default))


def cov():
    return os.environ.get("COV", "rank")


def npz_path(t, prefix=None):
    return os.path.join(RES, f"{prefix or PREFIX}_L{t:02d}.npz")


def load(t, prefix=None):
    """M (float32), counts, cap, genes (Ensembl, no version), contexts, full = count-balanced (context, gene)."""
    z = np.load(npz_path(t, prefix), allow_pickle=True)
    M = z["M"].astype(np.float32); counts = z["counts"]; cap = int(z["cap"])
    genes = np.array([g.split(".")[0] for g in z["genes"].astype(str)])
    if np.mean([g.startswith("ENSG") for g in genes]) < 0.9:
        sym2ens = {k.upper(): v for k, v in pickle.load(open(NAME_ID, "rb")).items()}
        print("[ctx_prefix] mapping gene symbols to Ensembl IDs", flush=True)
        genes = np.array([sym2ens.get(g.upper(), g) for g in genes])
    ctxs = (z["contexts"] if "contexts" in z else z["clusters"]).astype(str)
    full = (counts == cap).all(0)
    return M, counts, cap, genes, ctxs, full


def meta(t, prefix=None):
    """counts, cap, genes (Ensembl), contexts of an extraction, without loading the representation tensor M."""
    z = np.load(npz_path(t, prefix), allow_pickle=True)
    counts = z["counts"]; cap = int(z["cap"])
    genes = np.array([g.split(".")[0] for g in z["genes"].astype(str)])
    if np.mean([g.startswith("ENSG") for g in genes]) < 0.9:
        sym2ens = {k.upper(): v for k, v in pickle.load(open(NAME_ID, "rb")).items()}
        genes = np.array([sym2ens.get(g.upper(), g) for g in genes])
    ctxs = (z["contexts"] if "contexts" in z else z["clusters"]).astype(str)
    return counts, cap, genes, ctxs


def outtag():
    if os.environ.get("OUTTAG"):
        return os.environ["OUTTAG"]
    tag = PREFIX
    if PREFIX == "ctx_maxtoki" and (os.environ.get("TAPS") or os.environ.get("TAP")):   # never overwrite the headline
        tag += f"__taps{os.environ.get('TAPS') or os.environ.get('TAP')}".replace(",", "-")
    if os.environ.get("MIN_CTX"):
        tag += f"__minctx{os.environ['MIN_CTX']}"
    if os.environ.get("COV") and os.environ["COV"] != "rank":
        tag += f"__cov-{os.environ['COV']}"
    return tag


def out(name):
    """results/<name>.json for the default run; results/<name>__<tag>.json otherwise."""
    return os.path.join(RES, f"{name}.json" if IS_DEFAULT else f"{name}__{outtag()}.json")


def provenance(**extra):
    d = dict(prefix=PREFIX, cells=None, covariate=cov(), min_ctx=min_ctx(), outtag=outtag())
    try:
        d["cells"] = cells()
    except SystemExit:
        pass
    d.update(extra)
    return d


def rng_for(seed, tap_index, shared):
    """shared generator by default; a fresh one per tap when RESEED_PER_TAP=1 (non-default runs only)"""
    if os.environ.get("RESEED_PER_TAP") == "1" and not IS_DEFAULT:
        return np.random.default_rng(seed + 1000 * tap_index)
    return shared


def covariate_matrix(ctxs, genes):
    """(n_ctx, n_gene) covariate for the rank control on the extraction's own cells: MaxToki mean rank (COV=rank) or
    mean log1p CP10k over the cells (COV=expr). Cached in results/ctx_cov_cells{N}.npz."""
    import ctx_position_confound as CP
    N = cells()
    cache = os.path.join(RES, f"ctx_cov_cells{N}.npz")
    if not os.path.exists(cache):
        CP.build_covariate_cache(N, cache)
    c = np.load(cache, allow_pickle=True)
    cctx = list(c["contexts"].astype(str)); cg = {g: i for i, g in enumerate(c["genes"].astype(str))}
    key = "mean_rank" if cov() == "rank" else "mean_logcp10k_all"
    src = c[key]
    out_ = np.full((len(ctxs), len(genes)), np.nan)
    for ci, ct in enumerate(ctxs):
        if ct not in cctx:
            continue
        k = cctx.index(ct)
        for gi, g in enumerate(genes):
            j = cg.get(g)
            if j is not None:
                out_[ci, gi] = src[k, j]
    return out_

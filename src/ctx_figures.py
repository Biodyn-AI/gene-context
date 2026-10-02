"""Publication figures for the contextual-gene-representation paper. Reads results/ctx_*.json only.
Saved to $FIGDIR (default figures/) as ctx_fig{1..5}.pdf.

F1 contextualisation by layer: full depth scan (EXCESS same vs different gene) + the main extraction's layers
F2 rank control (rank-stable vs rank-moving vs rank-residualised EXCESS)
F3 functional-axis interaction power vs random / co-expression-matched / strong co-expression nulls
F4 causal steering dose-response (signed swing vs norm-matched random push)
F5 cross-model + scaling + random-weights control (EXCESS and functional-z)

Journal layout (BMC): single-column figures are exactly 85 mm wide, two-column figures 170 mm; all text >= 7.5 pt;
TrueType fonts embedded; FIG_TITLES=0 leaves titles out of the graphic (they go in the manuscript).
"""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
if not os.path.isdir(RES):                       # in the public repository the scripts sit in src/, results in ../results
    RES = os.path.join(os.path.dirname(HERE), "results")
FIG = os.environ.get("FIGDIR", os.path.join(os.path.dirname(RES), "figures")); os.makedirs(FIG, exist_ok=True)
TITLES = os.environ.get("FIG_TITLES", "1") == "1"      # FIG_TITLES=0: no titles inside the graphic (BMC style)
MM = 1 / 25.4
W1, W2 = 85 * MM, 170 * MM                               # BMC single- and double-column widths (inches)
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "legend.fontsize": 7.5, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "pdf.fonttype": 42, "ps.fonttype": 42})   # TrueType fonts embedded
BLUE, RED, GREY, GREEN, LBLUE = "#2c6fbb", "#c0392b", "#95a5a6", "#e69f00", "#a9cbe8"   # orange, not green, next to red
ORANGE = GREEN
ANN = 7.5                                                # minimum annotation size (pt)


def L(f): return json.load(open(os.path.join(RES, f)))


def title(ax, *a, **k):
    if TITLES:
        ax.set_title(*a, **k)


def save(fig, name):
    # constrained layout keeps the exact figure width (no bbox="tight", which would change it)
    fig.savefig(os.path.join(FIG, name)); plt.close(fig); print(name)


def fig1():
    # full depth scan (ctx_layer_curve.py on the ctxscan extraction, 4000-gene panel) + the four layers of the
    # main 6000-gene extraction (ctx_polysemy.json) as open markers
    sc = L("ctx_layer_curve.json")["layers"]
    ls = sorted(int(k[1:]) for k in sc)
    same = [sc[f"L{l:02d}"]["same"] for l in ls]; diff = [sc[f"L{l:02d}"]["diff"] for l in ls]
    d = L("ctx_polysemy.json")["taps"]
    lm = sorted(int(k[1:]) for k in d)
    fig, ax = plt.subplots(figsize=(W1, 2.6), layout="constrained")
    ax.plot(ls, same, "-o", color=BLUE, ms=3.5, label="same gene, independent cells")
    ax.plot(ls, diff, "-s", color=GREY, ms=3.5, label="different gene (null)")
    ax.plot(lm, [d[f"L{l:02d}"]["same"] for l in lm], "D", mfc="none", color="k", ms=5.5,
            label="main extraction")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("layer (0 = embedding output)"); ax.set_ylabel("agreement of the gene's\ncontext shift (cosine)")
    ax.set_xticks(ls)
    title(ax, "Gene-specific context response by layer", fontsize=8.5)
    ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.12))
    save(fig, "ctx_fig1.pdf")


def fig2():
    d = L("ctx_position_confound.json")["taps"]
    layers = sorted(int(k[1:]) for k in d)
    st = [d[f"L{l:02d}"]["excess_rank_stable"] for l in layers]
    mv = [d[f"L{l:02d}"]["excess_rank_moving"] for l in layers]
    rs = [d[f"L{l:02d}"]["excess_residualised"] for l in layers]
    x = np.arange(len(layers)); w = 0.26
    fig, ax = plt.subplots(figsize=(W1, 2.6), layout="constrained")
    ax.bar(x - w, st, w, color=BLUE, label="rank-stable genes")
    ax.bar(x, mv, w, color=RED, label="rank-moving genes")
    ax.bar(x + w, rs, w, color=GREEN, label="after rank regression")
    ax.set_xticks(x); ax.set_xticklabels([str(l) for l in layers]); ax.set_xlabel("layer")
    ax.set_ylabel("EXCESS (gene-specific\ncontext response)")
    ax.set_ylim(0, max(st + mv + rs) * 1.35)
    title(ax, "The effect is not the gene's rank position", fontsize=8.5)
    ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.02), columnspacing=0.8,
              handlelength=1.2)
    save(fig, "ctx_fig2.pdf")


def fig3():
    """Per functional axis: the axis's rank-controlled reproducible interaction power (red marker + line) against
    three null distributions with the functional poles' sizes: random genes (v2), gene sets matched on size AND
    average co-expression (v3), and modules of the genes most co-expressed with a random seed gene (v2 null b)."""
    d2 = L("ctx_coexpr_null_v2.json")["axes"]; d3 = L("ctx_coexpr_null_v3.json")["axes"]
    nice = {"nuclear_vs_surface": "nuclear / surface", "mito_vs_cytoskeleton": "mitochondrion / cytoskeleton",
            "transcription_vs_transport": "transcription / transport"}
    order = ["nuclear_vs_surface", "mito_vs_cytoskeleton", "transcription_vs_transport"]
    cols = [GREY, LBLUE, BLUE]
    fig, ax = plt.subplots(figsize=(W2, 3.3), layout="constrained")
    for i, name in enumerate(order):
        r2, r3 = d2[name], d3[name]
        data = [np.array(r2["null_random"]), np.array(r3["null_matched"]), np.array(r2["null_indep_coexpr"])]
        bp = ax.boxplot(data, positions=[i - 0.15, i + 0.1, i + 0.35], widths=0.2, patch_artist=True,
                        showfliers=False, medianprops=dict(color="k"))
        for patch, c in zip(bp["boxes"], cols):
            patch.set_facecolor(c); patch.set_alpha(0.9)
        pw = r3["power"]
        ax.plot([i - 0.3, i + 0.48], [pw, pw], color=RED, lw=1.2, zorder=5)
        ax.scatter([i - 0.38], [pw], s=40, color=RED, zorder=6, edgecolor="k", lw=0.6)
        pm_txt = f"p < 1/{r3['matched_n']}" if r3["matched_n_ge"] == 0 else f"p = {r3['matched_p']:.3f}"
        ps = r2["indep_coexpr_p"]
        ps_txt = f"p < 1/{len(r2['null_indep_coexpr'])}" if ps == 0 else f"p = {ps:.2f}"
        top = max([pw] + [float(np.max(x)) for x in data])
        ax.text(i + 0.1, top * 1.04 + 0.05, f"vs matched: {pm_txt}\nvs strong: {ps_txt}", ha="center", va="bottom",
                fontsize=ANN)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([nice[n] for n in order])
    ax.set_xlim(-0.6, len(order) - 0.4)
    ymax = max(max(float(np.max(d2[n]["null_indep_coexpr"])), float(np.max(d3[n]["null_matched"])), d3[n]["power"])
               for n in order)
    ax.set_ylim(0, ymax * 1.42)
    ax.set_ylabel("reproducible interaction power\n(rank-controlled)")
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    ax.legend(handles=[Patch(facecolor=cols[0], label="random genes"),
                       Patch(facecolor=cols[1], label="matched on size and average co-expression"),
                       Patch(facecolor=cols[2], label="most co-expressed genes (strong null)"),
                       Line2D([0], [0], marker="o", color=RED, markerfacecolor=RED, markeredgecolor="k",
                              markersize=5, label="functional axis")],
              frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.0))
    title(ax, "Functional organisation against random and co-expression nulls", fontsize=8.5)
    save(fig, "ctx_fig3.pdf")


def fig4():
    d = L("ctx_causal.json")
    al = d["alphas_xResidNorm"]
    fs = [d["signed"][f"alpha_{a}"]["func_swing"] for a in al]
    se = [d["signed"][f"alpha_{a}"]["func_swing_sem"] for a in al]
    rs = [d["signed"][f"alpha_{a}"]["rand_swing"] for a in al]
    fig, ax = plt.subplots(figsize=(W1, 2.6), layout="constrained")
    ax.errorbar(al, fs, yerr=se, fmt="-o", color=BLUE, ms=3.5, label="functional direction", capsize=2.5)
    ax.plot(al, rs, "-s", color=GREY, ms=3.5, label="norm-matched random direction")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(al); ax.set_xticklabels([f"{a:g}" for a in al])
    ax.set_xlabel("steering strength (× mean residual norm)")
    ax.set_ylabel("swing in nuclear − surface logit\nat untouched genes")
    title(ax, "The functional-context direction is causally used", fontsize=8.5)
    ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.1))
    save(fig, "ctx_fig4.pdf")


def fig5():
    """All representations on the same 600 cells per cell type, gene panel and cap (results/ctx_cross_model__all600.json);
    each model at its depth-matched layer. Falls back to the older mixed comparison if that file is absent."""
    p = os.path.join(RES, "ctx_cross_model__all600.json")
    if not os.path.exists(p) and os.environ.get("FIG5_LEGACY") != "1":
        raise SystemExit("fig5: results/ctx_cross_model__all600.json missing (set FIG5_LEGACY=1 for the old figure)")
    if os.path.exists(p):
        d = json.load(open(p))
        order = [("MaxToki-217M", "MaxToki-217M (layer 4)", BLUE), ("MaxToki-1B-L7", "MaxToki-1B (layer 7)", BLUE),
                 ("scGPT-L4", "scGPT (layer 4)", BLUE), ("STATE-L6", "STATE (layer 6)", BLUE),
                 ("MaxToki-217M-random", "MaxToki-217M, untrained", GREY), ("scGPT-random", "scGPT, untrained", GREY),
                 ("STATE-random", "STATE, untrained", GREY),
                 ("expression-landmarks", "Expression only (landmarks)", ORANGE),
                 ("expression-PCs", "Expression only (PCs)", ORANGE),
                 ("expression-centroids", "Expression only (centroids)", ORANGE)]
    else:
        d = L("ctx_cross_model.json")
        order = [("scGPT", "scGPT", BLUE), ("STATE-SE", "STATE (SE-600M)", BLUE), ("MaxToki-217M", "MaxToki-217M\n(layer 4)", BLUE),
                 ("MaxToki-1B-L7", "MaxToki-1B\n(layer 7)", BLUE), ("MaxToki-217M-random", "MaxToki-217M\n(random weights)", GREY)]
    rows = [(lab, d[k], c) for k, lab, c in order if k in d and isinstance(d[k], dict) and "excess" in d[k]]
    names = [r[0] for r in rows]; cols = [r[2] for r in rows]
    ex = [r[1]["excess"] for r in rows]
    fz = [r[1]["func_z"]["nuclear_vs_surface"]["z"] if r[1]["func_z"].get("nuclear_vs_surface") else 0 for r in rows]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 0.32 * len(rows) + 1.0), layout="constrained", sharey=True)
    y = np.arange(len(names))[::-1]
    a1.barh(y, ex, color=cols); a1.set_yticks(y); a1.set_yticklabels(names)
    a1.set_xlabel("EXCESS (gene-specific context shift)")
    title(a1, "Contextualisation strength", fontsize=8.5)
    a2.barh(y, fz, color=cols); a2.set_xlabel("functional-z (nuclear/surface axis)")
    a2.axvline(0, color="k", lw=0.5)
    a2.axvline(3, color="k", ls="--", lw=0.7); a2.text(3.3, y.min() - 0.3, "z = 3", fontsize=ANN, va="center")
    a2.set_ylim(y.min() - 0.7, y.max() + 0.5)
    title(a2, "Functional organisation", fontsize=8.5)
    save(fig, "ctx_fig5.pdf")


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4", "5"]
    fns = {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5}
    for k in which:
        try: fns[k]()
        except Exception as e: print(f"F{k} ERR {repr(e)[:160]}")
    print(f"[done] -> {FIG}")

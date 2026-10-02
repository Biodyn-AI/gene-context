"""Collect every number quoted in the manuscript from results/ctx_*.json (and the extraction counts), format it, and
write build/numbers.json. fill_template.py then renders manuscript_bmc.template.md -> manuscript_bmc.md.
A key whose source is missing is simply absent, so the fill step fails loudly on any placeholder it cannot resolve."""
import os, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RG = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.environ.get("NUMBERS_RES", os.path.join(RG, "results"))
N = {}
AX = {"nuclear_vs_surface": "nuc", "mito_vs_cytoskeleton": "mito", "transcription_vs_transport": "trans",
      "cellcycle_vs_diff": "cc"}


def J(name):
    p = os.path.join(RES, name)
    return json.load(open(p)) if os.path.exists(p) else None


def sg(x, d=3):            # signed fixed
    return f"{x:+.{d}f}".replace("-", "−")


def f(x, d=3):             # unsigned fixed
    return f"{x:.{d}f}".replace("-", "−")


def pfmt(p, n=None):       # Monte Carlo p: show "< 1/(n+1)" style only when it is the floor
    if n is not None and abs(p - 1 / (n + 1)) < 1e-12:
        return f"1/{n + 1}"
    return f"{p:.3f}" if p >= 0.001 else f"{p:.1e}".replace("e-0", " × 10^−").replace("e-", " × 10^−") + "^"


def sci(p):                # 7.2e-04 -> 7 × 10^−4^ (pandoc superscript)
    m, e = f"{p:.1e}".split("e")
    m = m.rstrip("0").rstrip(".") if "." in m else m
    return f"{float(m):.0f} × 10^−{int(e[1:])}^" if abs(float(m) - round(float(m))) < 0.05 else f"{m} × 10^−{int(e[1:])}^"


def put(k, v):
    N[k] = v


def rel(x, ref):
    return f"{100 * x / ref:.0f}"


# ---------------- extraction-level counts ----------------
# These few numbers come from the large .npz tensors, which are not committed; when the tensors are present they are
# also written to results/ctx_panel_stats.json, and a fresh clone without the tensors reads them from there.
PANEL_KEYS = ("n_genes_all_ctx", "n_genes_pair_median", "n_genes_panel", "state_panel_genes", "scgpt_panel_genes")


def _panel_stats_fallback():
    p = os.path.join(RES, "ctx_panel_stats.json")
    if os.path.exists(p):
        for k, v in json.load(open(p)).items():
            if k not in N: put(k, v)


def _panel_stats_save():
    have = {k: N[k] for k in PANEL_KEYS if k in N}
    if len(have) == len(PANEL_KEYS):
        json.dump(have, open(os.path.join(RES, "ctx_panel_stats.json"), "w"), indent=1)


def _sec0():
    global N
    z = os.path.join(RES, "ctx_maxtoki_L04.npz")
    if os.path.exists(z):
        d = np.load(z, allow_pickle=True)
        full = (d["counts"] == int(d["cap"])).all(0)          # (ctx, gene)
        put("n_genes_all_ctx", str(int(full.all(0).sum())))
        nC = full.shape[0]
        pair = [int((full[i] & full[j]).sum()) for i in range(nC) for j in range(i + 1, nC)]
        put("n_genes_pair_median", str(int(np.median(pair))))
        put("n_genes_panel", str(int(full.shape[1])))


try:
    _sec0()
except Exception as e:
    print('section extraction-level failed:', repr(e)[:150])


# ---------------- EXCESS by layer ----------------
def _sec1():
    global N
    pol = J("ctx_polysemy.json")
    if pol:
        for L, r in pol["taps"].items():
            l = int(L[1:])
            put(f"excess_L{l}", sg(r["excess"])); put(f"same_L{l}", sg(r["same"])); put(f"diff_L{l}", sg(r["diff"]))
            put(f"mainrep_L{l}", sg(r["main_effect_replication"])); put(f"npairs_L{l}", str(r["n_pairs"]))
            put(f"excess_ci_L{l}", f"{sg(r['ci'][0])} to {sg(r['ci'][1])}")
            put(f"excess_ci_halfwidth_L{l}", f"{(r['ci'][1] - r['ci'][0]) / 2:.3f}")
            if r.get("pairs"):
                put(f"pair_genes_median_L{l}", str(int(np.median([p["n"] for p in r["pairs"]]))))
    lc = J("ctx_layer_curve.json")
    if lc:
        lay = lc["layers"]; ex = {int(k[1:]): v["excess"] for k, v in lay.items()}
        pk = max(ex, key=ex.get)
        put("scan_peak_layer", str(pk)); put("scan_peak_excess", sg(ex[pk]))
        for l, v in ex.items():
            put(f"scan_L{l}", sg(v)); put(f"scan_pct_of_peak_L{l}", rel(v, ex[pk]))
        put("scan_last_layer", str(max(ex))); put("scan_last_excess", sg(ex[max(ex)]))
        put("scan_mainrep_min", f"{min(v['main_effect_replication'] for k, v in lay.items() if int(k[1:]) > 0):.3f}")


try:
    _sec1()
except Exception as e:
    print('section EXCESS failed:', repr(e)[:150])


# ---------------- anisotropy ----------------
def _sec2():
    global N
    an = J("ctx_anisotropy.json")
    if an:
        lay = an["layers"]; Ls = sorted(k for k in lay if k != "L00")
        put("aniso_raw_range", f"{sg(min(lay[L]['anisotropy_raw'] for L in Ls), 2)} to {sg(max(lay[L]['anisotropy_raw'] for L in Ls), 2)}")
        put("aniso_layers", f"{int(Ls[0][1:])}–{int(Ls[-1][1:])}")
        put("selfsim_raw_seq", " → ".join(sg(lay[L]["selfsim_raw"], 2) for L in Ls))
        put("selfsim_z_seq", "/".join(sg(lay[L]["corrected_selfsim_zscored"], 2) for L in Ls))
        put("selfsim_z_list", ", ".join(sg(lay[L]["corrected_selfsim_zscored"], 2) for L in Ls[:-1]) + " and " + sg(lay[Ls[-1]]["corrected_selfsim_zscored"], 2))
        put("selfsim_raw_list", ", ".join(sg(lay[L]["selfsim_raw"], 2) for L in Ls[:-1]) + " and " + sg(lay[Ls[-1]]["selfsim_raw"], 2))
        if all("corrected_selfsim_raw" in lay[L] for L in Ls):
            put("selfsim_rawcorr_list", ", ".join(sg(lay[L]["corrected_selfsim_raw"], 2) for L in Ls[:-1]) + " and " + sg(lay[Ls[-1]]["corrected_selfsim_raw"], 2))
        put("aniso_z_max", f"{max(abs(lay[L]['anisotropy_zscored']) for L in Ls):.2f}")
        put("toppc_raw_last", f"{lay[Ls[-1]]['top_pc_raw']:.2f}"); put("toppc_z_last", f"{lay[Ls[-1]]['top_pc_zscored']:.2f}")


try:
    _sec2()
except Exception as e:
    print('section anisotropy failed:', repr(e)[:150])


# ---------------- rank control ----------------
def _sec3():
    global N
    pc = J("ctx_position_confound.json"); pol = J("ctx_polysemy.json")
    if pc and pol:
        for L, r in pc["taps"].items():
            l = int(L[1:])
            put(f"rank_stable_L{l}", sg(r["excess_rank_stable"])); put(f"rank_moving_L{l}", sg(r["excess_rank_moving"]))
            put(f"rank_resid_L{l}", sg(r["excess_residualised"])); put(f"rank_rho_L{l}", f"{r['rho_mag_vs_rank']:.2f}")
            if L in pol["taps"]:
                put(f"rank_survive_L{l}", rel(r["excess_residualised"], pol["taps"][L]["excess"]))
            if "excess_residualised_curved" in r:
                put(f"rank_curved_L{l}", sg(r["excess_residualised_curved"]))
                if L in pol["taps"]: put(f"rank_curved_survive_L{l}", rel(r["excess_residualised_curved"], pol["taps"][L]["excess"]))
                put(f"rank_matchdiff_L{l}", sg(r["diff_rank_matched"])); put(f"rank_vs_matched_L{l}", sg(r["excess_vs_rank_matched"]))
                put(f"rank_med_stable_L{l}", f"{r['median_abs_rank_change_stable']:.0f}"); put(f"rank_med_moving_L{l}", f"{r['median_abs_rank_change_moving']:.0f}")


try:
    _sec3()
except Exception as e:
    print('section rank failed:', repr(e)[:150])


# ---------------- independent dataset ----------------
def _sec4():
    global N
    ind = J("ctx_independent.json")
    if ind:
        for L, r in ind["taps"].items():
            l = int(L[1:])
            put(f"setty_excess_L{l}", sg(r["excess"])); put(f"setty_diff_L{l}", sg(r["diff_null"], 4))
            put(f"setty_ci_L{l}", f"{sg(r['excess_ci'][0])} to {sg(r['excess_ci'][1])}")
            fz = r["func_z"].get("nuclear_vs_surface")
            if fz: put(f"setty_fz_L{l}", sg(fz["z"], 1))
            put("setty_nctx", str(r["n_ctx"]))
    ind2 = J("ctx_independent_ctx_devel_spliced.json")
    if ind2:
        for L, r in ind2["taps"].items():
            l = int(L[1:]); put(f"settysp_excess_L{l}", sg(r["excess"]))
            fz = r["func_z"].get("nuclear_vs_surface")
            if fz: put(f"settysp_fz_L{l}", sg(fz["z"], 1))


try:
    _sec4()
except Exception as e:
    print('section independent failed:', repr(e)[:150])


# ---------------- functional axes ----------------
def _sec5():
    global N
    fa = J("ctx_functional_axes.json")
    if fa:
        for L, axes in fa["taps"].items():
            l = int(L[1:])
            for a, r in axes.items():
                k = AX[a]
                put(f"fz_raw_{k}_L{l}", sg(r["z_raw"], 1)); put(f"fz_rc_{k}_L{l}", sg(r["z_rank_controlled"], 1))
                put(f"auc_{k}_L{l}", f"{r['validity_auc']:.2f}")
                put(f"poles_{k}", f"{r['poleA']}/{r['poleB']}")
    lo = J("ctx_celltype_loadings.json")
    if lo:
        ld = lo["per_context_loading"]
        for ct, key in [("CD8-positive, alpha-beta T cell", "cd8"), ("CD4-positive, alpha-beta T cell", "cd4"),
                        ("B cell", "bcell"), ("pulmonary alveolar type 2 cell", "at2"),
                        ("pulmonary alveolar type 1 cell", "at1"), ("macrophage", "macro"), ("neutrophil", "neut")]:
            if ct in ld: put(f"load_{key}", sg(ld[ct], 1))
        order = sorted(ld, key=lambda c: -ld[c]); put("load_top3", ", ".join(order[:3])); put("load_bottom3", ", ".join(order[-3:]))


try:
    _sec5()
except Exception as e:
    print('section functional failed:', repr(e)[:150])


# ---------------- co-expression / tightness nulls ----------------
def _sec6():
    global N
    v2 = J("ctx_coexpr_null_v2.json"); v3 = J("ctx_coexpr_null_v3.json"); ti = J("ctx_tightness_null.json")
    if v2 and v3:
        for a, r3 in v3["axes"].items():
            k = AX[a]; r2 = v2["axes"][a]
            put(f"pow_{k}", f"{r3['power']:.2f}")
            put(f"coh_{k}_A", f"{r3['coh_A']:.3f}"); put(f"coh_{k}_B", f"{r3['coh_B']:.3f}")
            put(f"cohrand_{k}", f"{r3['random_set_cohA_mean']:.3f}")
            put(f"strongcoh_{k}_A", f"{r3['strong_module_cohA_mean']:.3f}"); put(f"strongcoh_{k}_B", f"{r3['strong_module_cohB_mean']:.3f}")
            put(f"pm_{k}", pfmt(r3["matched_p"], r3["matched_n"])); put(f"pm_nge_{k}", f"{r3['matched_n_ge']}/{r3['matched_n']}")
            put(f"pm_mean_{k}", f"{r3['matched_mean']:.2f}"); put(f"pm_p95_{k}", f"{r3['matched_p95']:.2f}")
            put(f"pm_below_{k}", f"{r3['matched_poles_below_random']}")
            ns = len(r2["null_indep_coexpr"]); nge_s = int(sum(x >= r2["power"] for x in r2["null_indep_coexpr"]))
            put(f"ps_{k}", pfmt(r2["indep_coexpr_p"], ns)); put(f"ps_nge_{k}", f"{nge_s}/{ns}")
            put(f"ps_mean_{k}", f"{r2['indep_coexpr_mean']:.2f}"); put(f"ps_p95_{k}", f"{np.percentile(r2['null_indep_coexpr'], 95):.2f}")
            put(f"ps_pct_{k}", f"{100 * nge_s / ns:.0f}")
            put(f"pr_{k}", pfmt(r2["random_p"], len(r2["null_random"]))); put(f"pr_mean_{k}", f"{r2['random_mean']:.2f}")
            cur = r3["curve"]; fr = sorted({c["f"] for c in cur})
            put(f"curve_max_{k}", f"{max(c['power'] for c in cur):.2f}")
            put(f"curve_full_mean_{k}", f"{np.mean([c['power'] for c in cur if c['f'] == 1.0]):.2f}")
            put(f"curve_nge_{k}", str(sum(c["power"] >= r3["power"] for c in cur)))
            put(f"curve_n_{k}", str(len(cur)))
        rat = [r3["matched_p95"] / r3["matched_mean"] for r3 in v3["axes"].values() if r3.get("matched_mean")]
        if rat: put("pm_p95_over_mean_range", f"{min(rat):.1f}–{max(rat):.1f}")
        rng_ = v3.get("strong_over_functional_coherence_range")
        if rng_: put("strong_ratio_range", f"{rng_[0]:.0f}–{rng_[1]:.0f}" if rng_[0] >= 1.5 else f"{rng_[0]:.1f}–{rng_[1]:.1f}")
    if ti:
        for L, tv in ti["taps"].items():
            l = int(L[1:])
            for a, r in tv["axes"].items():
                k = AX[a]
                put(f"pt_{k}_L{l}", pfmt(r["matched_p"], r["matched_n"])); put(f"pt_nge_{k}_L{l}", f"{r['matched_n_ge']}/{r['matched_n']}")
                put(f"pt_mean_{k}_L{l}", f"{r['matched_mean']:.2f}"); put(f"tight_{k}_L{l}", f"{r['tightness_A']:.3f}/{r['tightness_B']:.3f}")
                put(f"tightrand_{k}_L{l}", f"{r['random_set_tightness']:.3f}")


try:
    _sec6()
except Exception as e:
    print('section co-expression failed:', repr(e)[:150])


# ---------------- causal ----------------
def _sec7():
    global N
    for fn, tag in [("ctx_causal.json", "nucL4"), ("ctx_causal_nuc_surf_L08.json", "nucL8"),
                    ("ctx_causal_mito_cyto_L04.json", "mitoL4"), ("ctx_causal_trans_transport_L04.json", "transL4")]:
        c = J(fn)
        if not c: continue
        m = c["signed"]["alpha_0.5"]
        put(f"cz_{tag}_swing", sg(m["func_swing"])); put(f"cz_{tag}_rand", sg(m["rand_swing"]))
        put(f"cz_{tag}_spp", sg(m["spec_plus"])); put(f"cz_{tag}_spm", sg(m["spec_minus"]))
        put(f"cz_{tag}_k", f"{m['func_gt_rand']} of {m['n']}"); put(f"cz_{tag}_p", sci(m["sign_p"]) if m["sign_p"] < 1e-3 else f"{m['sign_p']:.2g}")
        put(f"cz_{tag}_p_raw", m["sign_p"])
        put(f"cz_{tag}_swings", ", ".join(sg(c["signed"][f"alpha_{a}"]["func_swing"]) for a in c["alphas_xResidNorm"]))
        # two-sided sign-test p for a reversal (func swing < random swing)
        n = m["n"]; kk = n - m["func_gt_rand"]
        p_rev = sum(math.comb(n, i) for i in range(kk, n + 1)) / 2 ** n
        pr2 = min(1, 2 * p_rev); put(f"cz_{tag}_p_rev2", sci(pr2) if pr2 < 1e-3 else f"{pr2:.2g}")
        nat = c.get("natural_movement_along_axis")
        if nat:
            put(f"cz_{tag}_nat_sd", f"{nat['median_sd_across_contexts']:.2f}"); put(f"cz_{tag}_nat_range", f"{nat['median_range_across_contexts']:.2f}")
            put(f"cz_{tag}_nat_range95", f"{nat['p95_range_across_contexts']:.2f}")
            put(f"cz_{tag}_pole_dist", f"{nat['pole_centroid_distance']:.2f}")
            put(f"cz_{tag}_push", ", ".join(f"{x:.1f}" for x in c["push_size_raw"]))
            put(f"cz_{tag}_push_mid_over_range", f"{c['push_size_raw'][1] / nat['median_range_across_contexts']:.0f}")
            put(f"cz_{tag}_resid_norm", f"{c['resid_norm']:.1f}")
    ps_ = [J(f) for f in ["ctx_causal.json", "ctx_causal_nuc_surf_L08.json", "ctx_causal_mito_cyto_L04.json", "ctx_causal_trans_transport_L04.json"]]


try:
    _sec7()
except Exception as e:
    print('section causal failed:', repr(e)[:150])


# ---------------- Level 2 ----------------
def _sec8():
    global N
    dp = J("ctx_directional_probe.json")
    if dp:
        for L, tv in dp["taps"].items():
            l = int(L[1:])
            for a, r in tv["axes"].items():
                k = AX[a]
                put(f"beta_{k}_L{l}", sg(r["beta"], 2)); put(f"beta_p_{k}_L{l}", pfmt(r["p"], r["null_n"]))
                put(f"beta_nullmean_{k}_L{l}", sg(r["null_mean"], 3)); put(f"beta_nullsd_{k}_L{l}", f"{r['null_sd']:.2f}")
                put(f"beta_nullp95_{k}_L{l}", sg(r["null_p95"], 2))
        allp = [(r["p"], r["null_n"]) for tv in dp["taps"].values() for r in tv["axes"].values()]
        put("beta_min_p", pfmt(*min(allp))); put("beta_ntests", str(dp.get("holm_n_tests", len(allp))))
        put("beta_holm_sig", str(len(dp.get("holm_significant", []))))
        put("beta_nullsd_range", f"{min(r['null_sd'] for tv in dp['taps'].values() for r in tv['axes'].values()):.2f}–"
                                 f"{max(r['null_sd'] for tv in dp['taps'].values() for r in tv['axes'].values()):.2f}")
    ct = J("ctx_curated_targets.json")
    if ct:
        r = ct["taps"]["L04"]
        put("tf_n", str(r["n_tfs"])); put("tf_static", sg(r["static_excess_mean"])); put("tf_static_pos", f"{round(r['static_frac_pos'] * r['n_tfs'])} of {r['n_tfs']}")
        put("tf_static_p_2dp", f"{r['static_sign_p']:.2f}"); put("tf_mod", sg(r["modulation_mean_rho"]))
        put("tf_mod_pos", f"{round(r['modulation_frac_pos'] * r['n_modul'])} of {r['n_modul']}"); put("tf_mod_p", f"{r['modulation_sign_p']:.2f}")
    sw = J("ctx_switcher_test.json")
    if sw:
        put("switch_on_panel", ", ".join(sw.get("switchers_on_panel", []))); put("switch_n_on_panel", str(len(sw.get("switchers_on_panel", []))))


try:
    _sec8()
except Exception as e:
    print('section Level failed:', repr(e)[:150])


# ---------------- cross-model ----------------
def _sec9():
    global N
    cm = J("ctx_cross_model.json")
    if cm:
        key = {"scGPT": "scgpt", "STATE-SE": "state", "MaxToki-217M": "m217", "MaxToki-1B": "m1b", "MaxToki-217M-1k": "m217k",
               "MaxToki-217M-random": "rand", "MaxToki-217M-cap20": "cap20", "MaxToki-217M-L2": "m217L2", "MaxToki-1B-L7": "m1bL7"}
        for lab, k in key.items():
            r = cm.get(lab)
            if not r or "error" in r: continue
            put(f"cm_{k}_excess", sg(r["excess"])); put(f"cm_{k}_ctx", str(r["n_ctx"])); put(f"cm_{k}_genes", str(r["n_genes"]))
            put(f"cm_{k}_dim", str(r["dim"])); put(f"cm_{k}_mainrep", f"{r['main_effect_replication']:.3f}")
            put(f"cm_{k}_pairs", str(r.get("pairs_scored", "")))
            fz = r["func_z"].get("nuclear_vs_surface")
            if fz:
                put(f"cm_{k}_fz", sg(fz["z"], 1)); put(f"cm_{k}_poles", f"{fz['poleA']}/{fz['poleB']}")
        ok = [cm[l] for l in ("scGPT", "STATE-SE", "MaxToki-217M", "MaxToki-1B") if l in cm and "error" not in cm[l]]
        if ok:
            put("cm_excess_range", f"{min(r['excess'] for r in ok):.2f}–{max(r['excess'] for r in ok):.2f}")
            put("cm_mainrep_min", f"{min(r['main_effect_replication'] for r in ok):.2f}")


try:
    _sec9()
except Exception as e:
    print('section cross-model failed:', repr(e)[:150])


# ---------------- prediction link ----------------
def _sec10():
    global N
    pl = J("ctx_prediction_link.json")
    if pl:
        put("pl_benefit", f"{pl['mean_context_benefit']:+.1f}".replace("-", "−")); put("pl_benefit_median", f"{pl.get('median_context_benefit', float('nan')):+.1f}")
        put("pl_fold", f"{math.exp(pl['mean_context_benefit']):,.0f}".replace(",", " "))
        put("pl_ngenes", str(pl["n_genes"])); put("pl_rho", sg(pl["partial_rho"])); put("pl_ci", f"{sg(pl['partial_ci'][0])} to {sg(pl['partial_ci'][1])}")
        put("pl_skipped", str(pl.get("n_positions_skipped_no_clean_donor", ""))); put("pl_ncomp", str(pl.get("n_comparisons", "")))


try:
    _sec10()
except Exception as e:
    print('section prediction failed:', repr(e)[:150])


# ---------------- multi-model: the full pipeline on every 600-cell representation (2 Oct 2026) ----------------
MM = [("m217", "ctx217m600", 4), ("m1b", "ctx1b", 7), ("scgpt", "ctx_scgpt600", 4), ("state", "ctx_state600", 6),
      ("exprL", "ctxexpr", 0), ("exprP", "ctxexprpc", 0), ("exprC", "ctxexprcen50c", 0), ("exprCall", "ctxexprcen50", 0),
      ("exprC512", "ctxexprcen", 0)]


def _sec11():
    global N
    zs = os.path.join(RES, "ctx_state600_L06.npz")
    if os.path.exists(zs):
        put("state_panel_genes", str(len(np.load(zs, allow_pickle=True)['genes'])))
    zc = os.path.join(RES, "ctx_scgpt600_L04.npz")
    if os.path.exists(zc):
        put("scgpt_panel_genes", str(int((np.load(zc, allow_pickle=True)["counts"].sum((0, 1)) > 0).sum())))
    for k, pre, tap in MM:
        L = f"L{tap:02d}"
        pol = J(f"ctx_polysemy__{pre}.json")
        if pol and L in pol["taps"]:
            r = pol["taps"][L]
            put(f"mm_{k}_excess", sg(r["excess"])); put(f"mm_{k}_same", sg(r["same"])); put(f"mm_{k}_diff", sg(r["diff"]))
            put(f"mm_{k}_npairs", str(r["n_pairs"])); put(f"mm_{k}_mainrep", sg(r["main_effect_replication"], 2))
        pc = J(f"ctx_position_confound__{pre}.json")
        if pc and pol and L in pc["taps"]:
            put(f"mm_{k}_rank_survive", rel(pc["taps"][L]["excess_residualised"], pol["taps"][L]["excess"]))
        for mc in ("6", "9"):
            sfx = f"__{pre}__minctx{mc}"; t = "" if mc == "6" else "_mc9"
            fa = J(f"ctx_functional_axes{sfx}.json")
            if fa and L in fa["taps"]:
                for a, r in fa["taps"][L].items():
                    put(f"mm_{k}_fz_{AX[a]}{t}", sg(r["z_raw"], 1)); put(f"mm_{k}_fzrc_{AX[a]}{t}", sg(r["z_rank_controlled"], 1))
                    put(f"mm_{k}_auc_{AX[a]}{t}", f"{r['validity_auc']:.2f}")
            v3 = J(f"ctx_coexpr_null_v3{sfx}.json")
            if v3:
                for a, r in v3["summary"].items():
                    put(f"mm_{k}_pm_{AX[a]}{t}", f"{r['p_matched']:.2f}" if r["p_matched"] >= 0.01 else pfmt(r["p_matched"], 300))
                    put(f"mm_{k}_ps_{AX[a]}{t}", f"{r['p_strong_modules_v2']:.2f}" if r["p_strong_modules_v2"] >= 0.01 else pfmt(r["p_strong_modules_v2"], 300))
                pm = [r["p_matched"] for r in v3["summary"].values()]
                put(f"mm_{k}_pm_min{t}", f"{min(pm):.2f}" if min(pm) >= 0.01 else pfmt(min(pm), 300))
                put(f"mm_{k}_pm_nsig{t}", str(sum(p < 0.05 for p in pm)))
            v2 = J(f"ctx_coexpr_null_v2{sfx}.json")
            if v2 and v3:
                rf, rm, rs = [], [], []
                for a, r3_ in v3["axes"].items():
                    r2_ = v2["axes"][a]; rnd = r2_["random_mean"]
                    rf.append(r3_["power"] / rnd); rm.append(r3_["matched_mean"] / rnd); rs.append(r2_["indep_coexpr_mean"] / rnd)
                    put(f"mm_{k}_ratio_func_{AX[a]}{t}", f"{r3_['power'] / rnd:.1f}")
                    put(f"mm_{k}_ratio_matched_{AX[a]}{t}", f"{r3_['matched_mean'] / rnd:.1f}")
                    put(f"mm_{k}_ratio_strong_{AX[a]}{t}", f"{r2_['indep_coexpr_mean'] / rnd:.1f}")
                    put(f"mm_{k}_pr_{AX[a]}{t}", pfmt(r2_["random_p"], len(r2_["null_random"])))
                cend = []
                for a_, r3_ in v3["axes"].items():
                    cur = r3_.get("curve") or []
                    full_ = [c_["power"] for c_ in cur if c_["f"] == 1.0]
                    if full_: cend.append(float(np.mean(full_)) / v2["axes"][a_]["random_mean"])
                rng_ = lambda v: f"{min(v):.1f}–{max(v):.1f}" if f"{min(v):.1f}" != f"{max(v):.1f}" else f"{min(v):.1f}"
                if cend: put(f"mm_{k}_curve_full{t}", rng_(cend))
                put(f"mm_{k}_ratio_func{t}", rng_(rf)); put(f"mm_{k}_ratio_matched{t}", rng_(rm)); put(f"mm_{k}_ratio_strong{t}", rng_(rs))
            ti = J(f"ctx_tightness_null{sfx}.json")
            if ti and L in ti["taps"]:
                for a, r in ti["taps"][L]["axes"].items():
                    put(f"mm_{k}_pt_{AX[a]}{t}", pfmt(r["matched_p"], r["matched_n"]))
            v3e = J(f"ctx_coexpr_null_v3{sfx}__cov-expr.json")
            if v3e:
                for a, r in v3e["summary"].items():
                    put(f"mm_{k}_pmE_{AX[a]}{t}", f"{r['p_matched']:.2f}" if r["p_matched"] >= 0.01 else pfmt(r["p_matched"], 300))
            fae = J(f"ctx_functional_axes{sfx}__cov-expr.json")
            if fae and L in fae["taps"]:
                for a, r in fae["taps"][L].items():
                    put(f"mm_{k}_fzrcE_{AX[a]}{t}", sg(r["z_rank_controlled"], 1))
        dp = J(f"ctx_directional_probe__{pre}__minctx6.json")
        if dp and L in dp["taps"]:
            ps = [(r["p"], r["null_n"]) for r in dp["taps"][L]["axes"].values()]
            for a, r in dp["taps"][L]["axes"].items():
                put(f"mm_{k}_beta_{AX[a]}", sg(r["beta"], 2)); put(f"mm_{k}_beta_p_{AX[a]}", pfmt(r["p"], r["null_n"]))
            put(f"mm_{k}_beta_pmin", pfmt(*min(ps)))
            put(f"mm_{k}_beta_holm", str(len(dp.get("holm_significant", []))))
            put(f"mm_{k}_beta_ntests", str(dp.get("holm_n_tests", "")))
        ct = J(f"ctx_curated_targets__{pre}.json")
        if ct and L in ct["taps"]:
            r = ct["taps"][L]
            put(f"mm_{k}_tf_n", str(r["n_tfs"])); put(f"mm_{k}_tf_static", sg(r["static_excess_mean"]))
            put(f"mm_{k}_tf_static_p", f"{r['static_sign_p']:.3f}" if r["static_sign_p"] >= 0.001 else sci(r["static_sign_p"]))
            put(f"mm_{k}_tf_mod", sg(r["modulation_mean_rho"])); put(f"mm_{k}_tf_mod_p", f"{r['modulation_sign_p']:.2f}")
            if "specific_excess_mean" in r and r.get("specific_n"):
                put(f"mm_{k}_tf_spec", sg(r["specific_excess_mean"])); put(f"mm_{k}_tf_spec_pos", f"{r['specific_n_pos']} of {r['specific_n']}")
                put(f"mm_{k}_tf_spec_p", f"{r['specific_sign_p']:.3f}" if r["specific_sign_p"] >= 0.001 else sci(r["specific_sign_p"]))
    ex = [float(N[f"mm_{k}_excess"].replace("−", "-")) for k in ("m217", "m1b", "scgpt", "state") if f"mm_{k}_excess" in N]
    if ex: put("mm_excess_range", f"{min(ex):.2f} to {max(ex):.2f}")
    pm = [float(N[f"mm_{k}_pm_min"]) for k in ("m217", "m1b", "scgpt", "state") if f"mm_{k}_pm_min" in N and "/" not in N[f"mm_{k}_pm_min"]]
    if pm: put("mm_pm_min_all", f"{min(pm):.2f}")
    def spec_sentence(spec, pos, pp, praw, who):
        if praw is None: return ""
        b3 = min(1.0, 3 * praw)
        if b3 < 0.05:
            return (f"Compared with the targets of the other TFs tested, {who} is closer to its own targets ({spec}; {pos} TFs "
                    f"positive; p = {pp}; {b3:.3f} after Bonferroni over the three curated-target questions), so the closeness "
                    "is partly specific to each TF's own targets.")
        if praw < 0.05:
            return (f"Compared with the targets of the other TFs tested, {who} is slightly closer to its own targets ({spec}; "
                    f"{pos} TFs positive; p = {pp}, not corrected; {b3:.3f} after Bonferroni over the three curated-target "
                    "questions). So we cannot say that the closeness is specific to each TF's own targets.")
        return (f"Compared with the targets of the other TFs tested, {who} is only slightly closer to its own targets "
                f"({spec}; {pos} TFs positive; p = {pp}), so the closeness is not shown to be specific to each TF's own "
                "targets; much of it may reflect where regulators and their usual targets sit in general.")
    ct = J("ctx_curated_targets.json")
    if ct and "specific_sign_p" in ct["taps"].get("L04", {}):
        r = ct["taps"]["L04"]
        pp_ = f"{r['specific_sign_p']:.3f}" if r["specific_sign_p"] >= 0.001 else sci(r["specific_sign_p"])
        put("tf_spec_sentence", spec_sentence(sg(r["specific_excess_mean"]), f"{r['specific_n_pos']} of {r['specific_n']}",
                                              pp_, r["specific_sign_p"], "a TF"))
    cs = J("ctx_curated_targets__ctx_state600.json")
    if cs and "specific_sign_p" in cs["taps"].get("L06", {}):
        r = cs["taps"]["L06"]
        put("state_spec_sentence", spec_sentence(N["mm_state_tf_spec"], N["mm_state_tf_spec_pos"], N["mm_state_tf_spec_p"], r["specific_sign_p"], "a STATE TF"))
        r0 = ct["taps"]["L04"] if ct else None
        both = [x for x in (r["specific_sign_p"], r0 and r0.get("specific_sign_p")) if x is not None]
        if all(min(1.0, 3 * x) >= 0.05 for x in both):
            put("disc_spec_sentence", "Compared with the targets of other TFs, the closeness is not clearly specific to each TF's own "
                f"targets in either model (STATE: {N['mm_state_tf_spec_pos']} TFs positive, p = {N['mm_state_tf_spec_p']} before "
                f"correction; MaxToki-217M: {r0['specific_n_pos']} of {r0['specific_n']} TFs, p = {r0['specific_sign_p']:.3f}), so it "
                "may partly reflect where regulators and their usual targets sit in general.")
        else:
            put("disc_spec_sentence", "Compared with the targets of other TFs, the closeness is specific to each TF's own targets in "
                + ("STATE." if min(1.0, 3 * r["specific_sign_p"]) < 0.05 else "MaxToki-217M."))
    p95 = []
    for k in ("m217", "m1b", "scgpt", "state"):
        v3_ = J(f"ctx_coexpr_null_v3__{dict((a, b) for a, b, _ in MM)[k]}__minctx6.json")
        if v3_: p95 += [r3["matched_p95"] / r3["matched_mean"] for r3 in v3_["axes"].values() if r3.get("matched_mean")]
    if p95: put("mm_pm_p95_over_mean_all", f"{min(p95):.1f} to {max(p95):.1f}")
    tfn = [int(N[f"mm_{k}_tf_n"]) for k, _, _ in MM if f"mm_{k}_tf_n" in N]
    if tfn: put("mm_tf_n_range", f"{min(tfn)} to {max(tfn)}" if min(tfn) != max(tfn) else str(min(tfn)))
    # is the static TF-target closeness learned? the same test on the untrained models (same cells, panel, layer)
    parts = []
    for pre, t, lab in (("ctx_state600rand", 6, "STATE"), ("ctx_scgpt600rand", 4, "scGPT"), ("ctxrand", 4, "MaxToki-217M")):
        c = J(f"ctx_curated_targets__{pre}.json")
        if c and f"L{t:02d}" in c["taps"]:
            r = c["taps"][f"L{t:02d}"]
            put(f"un_{pre}_tf_static", sg(r["static_excess_mean"])); put(f"un_{pre}_tf_n", str(r["n_tfs"]))
            put(f"un_{pre}_tf_static_pos", f"{round(r['static_frac_pos'] * r['n_tfs'])} of {r['n_tfs']}")
            pp = r["static_sign_p"]; put(f"un_{pre}_tf_static_p", f"{pp:.3f}" if pp >= 0.001 else sci(pp))
            parts.append(f"{lab} {sg(r['static_excess_mean'])} (p = {N[f'un_{pre}_tf_static_p']}, {r['n_tfs']} TFs)")
            if r.get("specific_n"):
                put(f"un_{pre}_tf_spec", sg(r["specific_excess_mean"])); put(f"un_{pre}_tf_spec_pos", f"{r['specific_n_pos']} of {r['specific_n']}")
                put(f"un_{pre}_tf_spec_p", f"{r['specific_sign_p']:.3f}" if r["specific_sign_p"] >= 0.001 else sci(r["specific_sign_p"]))
    if all(f"un_{x}_tf_static" in N for x in ("ctx_state600rand", "ctx_scgpt600rand", "ctxrand")):
        put("state_static_untrained_sentence",
            "The untrained STATE, which uses the same fixed protein-embedding gene table as the trained model, shows the "
            f"same closeness ({N['un_ctx_state600rand_tf_static']}, p = {N['un_ctx_state600rand_tf_static_p']})"
            + (f" and is about as specific to each TF's own targets ({N['un_ctx_state600rand_tf_spec']}; "
               f"{N['un_ctx_state600rand_tf_spec_pos']} TFs positive; p = {N['un_ctx_state600rand_tf_spec_p']})"
               if "un_ctx_state600rand_tf_spec" in N else "") +
            f". Untrained scGPT and MaxToki-217M show no significant closeness ({N['un_ctx_scgpt600rand_tf_static']}, p = "
            f"{N['un_ctx_scgpt600rand_tf_static_p']}; {N['un_ctxrand_tf_static']}, with {N.get('un_ctxrand_tf_static_pos', '')} TFs "
            f"positive, p = {N['un_ctxrand_tf_static_p']}; this {N['un_ctxrand_tf_n']}-TF test has little power). So "
            "STATE's closeness does not come from training on expression. It most likely comes from the protein-sequence "
            "embeddings; we did not run STATE with a random gene table, which would test this directly.")
    for tag in ("all600", "all600_common"):
        cm = J(f"ctx_cross_model__{tag}.json")
        if not cm: continue
        t = "" if tag == "all600" else "c"
        if "common_genes" in cm: put("cmc_common_genes", str(cm["common_genes"]))
        if "common_entries" in cm: put("cmc_common_entries", f"{cm['common_entries']:,}")
        for r_ in cm.values():
            if isinstance(r_, dict) and (r_.get("func_z") or {}).get("nuclear_vs_surface"):
                put("cm_n_random", str(r_["func_z"]["nuclear_vs_surface"].get("n_random", 200))); break
        for lab, r in cm.items():
            if not isinstance(r, dict) or "excess" not in r: continue
            key = lab.replace("MaxToki-", "mt").replace("expression-", "expr").replace("-", "").replace("_", "")
            put(f"cm{t}_{key}_excess", sg(r["excess"])); put(f"cm{t}_{key}_genes", str(r.get("n_genes_balanced", r["n_genes"])))
            put(f"cm{t}_{key}_pairs", str(r.get("pairs_scored", ""))); put(f"cm{t}_{key}_dim", str(r["dim"]))
            put(f"cm{t}_{key}_ctx", str(r["n_ctx"]))
            fz = r["func_z"].get("nuclear_vs_surface")
            if fz: put(f"cm{t}_{key}_fz", sg(fz["z"], 1))
    cm = J("ctx_cross_model__all600.json")
    if cm:
        for k, tr, un in (("mt", "MaxToki-217M", "MaxToki-217M-random"), ("sc", "scGPT-L4", "scGPT-random"),
                          ("st", "STATE-L6", "STATE-random")):
            if tr in cm and un in cm and "excess" in cm[un]:
                put(f"cm_randpct_{k}", f"{100 * cm[un]['excess'] / cm[tr]['excess']:.0f}")
        rp_ = [int(N[f"cm_randpct_{k}"]) for k in ("mt", "sc", "st") if f"cm_randpct_{k}" in N]
        if rp_: put("cm_randpct_range", f"{min(rp_)}–{max(rp_)}")
    # steering: 217M vs 1B at matched depth (centred push unit) + the headline configs, each with a FIXED random-split
    # null (30 random splits of the pooled pole genes, same cells/positions/push; p = (1 + #{|split| >= |func|})/31)
    STEER = [("ctx_causal__217m_ctx217m600_centred_splits.json", "m217"),
             ("ctx_causal_mito_cyto_L04__217m_ctx217m600_centred_splits.json", "m217_mito"),
             ("ctx_causal_trans_transport_L04__217m_ctx217m600_centred_splits.json", "m217_trans"),
             ("ctx_causal_nuc_surf_L07__1b_ctx1b_centred_splits.json", "m1b"),
             ("ctx_causal_nuc_surf_L07__1b_ctx1b_centred_nomassive_splits.json", "m1b_nm"),
             ("ctx_causal__217m_ctx217m600_centred_nomassive_splits.json", "m217_nm"),
             ("ctx_causal_mito_cyto_L07__1b_ctx1b_centred_splits.json", "m1b_mito"),
             ("ctx_causal_trans_transport_L07__1b_ctx1b_centred_splits.json", "m1b_trans"),
             ("ctx_causal_splits.json", "nucL4"), ("ctx_causal_nuc_surf_L08_splits.json", "nucL8"),
             ("ctx_causal_mito_cyto_L04_splits.json", "mitoL4"), ("ctx_causal_trans_transport_L04_splits.json", "transL4")]
    for fn, k in STEER:
        c = J(fn)
        if not c: continue
        m = c["signed"]["alpha_0.5"]; sn = c.get("split_null")
        put(f"st_{k}_swing", sg(m["func_swing"])); put(f"st_{k}_rand", sg(m["rand_swing"]))
        put(f"st_{k}_k", f"{m['func_gt_rand']} of {m['n']}")
        put(f"st_{k}_p", sci(m["sign_p"]) if m["sign_p"] < 1e-3 else f"{m['sign_p']:.2g}")
        put(f"st_{k}_swings", ", ".join(sg(c["signed"][f"alpha_{a}"]["func_swing"]) for a in c["alphas_xResidNorm"]))
        if "resid_norm_raw" in c:
            put(f"st_{k}_norm_raw", f"{c['resid_norm_raw']:.0f}"); put(f"st_{k}_norm_c", f"{c['resid_norm_centred']:.0f}")
        put(f"st_{k}_massive_share", f"{100 * c['natural_movement_along_axis'].get('axis_share_in_massive_dims', 0):.0f}")
        if sn:
            put(f"st_{k}_split_abs", f"{sn['split_abs_mean']:.3f}"); put(f"st_{k}_split_nge", f"{sn['n_splits_abs_ge_func']} of {sn['n_splits']}")
            put(f"st_{k}_split_p", pfmt(sn["p_abs"], sn["n_splits"])); put(f"st_{k}_split_z", sg(sn["z_abs"], 1))
            put(f"st_{k}_split_nge_n", str(sn["n_splits_abs_ge_func"]))
            put(f"st_{k}_split_beat", f"{sn['n_splits'] - sn['n_splits_abs_ge_func']} of {sn['n_splits']}")
            put("st_n_splits", str(sn["n_splits"])); put("st_split_pmin", f"{1 / (sn['n_splits'] + 1):.3f}")
        if "resid_norm_raw" in c and c.get("resid_norm_centred"):
            put(f"st_{k}_rawratio", f"{c['resid_norm_raw'] / c['resid_norm_centred']:.1f}")
    # prediction link, paired 217M vs 1B
    for fn, k in (("ctx_prediction_link_paired__217m_ctx217m600_L04.json", "m217"),
                  ("ctx_prediction_link_paired__1b_ctx1b_L07.json", "m1b"),
                  ("ctx_prediction_link_paired__217m_ctx_maxtoki_L04.json", "m217k")):
        p = J(fn)
        if not p: continue
        _eq = max(N.get("_plp_equiv_raw", 0.0), p["equivalence_max_abs_diff"]); N["_plp_equiv_raw"] = _eq
        put("plp_equiv", f"{_eq:.1e}".replace("e-0", " × 10^−").replace("e-", " × 10^−") + "^")
        for arm, a in (("random_donor", "r"), ("same_celltype_donor", "s")):
            r = p["arms"].get(arm)
            if not r or "mean_benefit" not in r: continue
            put(f"plp_{k}_{a}_benefit", f"{r['mean_benefit']:+.2f}".replace("-", "−")); put(f"plp_{k}_{a}_n", str(r["n_genes"]))
            put(f"plp_{k}_{a}_rho", sg(r["partial_rho"])); put(f"plp_{k}_{a}_ci", f"{sg(r['partial_ci'][0])} to {sg(r['partial_ci'][1])}")


try:
    _sec11()
except Exception as e:
    print('section multi-model failed:', repr(e)[:200])
_panel_stats_save(); _panel_stats_fallback()


# ---------------- second-review additions (3 Oct 2026) ----------------
def _sec12():
    global N
    AXK = {"nuclear_vs_surface": "nuc", "mito_vs_cytoskeleton": "mito", "transcription_vs_transport": "trans"}
    # tightness null: how often a functional pole was less tight than random genes (a random pole was then used)
    fb = {}
    for lab, fn, taps_ in (("main", "ctx_tightness_null.json", ("L04", "L08")),) + tuple(
            (k, f"ctx_tightness_null__{pre}__minctx{mc}.json", (f"L{t:02d}",)) for k, pre, t in MM for mc in ("6", "9")):
        d = J(fn)
        if not d: continue
        for L in taps_:
            for a, r in d["taps"].get(L, {}).get("axes", {}).items():
                fb[(lab, L, a)] = max(r.get("matched_poles_below_random", [0, 0]))
    models = [v for (lab, L, a), v in fb.items() if lab in ("main", "m217", "m1b", "scgpt", "state")]
    if models: put("tight_fb_max_models", str(max(models)))
    for k in ("exprL", "exprP", "exprC"):
        v = {a: fb.get((k, "L00", a), 0) for a in AXK}
        put(f"tight_fb_{k}_max", str(max(v.values())))
        for a, x in v.items(): put(f"tight_fb_{k}_{AXK[a]}", str(x))
    put("tight_fb_main_L8_trans", str(fb.get(("main", "L08", "transcription_vs_transport"), 0)))
    AXW = {"nuc": "nuclear/surface", "mito": "mitochondrion/cytoskeleton", "trans": "transcription/transport"}
    bits = [f"This happened in at most {N['tight_fb_max_models']} of 300 null axes for any axis in the four trained models"]
    for k, nm in (("exprL", "landmark"), ("exprP", "principal-component"), ("exprC", "centroid")):
        mx = int(N[f"tight_fb_{k}_max"])
        if mx == 0: continue
        big = [a for a in ("nuc", "mito", "trans") if int(N[f"tight_fb_{k}_{a}"]) > 50]
        if big:
            bits.append(f"In the {nm} expression-only version a random pole was used in " +
                        " and ".join(f"{N[f'tight_fb_{k}_{a}']} of 300 null axes for the {AXW[a]} axis" for a in big) +
                        ", so its tightness p values for " + ("that axis are upper limits" if len(big) == 1 else "those axes are upper limits"))
        else:
            bits.append(f"In the {nm} expression-only version it happened in at most {mx} of 300")
    put("tight_fb_sentence", ". ".join(bits) + ".")
    # rank control: EXCESS measured against the rank-matched gene instead of a random gene
    pc, pol = J("ctx_position_confound.json"), J("ctx_polysemy.json")
    if pc and pol:
        for L, r in pc["taps"].items():
            if "excess_vs_rank_matched" in r and L in pol["taps"]:
                l = int(L[1:]); put(f"rank_vs_matched_L{l}", sg(r["excess_vs_rank_matched"]))
                put(f"rank_vs_matched_pct_L{l}", f"{100 * r['excess_vs_rank_matched'] / pol['taps'][L]['excess']:.0f}")
    # depth pattern of EXCESS in the 600-cell representations (Table 3 rule: pairs with >= 200 genes)
    for k, pre in (("m217", "ctx217m600"), ("m1b", "ctx1b"), ("scgpt", "ctx_scgpt600"), ("state", "ctx_state600")):
        p_ = J(f"ctx_polysemy__{pre}.json")
        if p_:
            for L, r in p_["taps"].items(): put(f"dp_{k}_L{int(L[1:])}", sg(r["excess"]))
    # co-expression curve at fixed pole size (main extraction, nuclear/surface): null axes >= functional, per share
    v3 = J("ctx_coexpr_null_v3.json")
    if v3:
        a = v3["axes"]["nuclear_vs_surface"]; by = {}
        for e in a.get("curve", []): by.setdefault(e["f"], []).append(e)
        for f_, es in by.items():
            tag = f"{int(round(100 * f_))}"
            put(f"curve_nuc_f{tag}_nge", f"{sum(e['power'] >= a['power'] for e in es)} of {len(es)}")
            put(f"curve_nuc_f{tag}_coh", f"{np.mean([e['cohA'] for e in es]):.3f}")
            put(f"curve_nuc_f{tag}_cohB", f"{np.mean([e['cohB'] for e in es]):.3f}")
        put("curve_nuc_ndraws", str(len(next(iter(by.values())))) if by else "")
        put("curve_nuc_power", f"{a['power']:.2f}")
        put("coh_nuc_A3", f"{a['coh_A']:.3f}"); put("coh_nuc_B3", f"{a['coh_B']:.3f}")
    # steering: natural range of the projection; Bonferroni for the smallest push
    c = J("ctx_causal.json")
    if c:
        nat = c["natural_movement_along_axis"]
        put("cz_nucL4_nat_p95", f"{nat['p95_range_across_contexts']:.2f}")
        a0 = c["signed"]["alpha_0.25"]; put("cz_nucL4_a0_p_bonf", f"{min(1, 4 * a0['sign_p']):.3f}")
        put("cz_nucL4_push_hi_over_r", f"{c['push_size_raw'][2] / nat['median_range_across_contexts']:.0f}")
        put("cz_nucL4_push_lo_over_r", f"{c['push_size_raw'][0] / nat['median_range_across_contexts']:.1f}")
    # directional congruence: null spread per representation; tightness Holm over 3 axes; curated details
    for k, pre, tap in MM:
        L = f"L{tap:02d}"
        dp = J(f"ctx_directional_probe__{pre}__minctx6.json")
        if dp and L in dp["taps"]:
            sds = []
            for a, r in dp["taps"][L]["axes"].items():
                put(f"mm_{k}_beta_nullsd_{AXK.get(a, a)}", f"{r['null_sd']:.2f}"); put(f"mm_{k}_beta_nullmean_{AXK.get(a, a)}", sg(r["null_mean"], 2))
                sds.append(r["null_sd"])
            put(f"mm_{k}_beta_nullsd_range", f"{min(sds):.2f}–{max(sds):.2f}" if f"{min(sds):.2f}" != f"{max(sds):.2f}" else f"{min(sds):.2f}")
        ti = J(f"ctx_tightness_null__{pre}__minctx6.json")
        if ti and L in ti["taps"]:
            ps = sorted((r["matched_p"], a) for a, r in ti["taps"][L]["axes"].items() if a in AXK)
            run = 0.0
            for i, (p_, a) in enumerate(ps):
                run = max(run, min(1.0, (len(ps) - i) * p_)); put(f"mm_{k}_pt_holm_{AXK[a]}", f"{run:.2f}")
        ct = J(f"ctx_curated_targets__{pre}.json")
        if ct:
            for L_, r in ct["taps"].items():
                sfx = "" if L_ == L else f"_L{int(L_[1:])}"
                put(f"mm_{k}_tf_static_pos{sfx}", f"{round(r['static_frac_pos'] * r['n_tfs'])} of {r['n_tfs']}")
                put(f"mm_{k}_tf_mod_pos{sfx}", f"{round(r['modulation_frac_pos'] * r['n_modul'])} of {r['n_modul']}")
                put(f"mm_{k}_tf_mod_p_bonf3{sfx}", f"{min(1, 3 * r['modulation_sign_p']):.2f}")
                if r.get("specific_n"):
                    put(f"mm_{k}_tf_spec{sfx}", sg(r["specific_excess_mean"])); put(f"mm_{k}_tf_spec_pos{sfx}", f"{r['specific_n_pos']} of {r['specific_n']}")
                    pp = r["specific_sign_p"]
                    put(f"mm_{k}_tf_spec_p{sfx}", f"{pp:.3f}" if pp >= 0.001 else sci(pp))
                    put(f"mm_{k}_tf_spec_p_bonf3{sfx}", f"{min(1, 3 * pp):.3f}")
    for pre, t in (("ctx_state600rand", 6), ("ctx_scgpt600rand", 4), ("ctxrand", 4)):
        c_ = J(f"ctx_curated_targets__{pre}.json")
        if c_ and f"L{t:02d}" in c_["taps"]:
            r = c_["taps"][f"L{t:02d}"]
            put(f"un_{pre}_tf_static_pos", f"{round(r['static_frac_pos'] * r['n_tfs'])} of {r['n_tfs']}")
            if r.get("specific_n"):
                put(f"un_{pre}_tf_spec", sg(r["specific_excess_mean"])); put(f"un_{pre}_tf_spec_pos", f"{r['specific_n_pos']} of {r['specific_n']}")
                pp = r["specific_sign_p"]; put(f"un_{pre}_tf_spec_p", f"{pp:.3f}" if pp >= 0.001 else sci(pp))
        fa = J(f"ctx_functional_axes__{pre}__minctx6.json")
        if fa and f"L{t:02d}" in fa["taps"]:
            aucs = []
            for a, r in fa["taps"][f"L{t:02d}"].items():
                if a in AXK:
                    put(f"un_{pre}_fzrc_{AXK[a]}", sg(r["z_rank_controlled"], 1)); aucs.append(r["validity_auc"])
            put(f"un_{pre}_auc_range", f"{min(aucs):.2f}–{max(aucs):.2f}")
    # what else differs between cell types: dataset, donor, assay (the selected cells)
    comp = J("ctx_context_composition.json")
    if comp:
        m = comp["cells_1000"]; one = m["cell_types_one_dataset_by_dataset"]
        put("comp_one_ds_n", str(len(one))); put("comp_one_ds_lung", str(sum(v == "lung" for v in one.values())))
        per = m["per_cell_type"]
        put("comp_kidney_donors", str(per["kidney epithelial cell"]["n_donors"]))
        r = per["pulmonary alveolar type 1 cell"]; put("comp_at1_top", f"{r['top_donor_n']} of {r['n']}")
        put("comp_min_10x", f"{100 * m['min_share_10x_3v3']:.0f}")
    # anisotropy: genes behind the self-similarity mean
    an = J("ctx_anisotropy.json")
    if an and "selfsim_n_genes" in an["layers"].get("L04", {}):
        put("selfsim_n_genes_L4", str(an["layers"]["L04"]["selfsim_n_genes"]))
    cmc = J("ctx_cross_model__all600_common.json")
    if cmc and "MaxToki-217M" in cmc: put("cmc_pairs", str(cmc["MaxToki-217M"].get("pairs_scored", "")))
    ec = J("ctx_expression_cellcounts.json")
    if ec:
        e = ec["expressing_cells_per_entry"]
        put("expr_cells_median", f"{e['median']:.0f}"); put("expr_cells_p10", f"{e['p10']:.0f}"); put("expr_cells_p90", f"{e['p90']:.0f}")
    sw = J("ctx_switcher_test.json")
    if sw and "n_switchers_in_panel" in sw:
        put("switch_n_in_panel", str(sw["n_switchers_in_panel"]))


try:
    _sec12()
except Exception as e:
    print('section second-review failed:', repr(e)[:200])


# ---------------- derived / abstract-level values
def _derived():
    global N
    J_ = J
    pol = J_("ctx_polysemy.json"); cm = J_("ctx_cross_model.json"); v3 = J_("ctx_coexpr_null_v3.json")
    v2 = J_("ctx_coexpr_null_v2.json"); ind = J_("ctx_independent.json"); pl = J_("ctx_prediction_link.json")
    fa = J_("ctx_functional_axes.json"); ct = J_("ctx_curated_targets.json")
    if pol:
        put("excess_L4_r2", f"{pol['taps']['L04']['excess']:.2f}"); put("diff_L4_r3", f"{pol['taps']['L04']['diff']:.3f}")
        put("excess_ci_halfwidth_max", f"{max((r['ci'][1]-r['ci'][0])/2 for r in pol['taps'].values()):.3f}")
    if ind:
        put("setty_excess_L4_r2", f"{ind['taps']['L04']['excess']:.2f}")
        for L, r in ind["taps"].items():
            m = r["func_z"].get("mito_vs_cytoskeleton")
            if m: put(f"setty_fz_mito_L{int(L[1:])}", sg(m["z"], 1))
    if cm:
        put("cm_rand_excess_r2", f"{cm['MaxToki-217M-random']['excess']:.2f}")
        mx = [cm[k]["func_z"]["nuclear_vs_surface"]["z"] for k in ("MaxToki-217M", "MaxToki-217M-L2", "MaxToki-1B", "MaxToki-1B-L7")]
        put("cm_maxtoki_fz_range", f"{min(mx):.1f}–{max(mx):.1f}")
        put("cm_m217_excess_r2", f"{cm['MaxToki-217M']['excess']:.2f}")
        put("cm_rand_excess_pct", f"{100 * cm['MaxToki-217M-random']['excess'] / cm['MaxToki-217M']['excess']:.0f}")
        put("cm_rand_fz_pct", f"{100 * cm['MaxToki-217M-random']['func_z']['nuclear_vs_surface']['z'] / cm['MaxToki-217M']['func_z']['nuclear_vs_surface']['z']:.0f}")
        put("cm_scgpt_fz_r1", f"{cm['scGPT']['func_z']['nuclear_vs_surface']['z']:.1f}")
        put("cm_state_fz_r1", f"{cm['STATE-SE']['func_z']['nuclear_vs_surface']['z']:.1f}")
    if v3:
        pm = [v["p_matched"] for v in v3["summary"].values()]; psv = [v["p_strong_modules_v2"] for v in v3["summary"].values()]
        put("pm_range", f"{min(pm):.2f}–{max(pm):.2f}"); put("ps_range", f"{min(psv):.2f}–{max(psv):.2f}")
    if pl:
        put("pl_rho_r2", f"{pl['partial_rho']:.2f}")
        put("pl_fold_round", f"{float(f'{math.exp(pl[chr(109)+chr(101)+chr(97)+chr(110)+chr(95)+chr(99)+chr(111)+chr(110)+chr(116)+chr(101)+chr(120)+chr(116)+chr(95)+chr(98)+chr(101)+chr(110)+chr(101)+chr(102)+chr(105)+chr(116)]):.2g}'):,.0f}")
    if fa:
        put("fz_raw_nuc_L4_r0", f"{fa['taps']['L04']['nuclear_vs_surface']['z_raw']:.0f}")
    if ct:
        for L in ("L04", "L08"):
            r = ct["taps"][L]; l = int(L[1:])
            put(f"tf_static_p_L{l}", f"{r['static_sign_p']:.3f}" if r["static_sign_p"] >= 0.001 else f"{r['static_sign_p']:.1e}")
            put(f"tf_mod_p_L{l}", f"{r['modulation_sign_p']:.2f}"); put(f"tf_static_L{l}", sg(r["static_excess_mean"]))
            put(f"tf_mod_L{l}", sg(r["modulation_mean_rho"]))
        put("tf_static_p", N[f"tf_static_p_L4"]); put("tf_static_p_bonf4", f"{min(1, 4 * ct['taps']['L04']['static_sign_p']):.3f}")
        put("tf_mod_p_bonf4", f"{min(1, 4 * ct['taps']['L04']['modulation_sign_p']):.2f}")
        for L in ("L04", "L08"):
            r = ct["taps"][L]; l = int(L[1:])
            if "specific_excess_mean" in r and r.get("specific_n"):
                put(f"tf_spec_L{l}", sg(r["specific_excess_mean"])); put(f"tf_spec_pos_L{l}", f"{r['specific_n_pos']} of {r['specific_n']}")
                put(f"tf_spec_p_L{l}", f"{r['specific_sign_p']:.3f}" if r["specific_sign_p"] >= 0.001 else sci(r["specific_sign_p"]))
    for fn, tag in [("ctx_causal.json", "nucL4"), ("ctx_causal_nuc_surf_L08.json", "nucL8"),
                    ("ctx_causal_mito_cyto_L04.json", "mitoL4"), ("ctx_causal_trans_transport_L04.json", "transL4")]:
        c = J_(fn)
        if not c: continue
        m = c["signed"]["alpha_0.5"]; pb = min(1.0, 4 * m["sign_p"])
        put(f"cz_{tag}_p_bonf", f"{pb:.3f}" if pb >= 0.001 else sci(pb))
        put(f"cz_{tag}_push_mid", f"{c['push_size_raw'][1]:.1f}")
        a1 = c["signed"]["alpha_1.0"]; a0 = c["signed"]["alpha_0.25"]
        put(f"cz_{tag}_a1_swing", sg(a1["func_swing"])); put(f"cz_{tag}_a1_k", f"{a1['func_gt_rand']} of {a1['n']}")
        put(f"cz_{tag}_a1_p", f"{a1['sign_p']:.3f}"); put(f"cz_{tag}_a1_p_bonf", f"{min(1, 4 * a1['sign_p']):.2f}")
        put(f"cz_{tag}_a0_k", f"{a0['func_gt_rand']} of {a0['n']}"); put(f"cz_{tag}_a0_p", f"{a0['sign_p']:.2g}")
        put(f"cz_{tag}_a0_swing", sg(a0["func_swing"]))
        put(f"cz_{tag}_push_lo", f"{c['push_size_raw'][0]:.1f}"); put(f"cz_{tag}_push_hi", f"{c['push_size_raw'][2]:.1f}")
        nat = c["natural_movement_along_axis"]
        put(f"cz_{tag}_push_lo_over", f"{c['push_size_raw'][0] / nat['median_range_across_contexts']:.1f}")
        put(f"cz_{tag}_push_hi_over", f"{c['push_size_raw'][2] / nat['median_range_across_contexts']:.0f}")


try:
    _derived()
except Exception as e:
    print("section derived failed:", repr(e)[:200])

# steering sentences for the abstract and cover letter (fixed text; the numbers are in Results)
put("ABS_STEER", "Adding the nuclear/surface direction to half of a cell's genes shifts MaxToki's predictions for the other half more than random directions from the same genes; as the pushes exceed natural changes, this shows the model can read the direction, not that it uses it.")
put("COVER_STEER", "Steering along the nuclear/surface direction changes MaxToki's predictions for untouched genes in both sizes, more than directions built from random groups of the same genes; MaxToki-1B also reads the other two functional directions tested, while MaxToki-217M reads one of them in reverse. The pushes are larger than natural changes, so this shows that the models can read these directions, not that they use the natural changes.")
N.pop("_plp_equiv_raw", None)
json.dump(N, open(os.path.join(HERE, "numbers.json"), "w"), indent=1, ensure_ascii=False, default=str)
print(f"{len(N)} numbers -> build/numbers.json")

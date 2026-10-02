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
        put(f"cz_{tag}_a0_k", f"{a0['func_gt_rand']} of {a0['n']}"); put(f"cz_{tag}_a0_p", f"{a0['sign_p']:.4f}")
        put(f"cz_{tag}_push_lo", f"{c['push_size_raw'][0]:.1f}"); put(f"cz_{tag}_push_hi", f"{c['push_size_raw'][2]:.1f}")
        nat = c["natural_movement_along_axis"]
        put(f"cz_{tag}_push_lo_over", f"{c['push_size_raw'][0] / nat['median_range_across_contexts']:.1f}")
        put(f"cz_{tag}_push_hi_over", f"{c['push_size_raw'][2] / nat['median_range_across_contexts']:.0f}")


try:
    _derived()
except Exception as e:
    print("section derived failed:", repr(e)[:200])

json.dump(N, open(os.path.join(HERE, "numbers.json"), "w"), indent=1, ensure_ascii=False, default=str)
print(f"{len(N)} numbers -> build/numbers.json")

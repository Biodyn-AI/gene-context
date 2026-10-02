#!/bin/bash
# Every multi-model command run for the gene-context paper on 2-3 Oct 2026, in run order (generated from the
# job queue; each block is one job). Run from src/ with PY set to a Python that has the packages in requirements.txt
# (scGPT and STATE extraction need their own environments, see environment.md).
set -e
PY=${PY:-python}
# --- 098a_causal_217m600_L4_centred.sh
MODEL=217m PREFIX=ctx217m600 AXIS=nuc_surf SITE=3 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 098b_causal_1b_L7_centred.sh
MODEL=1b PREFIX=ctx1b AXIS=nuc_surf SITE=6 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 098c_predlink_paired_217m600.sh
PAIRED=1 MODEL=217m PREFIX=ctx217m600 TAP=4 $PY -u ctx_prediction_link.py
# --- 098d_predlink_paired_1b.sh
PAIRED=1 MODEL=1b PREFIX=ctx1b TAP=7 $PY -u ctx_prediction_link.py
# --- 098e_predlink_paired_217m_headline.sh
PAIRED=1 MODEL=217m PREFIX=ctx_maxtoki TAP=4 $PY -u ctx_prediction_link.py
# --- 098f_causal_1b_L7_nomassive.sh
MODEL=1b PREFIX=ctx1b AXIS=nuc_surf SITE=6 ALPHA_UNIT=centred ZERO_MASSIVE=1 $PY -u ctx_causal.py
# --- 098g_causal_217m600_L4_nomassive.sh
MODEL=217m PREFIX=ctx217m600 AXIS=nuc_surf SITE=3 ALPHA_UNIT=centred ZERO_MASSIVE=1 $PY -u ctx_causal.py
# --- 098h_causal_217m600_mito.sh
MODEL=217m PREFIX=ctx217m600 AXIS=mito_cyto SITE=3 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 098i_causal_217m600_trans.sh
MODEL=217m PREFIX=ctx217m600 AXIS=trans_transport SITE=3 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 098j_causal_1b_mito.sh
MODEL=1b PREFIX=ctx1b AXIS=mito_cyto SITE=6 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 098k_causal_1b_trans.sh
MODEL=1b PREFIX=ctx1b AXIS=trans_transport SITE=6 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 099a_scgpt600_extract.sh
$PY -u ctx_extract_scgpt_std.py
# --- 099b_scgpt600rand_extract.sh
RANDOM_INIT=1 $PY -u ctx_extract_scgpt_std.py
# --- 100_ctx217m600_polysemy.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 101_ctx217m600_anisotropy.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_anisotropy.py
# --- 102_ctx217m600_position.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_position_confound.py
# --- 103_ctx217m600_loadings.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_celltype_loadings.py
# --- 104_ctx217m600_mc9_funcaxes.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_functional_axes.py
# --- 105_ctx217m600_mc9_v2.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v2.py
# --- 106_ctx217m600_mc9_v3.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v3.py
# --- 107_ctx217m600_mc9_tightness.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_tightness_null.py
# --- 108_ctx217m600_mc6_funcaxes.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 109_ctx217m600_mc6_v2.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v2.py
# --- 110_ctx217m600_mc6_v3.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v3.py
# --- 111_ctx217m600_mc6_tightness.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_tightness_null.py
# --- 112_ctx217m600_curated.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 113_ctxexpr_polysemy.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 114_ctxexpr_anisotropy.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_anisotropy.py
# --- 115_ctxexpr_position.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_position_confound.py
# --- 116_ctxexpr_loadings.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_celltype_loadings.py
# --- 117_ctxexpr_mc9_funcaxes.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_functional_axes.py
# --- 118_ctxexpr_mc9_v2.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v2.py
# --- 119_ctxexpr_mc9_v3.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v3.py
# --- 120_ctxexpr_mc9_tightness.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_tightness_null.py
# --- 121_ctxexpr_mc6_funcaxes.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 122_ctxexpr_mc6_v2.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v2.py
# --- 123_ctxexpr_mc6_v3.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v3.py
# --- 124_ctxexpr_mc6_tightness.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_tightness_null.py
# --- 125_ctxexpr_curated.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 126_ctxexprpc_polysemy.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 127_ctxexprpc_anisotropy.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_anisotropy.py
# --- 128_ctxexprpc_position.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_position_confound.py
# --- 129_ctxexprpc_loadings.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_celltype_loadings.py
# --- 130_ctxexprpc_mc9_funcaxes.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_functional_axes.py
# --- 131_ctxexprpc_mc9_v2.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v2.py
# --- 132_ctxexprpc_mc9_v3.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v3.py
# --- 133_ctxexprpc_mc9_tightness.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_tightness_null.py
# --- 134_ctxexprpc_mc6_funcaxes.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 135_ctxexprpc_mc6_v2.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v2.py
# --- 136_ctxexprpc_mc6_v3.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v3.py
# --- 137_ctxexprpc_mc6_tightness.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_tightness_null.py
# --- 138_ctxexprpc_curated.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 139_ctx1b_polysemy.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 140_ctx1b_anisotropy.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 $PY -u ctx_anisotropy.py
# --- 141_ctx1b_position.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 $PY -u ctx_position_confound.py
# --- 142_ctx1b_loadings.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 $PY -u ctx_celltype_loadings.py
# --- 143_ctx1b_mc9_funcaxes.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_functional_axes.py
# --- 144_ctx1b_mc9_v2.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v2.py
# --- 145_ctx1b_mc9_v3.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v3.py
# --- 146_ctx1b_mc9_tightness.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_tightness_null.py
# --- 147_ctx1b_mc6_funcaxes.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 148_ctx1b_mc6_v2.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v2.py
# --- 149_ctx1b_mc6_v3.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v3.py
# --- 150_ctx1b_mc6_tightness.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_tightness_null.py
# --- 151_ctx1b_curated.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 1600_ctx_scgpt600_polysemy.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 1601_ctx_scgpt600_anisotropy.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_anisotropy.py
# --- 1602_ctx_scgpt600_position_confound.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_position_confound.py
# --- 1603_ctx_scgpt600_celltype_loadings.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_celltype_loadings.py
# --- 1604_ctx_scgpt600_mc9_functional_axes.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_functional_axes.py
# --- 1605_ctx_scgpt600_mc9_coexpr_null_v2.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v2.py
# --- 1606_ctx_scgpt600_mc9_coexpr_null_v3.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v3.py
# --- 1607_ctx_scgpt600_mc9_tightness_null.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_tightness_null.py
# --- 1608_ctx_scgpt600_mc6_functional_axes.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 1609_ctx_scgpt600_mc6_coexpr_null_v2.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v2.py
# --- 1610_ctx_scgpt600_mc6_coexpr_null_v3.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v3.py
# --- 1611_ctx_scgpt600_mc6_tightness_null.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_tightness_null.py
# --- 1612_ctx_scgpt600_curated.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 1613_ctx_scgpt600_mc6_directional.sh
env PREFIX=ctx_scgpt600 TAPS=4,2,9 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_directional_probe.py
# --- 176_ctx_scgpt600rand_polysemy.sh
env PREFIX=ctx_scgpt600rand TAPS=4,2 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 177_ctx_scgpt600rand_mc6_funcaxes.sh
env PREFIX=ctx_scgpt600rand TAPS=4 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 200_ctx217m600_mc6_directional.sh
env PREFIX=ctx217m600 TAPS=4,2 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_directional_probe.py
# --- 201_ctxexpr_mc6_directional.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_directional_probe.py
# --- 202_ctxexprpc_mc6_directional.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_directional_probe.py
# --- 203_ctx1b_mc6_directional.sh
env PREFIX=ctx1b TAPS=7,4 TAP=7 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_directional_probe.py
# --- 290a_state600_extract.sh
MATCH=ctx217m600 CELLSCTX=600 MAXGENES=5000 TAPS=3,6,12 BATCH=4 OUTPREFIX=ctx_state600 $PY -u ctx_extract_state_std.py
# --- 290b_state600rand_extract.sh
RANDOM_INIT=1 MATCH=ctx217m600 CELLSCTX=600 MAXGENES=5000 TAPS=6 BATCH=4 OUTPREFIX=ctx_state600rand $PY -u ctx_extract_state_std.py
# --- 314_ctx_state600_polysemy.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 315_ctx_state600_anisotropy.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 $PY -u ctx_anisotropy.py
# --- 316_ctx_state600_position_confound.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 $PY -u ctx_position_confound.py
# --- 317_ctx_state600_celltype_loadings.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 $PY -u ctx_celltype_loadings.py
# --- 318_ctx_state600_mc9_functional_axes.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_functional_axes.py
# --- 319_ctx_state600_mc9_coexpr_null_v2.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v2.py
# --- 320_ctx_state600_mc9_coexpr_null_v3.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v3.py
# --- 321_ctx_state600_mc9_tightness_null.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_tightness_null.py
# --- 322_ctx_state600_mc6_functional_axes.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 323_ctx_state600_mc6_coexpr_null_v2.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v2.py
# --- 324_ctx_state600_mc6_coexpr_null_v3.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v3.py
# --- 325_ctx_state600_mc6_tightness_null.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_tightness_null.py
# --- 326_ctx_state600_curated.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 327_ctx_state600_mc6_directional.sh
env PREFIX=ctx_state600 TAPS=6,3,12 TAP=6 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_directional_probe.py
# --- 330_ctx_state600rand_polysemy.sh
env PREFIX=ctx_state600rand TAPS=6 TAP=6 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 331_ctx_state600rand_mc6_funcaxes.sh
env PREFIX=ctx_state600rand TAPS=6 TAP=6 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 400_crossmodel_all600.sh
MODELS=MaxToki-217M-L2:ctx217m600:2,MaxToki-217M:ctx217m600:4,MaxToki-1B-L4:ctx1b:4,MaxToki-1B-L7:ctx1b:7,scGPT-L2:ctx_scgpt600:2,scGPT-L4:ctx_scgpt600:4,STATE-L3:ctx_state600:3,STATE-L6:ctx_state600:6,MaxToki-217M-random:ctxrand:4,scGPT-random:ctx_scgpt600rand:4,STATE-random:ctx_state600rand:6,expression-landmarks:ctxexpr:0,expression-PCs:ctxexprpc:0 OUTTAG=all600 $PY -u ctx_cross_model.py
# --- 401_crossmodel_all600_common.sh
MODELS=MaxToki-217M-L2:ctx217m600:2,MaxToki-217M:ctx217m600:4,MaxToki-1B-L4:ctx1b:4,MaxToki-1B-L7:ctx1b:7,scGPT-L2:ctx_scgpt600:2,scGPT-L4:ctx_scgpt600:4,STATE-L3:ctx_state600:3,STATE-L6:ctx_state600:6,MaxToki-217M-random:ctxrand:4,scGPT-random:ctx_scgpt600rand:4,STATE-random:ctx_state600rand:6,expression-landmarks:ctxexpr:0,expression-PCs:ctxexprpc:0 OUTTAG=all600_common COMMON=1 $PY -u ctx_cross_model.py
# --- 410_ctx217m600_covexpr_funcaxes.sh
env PREFIX=ctx217m600 TAPS=4 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_functional_axes.py
# --- 411_ctx217m600_covexpr_v2.sh
env PREFIX=ctx217m600 TAPS=4 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v2.py
# --- 412_ctx217m600_covexpr_v3.sh
env PREFIX=ctx217m600 TAPS=4 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v3.py
# --- 413_ctx1b_covexpr_funcaxes.sh
env PREFIX=ctx1b TAPS=7 TAP=7 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_functional_axes.py
# --- 414_ctx1b_covexpr_v2.sh
env PREFIX=ctx1b TAPS=7 TAP=7 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v2.py
# --- 415_ctx1b_covexpr_v3.sh
env PREFIX=ctx1b TAPS=7 TAP=7 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v3.py
# --- 416_ctx_scgpt600_covexpr_funcaxes.sh
env PREFIX=ctx_scgpt600 TAPS=4 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_functional_axes.py
# --- 417_ctx_scgpt600_covexpr_v2.sh
env PREFIX=ctx_scgpt600 TAPS=4 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v2.py
# --- 418_ctx_scgpt600_covexpr_v3.sh
env PREFIX=ctx_scgpt600 TAPS=4 TAP=4 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v3.py
# --- 419_ctx_state600_covexpr_funcaxes.sh
env PREFIX=ctx_state600 TAPS=6 TAP=6 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_functional_axes.py
# --- 420_ctx_state600_covexpr_v2.sh
env PREFIX=ctx_state600 TAPS=6 TAP=6 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v2.py
# --- 421_ctx_state600_covexpr_v3.sh
env PREFIX=ctx_state600 TAPS=6 TAP=6 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v3.py
# --- 422_ctxexpr_covexpr_funcaxes.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_functional_axes.py
# --- 423_ctxexpr_covexpr_v2.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v2.py
# --- 424_ctxexpr_covexpr_v3.sh
env PREFIX=ctxexpr TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v3.py
# --- 425_ctx_state600rand_curated.sh
env PREFIX=ctx_state600rand TAPS=6 TAP=6 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 426_ctx_scgpt600rand_curated.sh
env PREFIX=ctx_scgpt600rand TAPS=4 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 427_ctxrand_curated.sh
env PREFIX=ctxrand TAPS=4 TAP=4 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 428a_ctxexprcen_extract.sh
MODE=centroid $PY -u ctx_extract_expression.py
# --- 428b_ctxexprcen_polysemy.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 428c_ctxexprcen_position.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_position_confound.py
# --- 428d_ctxexprcen_mc6_funcaxes.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 428e_ctxexprcen_mc6_v2.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v2.py
# --- 428f_ctxexprcen_mc6_v3.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v3.py
# --- 428g_ctxexprcen_mc6_tightness.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_tightness_null.py
# --- 428h_ctxexprcen_mc6_directional.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_directional_probe.py
# --- 428i_ctxexprcen_curated.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 428j_ctxexprcen_covexpr.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_functional_axes.py && env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v2.py && env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v3.py
# --- 428k_ctxexprcen_mc9.sh
env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_functional_axes.py && env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v2.py && env PREFIX=ctxexprcen TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v3.py
# --- 428l_crossmodel_all600.sh
N_RANDOM=1000 DISJOINT=1 MODELS=MaxToki-217M-L2:ctx217m600:2,MaxToki-217M:ctx217m600:4,MaxToki-1B-L4:ctx1b:4,MaxToki-1B-L7:ctx1b:7,scGPT-L2:ctx_scgpt600:2,scGPT-L4:ctx_scgpt600:4,STATE-L3:ctx_state600:3,STATE-L6:ctx_state600:6,MaxToki-217M-random:ctxrand:4,scGPT-random:ctx_scgpt600rand:4,STATE-random:ctx_state600rand:6,expression-landmarks:ctxexpr:0,expression-PCs:ctxexprpc:0,expression-centroids:ctxexprcen:0 OUTTAG=all600 $PY -u ctx_cross_model.py
# --- 428m_crossmodel_all600_common.sh
N_RANDOM=1000 DISJOINT=1 MODELS=MaxToki-217M-L2:ctx217m600:2,MaxToki-217M:ctx217m600:4,MaxToki-1B-L4:ctx1b:4,MaxToki-1B-L7:ctx1b:7,scGPT-L2:ctx_scgpt600:2,scGPT-L4:ctx_scgpt600:4,STATE-L3:ctx_state600:3,STATE-L6:ctx_state600:6,MaxToki-217M-random:ctxrand:4,scGPT-random:ctx_scgpt600rand:4,STATE-random:ctx_state600rand:6,expression-landmarks:ctxexpr:0,expression-PCs:ctxexprpc:0,expression-centroids:ctxexprcen:0 OUTTAG=all600_common COMMON=1 $PY -u ctx_cross_model.py
# --- 428n_ctxexprpc_covexpr.sh
env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_functional_axes.py && env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v2.py && env PREFIX=ctxexprpc TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v3.py
# --- 429a_position_headline.sh
$PY -u ctx_position_confound.py
# --- 429b_position_600.sh
for P in ctx217m600:4 ctx1b:7 ctx_scgpt600:4 ctx_state600:6; do pre=${P%%:*}; t=${P#*:}; env PREFIX=$pre TAPS=$t TAP=$t RESEED_PER_TAP=1 $PY -u ctx_position_confound.py; done
# --- 429c_curated_headline.sh
$PY -u ctx_curated_targets.py
# --- 429d_curated_600.sh
for P in ctx217m600:4,2 ctx1b:7,4 ctx_scgpt600:4,2,9 ctx_state600:6,3,12 ctxexpr:0 ctxexprpc:0 ctxexprcen:0 ctx_state600rand:6 ctx_scgpt600rand:4 ctxrand:4; do pre=${P%%:*}; t=${P#*:}; t0=${t%%,*}; env PREFIX=$pre TAPS=$t TAP=$t0 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py; done
# --- 429e_split_217m600_nomassive.sh
SPLITS=30 MODEL=217m PREFIX=ctx217m600 AXIS=nuc_surf SITE=3 ALPHA_UNIT=centred ZERO_MASSIVE=1 $PY -u ctx_causal.py
# --- 429f_curated_headline_v2.sh
$PY -u ctx_curated_targets.py
# --- 429g_curated_600_v2.sh
for P in ctx217m600:4,2 ctx1b:7,4 ctx_scgpt600:4,2,9 ctx_state600:6,3,12 ctxexpr:0 ctxexprpc:0 ctxexprcen:0 ctx_state600rand:6 ctx_scgpt600rand:4 ctxrand:4; do pre=${P%%:*}; t=${P#*:}; t0=${t%%,*}; env PREFIX=$pre TAPS=$t TAP=$t0 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py; done
# --- 429h_curated_headline_repeat.sh
cp results/ctx_curated_targets.json results/_curated_run1.json && $PY -u ctx_curated_targets.py && cmp results/ctx_curated_targets.json results/_curated_run1.json && echo DETERMINISTIC
# --- 430_split_217m600_nuc.sh
SPLITS=30 MODEL=217m PREFIX=ctx217m600 AXIS=nuc_surf SITE=3 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 431_split_217m600_mito.sh
SPLITS=30 MODEL=217m PREFIX=ctx217m600 AXIS=mito_cyto SITE=3 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 432_split_217m600_trans.sh
SPLITS=30 MODEL=217m PREFIX=ctx217m600 AXIS=trans_transport SITE=3 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 433_split_1b_nuc.sh
SPLITS=30 MODEL=1b PREFIX=ctx1b AXIS=nuc_surf SITE=6 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 434_split_1b_nuc_nomassive.sh
SPLITS=30 MODEL=1b PREFIX=ctx1b AXIS=nuc_surf SITE=6 ALPHA_UNIT=centred ZERO_MASSIVE=1 $PY -u ctx_causal.py
# --- 435_split_1b_mito.sh
SPLITS=30 MODEL=1b PREFIX=ctx1b AXIS=mito_cyto SITE=6 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 436_split_1b_trans.sh
SPLITS=30 MODEL=1b PREFIX=ctx1b AXIS=trans_transport SITE=6 ALPHA_UNIT=centred $PY -u ctx_causal.py
# --- 437_split_head_nucL4.sh
SPLITS=30 AXIS=nuc_surf SITE=3 $PY -u ctx_causal.py
# --- 438_split_head_nucL8.sh
SPLITS=30 AXIS=nuc_surf SITE=7 $PY -u ctx_causal.py
# --- 439_split_head_mitoL4.sh
SPLITS=30 AXIS=mito_cyto SITE=3 $PY -u ctx_causal.py
# --- 440_split_head_transL4.sh
SPLITS=30 AXIS=trans_transport SITE=3 $PY -u ctx_causal.py
# --- 450a_cen50_extract.sh
MODE=centroid N_LANDMARK=50 OUTPREFIX=ctxexprcen50 $PY -u ctx_extract_expression.py
# --- 450b_cen50_polysemy.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 450c_cen50_position.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_position_confound.py
# --- 450d_cen50_mc6_funcaxes.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 450e_cen50_mc6_v2.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v2.py
# --- 450f_cen50_mc6_v3.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v3.py
# --- 450g_cen50_mc6_tightness.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_tightness_null.py
# --- 450h_cen50_mc6_directional.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_directional_probe.py
# --- 450i_cen50_curated.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 450j_cen50_covexpr.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_functional_axes.py && env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v2.py && env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v3.py
# --- 450k_cen50_mc9.sh
env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_functional_axes.py && env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v2.py && env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v3.py && env PREFIX=ctxexprcen50 TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_tightness_null.py
# --- 450l_crossmodel_all600.sh
N_RANDOM=1000 DISJOINT=1 MODELS=MaxToki-217M-L2:ctx217m600:2,MaxToki-217M:ctx217m600:4,MaxToki-1B-L4:ctx1b:4,MaxToki-1B-L7:ctx1b:7,scGPT-L2:ctx_scgpt600:2,scGPT-L4:ctx_scgpt600:4,STATE-L3:ctx_state600:3,STATE-L6:ctx_state600:6,MaxToki-217M-random:ctxrand:4,scGPT-random:ctx_scgpt600rand:4,STATE-random:ctx_state600rand:6,expression-landmarks:ctxexpr:0,expression-PCs:ctxexprpc:0,expression-centroids:ctxexprcen50:0,expression-centroids-512:ctxexprcen:0 OUTTAG=all600 $PY -u ctx_cross_model.py
# --- 450m_crossmodel_all600_common.sh
N_RANDOM=1000 DISJOINT=1 MODELS=MaxToki-217M-L2:ctx217m600:2,MaxToki-217M:ctx217m600:4,MaxToki-1B-L4:ctx1b:4,MaxToki-1B-L7:ctx1b:7,scGPT-L2:ctx_scgpt600:2,scGPT-L4:ctx_scgpt600:4,STATE-L3:ctx_state600:3,STATE-L6:ctx_state600:6,MaxToki-217M-random:ctxrand:4,scGPT-random:ctx_scgpt600rand:4,STATE-random:ctx_state600rand:6,expression-landmarks:ctxexpr:0,expression-PCs:ctxexprpc:0,expression-centroids:ctxexprcen50:0,expression-centroids-512:ctxexprcen:0 OUTTAG=all600_common COMMON=1 $PY -u ctx_cross_model.py
# --- 451a_anisotropy_allgenes.sh
$PY -u ctx_anisotropy.py
# --- 451b_switcher_rerun.sh
$PY -u ctx_switcher_test.py
# --- 452a_cen50c_extract.sh
MODE=centroid N_LANDMARK=50 CELL_CAP=50 OUTPREFIX=ctxexprcen50c $PY -u ctx_extract_expression.py
# --- 452b_cen50c_polysemy.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_polysemy.py
# --- 452c_cen50c_position.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_position_confound.py
# --- 452d_cen50c_mc6_funcaxes.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_functional_axes.py
# --- 452e_cen50c_mc6_v2.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v2.py
# --- 452f_cen50c_mc6_v3.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_coexpr_null_v3.py
# --- 452g_cen50c_mc6_tightness.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_tightness_null.py
# --- 452h_cen50c_mc6_directional.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 $PY -u ctx_directional_probe.py
# --- 452i_cen50c_curated.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 $PY -u ctx_curated_targets.py
# --- 452j_cen50c_covexpr.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_functional_axes.py && env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v2.py && env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=6 COV=expr $PY -u ctx_coexpr_null_v3.py
# --- 452k_cen50c_mc9.sh
env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_functional_axes.py && env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v2.py && env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_coexpr_null_v3.py && env PREFIX=ctxexprcen50c TAPS=0 TAP=0 RESEED_PER_TAP=1 MIN_CTX=9 $PY -u ctx_tightness_null.py
# --- 452l_crossmodel_all600.sh
N_RANDOM=1000 DISJOINT=1 MODELS=MaxToki-217M-L2:ctx217m600:2,MaxToki-217M:ctx217m600:4,MaxToki-1B-L4:ctx1b:4,MaxToki-1B-L7:ctx1b:7,scGPT-L2:ctx_scgpt600:2,scGPT-L4:ctx_scgpt600:4,STATE-L3:ctx_state600:3,STATE-L6:ctx_state600:6,MaxToki-217M-random:ctxrand:4,scGPT-random:ctx_scgpt600rand:4,STATE-random:ctx_state600rand:6,expression-landmarks:ctxexpr:0,expression-PCs:ctxexprpc:0,expression-centroids:ctxexprcen50c:0,expression-centroids-allcells:ctxexprcen50:0,expression-centroids-512:ctxexprcen:0 OUTTAG=all600 $PY -u ctx_cross_model.py
# --- 452m_crossmodel_all600_common.sh
N_RANDOM=1000 DISJOINT=1 MODELS=MaxToki-217M-L2:ctx217m600:2,MaxToki-217M:ctx217m600:4,MaxToki-1B-L4:ctx1b:4,MaxToki-1B-L7:ctx1b:7,scGPT-L2:ctx_scgpt600:2,scGPT-L4:ctx_scgpt600:4,STATE-L3:ctx_state600:3,STATE-L6:ctx_state600:6,MaxToki-217M-random:ctxrand:4,scGPT-random:ctx_scgpt600rand:4,STATE-random:ctx_state600rand:6,expression-landmarks:ctxexpr:0,expression-PCs:ctxexprpc:0,expression-centroids:ctxexprcen50c:0,expression-centroids-allcells:ctxexprcen50:0,expression-centroids-512:ctxexprcen:0 OUTTAG=all600_common COMMON=1 $PY -u ctx_cross_model.py

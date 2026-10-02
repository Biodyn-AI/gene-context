#!/bin/bash
# Full re-run after the MaxToki input-encoding fix (ctx_tokenise.py, 1 Oct 2026). Serial on purpose: the box has
# ~14 GB of real headroom and is shared. Old tensors were moved to results_pre_tokfix/npz/ so a failed extraction
# cannot be silently replaced by stale data. Logs: logs_tokfix/<step>.log; master log: logs_tokfix/master.log
cd "$(dirname "$0")"; PY=../../.venv_state/bin/python; mkdir -p logs_tokfix
M=logs_tokfix/master.log
avail_gb() { vm_stat | awk '/page size of/ {ps=$8} /Pages free/ {f=$3} /Pages inactive/ {i=$3} /Pages speculative/ {s=$3} END {printf "%d", (f+i+s)*ps/1e9}'; }
wait_mem() {  # wait (up to ~2 h) until at least $1 GB is free+inactive
  for k in $(seq 1 240); do a=$(avail_gb); [ "$a" -ge "$1" ] && return 0; echo "  [mem] ${a} GB available < $1 GB; waiting" >> $M; sleep 30; done; return 0; }
step() {  # step <name> <min_free_gb> <fatal 0/1> <env/cmd...>
  name=$1; need=$2; fatal=$3; shift 3
  wait_mem $need
  echo "[start] $name $(date '+%F %T') (avail $(avail_gb) GB)" >> $M
  env "$@" > logs_tokfix/$name.log 2>&1; rc=$?
  echo "[end]   $name $(date '+%F %T') exit $rc" >> $M
  if [ $rc -ne 0 ] && [ $fatal -eq 1 ]; then echo "[ABORT] fatal step $name failed" >> $M; exit 1; fi
}
echo "[chain start] $(date '+%F %T')" > $M
# A1: main extraction (217M, 12 cell types x 1000 cells, 6000 genes, layers 0/4/8/11) + cap-20 control at layer 4
step A1_extract_main 7 1 TAPS=0,4,8,11 MAXGENES=6000 CELLSCTX=1000 EXTRA_CAP=20 EXTRA_TAPS=4 EXTRA_PREFIX=ctx217m_cap20 $PY -u ctx_extract_maxtoki.py
# analyses on the main extraction
step B_polysemy 4 0 $PY -u ctx_polysemy.py
step B_anisotropy 4 0 $PY -u ctx_anisotropy.py
step B_position 4 0 $PY -u ctx_position_confound.py
step B_funcaxes 4 0 $PY -u ctx_functional_axes.py
step B_coexpr_v2 4 0 $PY -u ctx_coexpr_null_v2.py
step B_coexpr_v3 4 0 $PY -u ctx_coexpr_null_v3.py
step B_tightness 4 0 $PY -u ctx_tightness_null.py
step B_directional 4 0 $PY -u ctx_directional_probe.py
step B_curated 4 0 $PY -u ctx_curated_targets.py
step B_switcher 4 0 $PY -u ctx_switcher_test.py
step B_loadings 4 0 $PY -u ctx_celltype_loadings.py
step B_causal_nuc_L4 5 0 AXIS=nuc_surf SITE=3 $PY -u ctx_causal.py
step B_causal_nuc_L8 5 0 AXIS=nuc_surf SITE=7 $PY -u ctx_causal.py
step B_causal_mito_L4 5 0 AXIS=mito_cyto SITE=3 $PY -u ctx_causal.py
step B_causal_trans_L4 5 0 AXIS=trans_transport SITE=3 $PY -u ctx_causal.py
step B_predlink 5 0 $PY -u ctx_prediction_link.py
# scaling pair + random-weights control (same 600 cells per type, same 5000-gene panel)
step A3_extract_217m600 6 1 MODEL=MaxToki-217M-HF CELLSCTX=600 MAXGENES=5000 TAPS=2,4 OUTPREFIX=ctx217m600 $PY -u ctx_extract_maxtoki.py
step A4_extract_1b600 10 1 MODEL=MaxToki-1B-HF CELLSCTX=600 MAXGENES=5000 TAPS=4,7 OUTPREFIX=ctx1b $PY -u ctx_extract_maxtoki.py
step A5_extract_random 6 1 $PY -u ctx_extract_random.py
step B_crossmodel 5 0 $PY -u ctx_cross_model.py
# independent dataset (full Setty data)
step A6_extract_setty 6 1 $PY -u ctx_extract_devel.py
step B_independent 4 0 $PY -u ctx_independent.py
# full depth profile: all 12 hidden states, 4000 genes, in two passes of 6 layers (memory)
step A2_scan_a 7 1 TAPS=0,1,2,3,4,5 MAXGENES=4000 CELLSCTX=1000 OUTPREFIX=ctxscan $PY -u ctx_extract_maxtoki.py
step A2_scan_b 7 1 TAPS=6,7,8,9,10,11 MAXGENES=4000 CELLSCTX=1000 OUTPREFIX=ctxscan $PY -u ctx_extract_maxtoki.py
step B_layercurve 4 0 $PY -u ctx_layer_curve.py
echo "[chain done] $(date '+%F %T')" >> $M
# robustness (run separately on 2 Oct 2026): Setty with spliced counts only
step A6b_setty_spliced 6 0 SETTY_COUNTS=spliced $PY -u ctx_extract_devel.py
step B_independent_spliced 4 0 PREFIX=ctx_devel_spliced $PY -u ctx_independent.py

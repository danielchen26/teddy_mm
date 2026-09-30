#!/usr/bin/env bash
# Rerun the TEDDY -> head -> ANM chain on an embedding made with the official TEDDY-G
# preprocessing (scripts/03_embed_rna.py --preprocessing official). Nothing existing is
# overwritten: every output goes to a new *_official (or *_legacy_s0) folder.
#
# Resumable: each step is skipped when its output already exists, and the embedding is written in
# shards that a rerun skips. RERUN_BUDGET_SEC stops the chain (exit 75) before a step once that many
# seconds have passed, so it can be driven in short calls; rerun the same command to continue.
#
# usage: bash scripts/rerun_official.sh [step ...]      (default: all steps, in order)
# steps: medians legacy_check pilot embed mode_a train eval export hard_proof exactness scope_refine
#        depth nkt examples pooling_first
#   mode_a runs twice in the default order: after embed (layer probes, CPU only) and after train (adds
#   the head reference). It needs the embedding made with EMBED_EXTRA=--save-layer-means (see step_embed)
#   and can run alongside 'train' in a second shell: bash scripts/rerun_official.sh mode_a
#
# env (defaults in brackets):
#   MAIN      main checkout holding data/ and outputs/ [derived from the worktree's .git file, else repo root]
#   PY        python with torch, numpy, scipy, sklearn, tqdm, safetensors [$MAIN/.venv312/bin/python]
#   ANM_ROOT  ANM repo for the bridge proofs [$MAIN/../ANM]
#   CKPT      TEDDY-G 70M checkpoint folder [$MAIN/../teddy_mwe/ckpt/teddy_g_70M]
#   EXP_ROOT  folder with the scratch experiments redesign/experiment-2, redesign/experiment-6 and deck/
#             (depth, nkt and examples steps are skipped without it)
#   WORK      scratch folder for checks, pilot and experiment outputs [$MAIN/.tmp/official_rerun]
#   DEVICE    [mps]   THREADS [6]   WORKERS (ANM processes) [4]   SHARD [5000]
#   MODE_A_THREADS  CPU threads for mode_a [3]
#   P95_MANIFEST  export manifest whose p95_train the depth step uses
#                 [$MAIN/outputs/anm_cite_bridge/missing_modality/export_manifest.json]
#   SCHEMA_O0     O0 ANM schema for the nkt step [$MAIN/outputs/anm_cite_bridge/schema_O0.json]
#   TRAIN_EXTRA / DEPTH_EXTRA   extra args for 04_train.py / s1_thin_embed.py (dry runs: "--epochs 2", "--n 20")
#
# Sub-steps can be named on their own (for running lanes in parallel): train_official train_legacy_s0
# eval_official eval_legacy_s0 eval_repro depth_embed depth_anm nkt_predict nkt_k10 nkt_k10_nogdt nkt_diag
# nkt_k30 nkt_report.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ -z "${MAIN:-}" ]; then
  if [ -f "$ROOT/.git" ]; then
    gd="$(sed -n 's/^gitdir: //p' "$ROOT/.git")"; MAIN="${gd%%/.git/*}"
  else
    MAIN="$ROOT"
  fi
fi
PY="${PY:-$MAIN/.venv312/bin/python}"
export ANM_ROOT="${ANM_ROOT:-$MAIN/../ANM}"
CKPT="${CKPT:-$MAIN/../teddy_mwe/ckpt/teddy_g_70M}"
WORK="${WORK:-$MAIN/.tmp/official_rerun}"
DEVICE="${DEVICE:-mps}"
THREADS="${THREADS:-6}"
WORKERS="${WORKERS:-4}"
SHARD="${SHARD:-5000}"
MODE_A_THREADS="${MODE_A_THREADS:-3}"
BUDGET="${RERUN_BUDGET_SEC:-0}"
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" VECLIB_MAXIMUM_THREADS="$THREADS"
export PYTHONUNBUFFERED=1

MEDIANS="$MAIN/data/reference/teddy_gene_medians.json"
L="$MAIN/data/processed/cite"                  # legacy embedding (stored, read-only here)
P="$MAIN/data/processed/cite_official"         # official embedding
H="$MAIN/outputs/cite_phase1_official"         # head trained on the official z, seed 0
HL="$MAIN/outputs/cite_phase1_legacy_s0"       # head trained on the legacy z, seed 0 (control)
HOLD="$MAIN/outputs/cite_phase1"               # legacy head behind the current pages (unseeded)
B="$MAIN/outputs/anm_cite_bridge_official"     # bridge export + proofs on the official head
MA="$MAIN/outputs/mode_a_official"             # Mode A layer probes on the official embedding
# inputs the scratch experiments otherwise read from hard-coded /Users/.../teddy_mm/outputs paths
# (the main checkout's copy when present, else the committed copy in this worktree)
first_existing() { for f in "$@"; do if [ -e "$f" ]; then echo "$f"; return; fi; done; echo "$1"; }
P95_MANIFEST="${P95_MANIFEST:-$(first_existing "$MAIN/outputs/anm_cite_bridge/missing_modality/export_manifest.json" \
    "$ROOT/outputs/anm_cite_bridge/missing_modality/export_manifest.json")}"
SCHEMA_O0="${SCHEMA_O0:-$(first_existing "$MAIN/outputs/anm_cite_bridge/schema_O0.json" "$ROOT/outputs/anm_cite_bridge/schema_O0.json")}"
read -r -a TRAIN_X <<< "${TRAIN_EXTRA:-}"
read -r -a DEPTH_X <<< "${DEPTH_EXTRA:-}"
LOG="$WORK/logs"
mkdir -p "$WORK" "$LOG"
T0=$(date +%s)

say() { echo "[$(date '+%F %T')] $*"; }
budget_left() {
  if [ "$BUDGET" = "0" ]; then echo 0; return; fi
  echo $(( BUDGET - ($(date +%s) - T0) ))
}
check_budget() {
  if [ "$BUDGET" != "0" ] && [ "$(budget_left)" -le 0 ]; then
    say "budget of ${BUDGET}s used; rerun to resume"; exit 75
  fi
}
# run NAME MARKER cmd... : skip when MARKER exists, else run and log; stop the chain on failure
run() {
  local name="$1" marker="$2"; shift 2
  if [ -e "$marker" ]; then say "skip $name ($marker exists)"; return 0; fi
  check_budget
  say "run  $name"
  local t=$(date +%s)
  "$@" >>"$LOG/$name.log" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then say "FAIL $name (exit $rc, see $LOG/$name.log)"; tail -20 "$LOG/$name.log"; exit $rc; fi
  if [ ! -e "$marker" ]; then say "FAIL $name: $marker not written (see $LOG/$name.log)"; exit 1; fi
  say "done $name in $(( $(date +%s) - t ))s"
}
# embed OUT_DIR args... : call 03_embed_rna.py until every shard exists (exit 75 = budget, resume)
embed() {
  local out="$1"; shift
  if [ -e "$out/z_rna.npy" ]; then say "skip embed -> $out (z_rna.npy exists)"; return 0; fi
  check_budget
  say "run  embed -> $out"
  local extra=()
  if [ "$BUDGET" != "0" ]; then extra=(--time-budget-sec "$(budget_left)"); fi
  "$PY" "$ROOT/scripts/03_embed_rna.py" --out-dir "$out" --ckpt "$CKPT" --device "$DEVICE" --threads "$THREADS" ${extra[@]+"${extra[@]}"} "$@" \
      >>"$LOG/embed_$(basename "$out").log" 2>&1
  local rc=$?
  if [ $rc -eq 75 ]; then say "embed -> $out paused (budget); rerun to resume"; exit 75; fi
  if [ $rc -ne 0 ]; then say "FAIL embed -> $out (exit $rc)"; tail -20 "$LOG/embed_$(basename "$out").log"; exit $rc; fi
  say "done embed -> $out"
}

step_medians() { run medians "$MEDIANS" bash "$ROOT/scripts/01b_fetch_teddy_medians.sh" "$MEDIANS"; }

step_legacy_check() {
  # legacy path must reproduce the stored z: 208 contiguous cells (13 batches of 16, the old batching)
  export PYTORCH_ENABLE_MPS_FALLBACK=1
  embed "$WORK/legacy_check" --preprocessing legacy --processed "$L" --start 4000 --stop 4208 --shard-size 0
  unset PYTORCH_ENABLE_MPS_FALLBACK
  run legacy_check_cmp "$WORK/legacy_check/compare.json" "$PY" "$ROOT/scripts/rerun_checks.py" legacy \
      --run "$WORK/legacy_check" --legacy "$L"
}

step_pilot() {
  # 3,000 random cells, fixed padding to 2048 so all three poolings come from one forward pass
  embed "$WORK/pilot" --preprocessing official --processed "$L" --medians "$MEDIANS" --max-cells 3000 --cell-seed 0 \
      --save-poolings first,all-positions --shard-size 1000
  # descriptive comparison with the legacy z, plus a length-bucketing check on 64 pilot cells
  run pilot_cmp "$WORK/pilot/compare.json" "$PY" "$ROOT/scripts/rerun_checks.py" pilot --run "$WORK/pilot" --legacy "$L" \
      --medians "$MEDIANS" --ckpt "$CKPT" --device "$DEVICE"
}

step_embed() {
  # shards live in $WORK (scratch); only the assembled z_rna*.npy + manifest go to $P.
  # These are the settings of the official run: fp16 (matches fp32, cosine >= 0.999999), batch 32, and
  # --save-layer-means for Mode A (about 1.1 GB float16 for 12 layers x 90,261 cells).
  embed "$P" --preprocessing official --processed "$L" --medians "$MEDIANS" --seq-len 2048 --pooling gene-mean \
      --save-poolings first --length-buckets --batch-size 32 --autocast fp16 --save-layer-means --shard-size "$SHARD" --shard-dir "$WORK/shards_official" \
      ${EMBED_EXTRA:-}
}

step_mode_a() {
  # Mode A: linear probes at every TEDDY depth (input embedding, layers 1..L, z) for the NK-T difference,
  # read against the head (scripts/mode_a_layer_probes.py; CPU only, MODE_A_THREADS threads).
  # 1st call (after embed): the probes -> $MA/mode_a_probes.json. Once the official head has finished
  # training ($H/metrics.json), a later call adds the head reference and reuses the cached probes.
  if [ ! -e "$P/z_rna_layer_means.npy" ]; then
    say "FAIL mode_a: $P/z_rna_layer_means.npy missing (the embedding needs EMBED_EXTRA=--save-layer-means)"; return 1
  fi
  local args=(--processed "$L" --embed-dir "$P" --ckpt "$CKPT" --medians "$MEDIANS" --out-dir "$MA" --threads "$MODE_A_THREADS")
  run mode_a_probes "$MA/mode_a_probes.json" "$PY" "$ROOT/scripts/mode_a_layer_probes.py" "${args[@]}"
  if [ -e "$H/metrics.json" ]; then
    run mode_a_head "$MA/mode_a_head_reference.json" "$PY" "$ROOT/scripts/mode_a_layer_probes.py" "${args[@]}" \
        --head-ckpt "$H/best.pt"
  else
    say "mode_a: head $H not trained yet; run 'mode_a' again after 'train' to add the head reference"
  fi
}

step_train_official() {
  run train_official "$H/metrics.json" "$PY" "$ROOT/scripts/04_train.py" --processed "$P" --out "$H" --seed 0 --device "$DEVICE" \
      ${TRAIN_X[@]+"${TRAIN_X[@]}"}
}
step_train_legacy_s0() {
  run train_legacy_s0 "$HL/metrics.json" "$PY" "$ROOT/scripts/04_train.py" --processed "$L" --out "$HL" --seed 0 --device "$DEVICE" \
      ${TRAIN_X[@]+"${TRAIN_X[@]}"}
}
step_train() { step_train_official; step_train_legacy_s0; }

step_eval_official() {
  run eval_official "$H/test_eval_meta.json" "$PY" "$ROOT/scripts/05_eval.py" --processed "$P" --ckpt "$H/best.pt" --out-dir "$H" \
      --size-factor train-median --seed 0 --device "$DEVICE"
}
step_eval_legacy_s0() {
  run eval_legacy_s0 "$HL/test_eval_meta.json" "$PY" "$ROOT/scripts/05_eval.py" --processed "$L" --ckpt "$HL/best.pt" --out-dir "$HL" \
      --size-factor train-median --seed 0 --device "$DEVICE"
}
step_eval_repro() {
  # sanity: the legacy head behind the pages must still give mean MLP Pearson 0.6033
  run eval_repro "$WORK/eval_repro/test_eval_meta.json" "$PY" "$ROOT/scripts/05_eval.py" --processed "$L" --ckpt "$HOLD/best.pt" \
      --out-dir "$WORK/eval_repro" --size-factor train-median --seed 0 --device "$DEVICE"
}
step_eval() { step_eval_official; step_eval_legacy_s0; step_eval_repro; }

step_export() {
  # --device: the exporter's own default is auto (MPS when present), so pass DEVICE through
  run export "$B/export_manifest.json" "$PY" "$ROOT/bridge_anm/export_cite_events.py" --processed "$P" --ckpt "$H/best.pt" \
      --per-protein "$H/test_per_protein.json" --metrics "$H/metrics.json" --out-dir "$B" --device "$DEVICE"
}

step_hard_proof() {
  run hard_proof "$B/hard_proof/hard_proof_results.json" "$PY" "$ROOT/bridge_anm/run_hard_proof.py" --bridge-dir "$B" \
      --out-dir "$B/hard_proof" --workers "$WORKERS"
}

step_exactness() {
  # ANM vs the fixed rule, cell by cell, and out-of-scope call shares (same computation as the v2 run)
  run exactness "$B/hard_proof/exactness_matched_coverage_official.json" "$PY" "$ROOT/bridge_anm/run_exactness_check.py" \
      "$B" "$B/hard_proof/exactness_matched_coverage_official.json" --workers "$WORKERS"
}

step_scope_refine() {
  run scope_refine "$B/scope_refine/scope_refine_results.json" "$PY" "$ROOT/bridge_anm/run_scope_refine_proof.py" \
      --bridge-dir "$B" --out-dir "$B/scope_refine" --skip-adt-only --workers "$WORKERS"
}

need_exp() { if [ -z "${EXP_ROOT:-}" ]; then say "skip $1 (EXP_ROOT not set)"; return 1; fi; }

depth_envs() {
  E2="$EXP_ROOT/redesign/experiment-2"; O2="$WORK/experiment-2"
  envs2=(env TEDDY_MM_CODE="$ROOT" TEDDY_PROCESSED="$P" TEDDY_PHASE1="$H" EXP2_OUT="$O2" TEDDY_CKPT="$CKPT"
         P95_MANIFEST="$P95_MANIFEST")
  mkdir -p "$O2"
}
step_depth_embed() {
  # re-embeds thinned RNA with TEDDY: same preprocessing, medians, seq_len and autocast as the
  # official z (s1_thin_embed.py --autocast match reads them from $P/z_rna_manifest.json and stops on
  # a mismatch); length-bucketed batches as in 03_embed_rna.py --length-buckets
  need_exp depth || return 0
  depth_envs
  run depth_embed "$O2/cells_val.npz" "${envs2[@]}" "$PY" "$E2/s1_thin_embed.py" --preprocessing official --seq-len 2048 \
      --medians "$MEDIANS" --autocast match --device "$DEVICE" --batch 16 ${DEPTH_X[@]+"${DEPTH_X[@]}"}
}
step_depth_anm() {
  need_exp depth || return 0
  depth_envs
  run depth_trust "$O2/results.json" "${envs2[@]}" "$PY" "$E2/s2_trust_anm.py" --skip-phase2
  # stdout goes to a temp file first so a failed run leaves no marker behind
  run depth_summary "$O2/summary.txt" bash -c "cd '$E2' && ${envs2[*]} '$PY' s3_summary.py > '$O2/summary.txt.tmp' && mv '$O2/summary.txt.tmp' '$O2/summary.txt'"
  run depth_boot "$O2/pooled_bootstrap.json" bash -c "cd '$E2' && ${envs2[*]} '$PY' s4_pooled_boot.py"
}
step_depth() { step_depth_embed; step_depth_anm; }

nkt_envs() {
  E6="$EXP_ROOT/redesign/experiment-6"; O6="$WORK/experiment-6"
  envs6=(env TEDDY_MM_CODE="$ROOT" TEDDY_PROCESSED="$P" TEDDY_PHASE1="$H" BRIDGE_DIR="$B" EXP6_OUT="$O6"
         SCHEMA_O0="$SCHEMA_O0")
  mkdir -p "$O6"
}
# nkt sub-steps: after nkt_predict, nkt_k10 / nkt_k10_nogdt / nkt_k30 / nkt_diag are independent (single-core
# ANM loops) and can run as parallel lanes; nkt_report needs nkt_k10.
step_nkt_predict() { need_exp nkt || return 0; nkt_envs
  run nkt_predict "$O6/pred_all.npy" "${envs6[@]}" "$PY" "$E6/s1_predict.py"; }
step_nkt_k10() { need_exp nkt || return 0; nkt_envs
  run nkt_k10 "$O6/exp6_results_k10.json" "${envs6[@]}" "$PY" "$E6/s2_lookalike.py" 10 main; }
step_nkt_k10_nogdt() { need_exp nkt || return 0; nkt_envs
  run nkt_k10_nogdt "$O6/exp6_results_k10_no_gdT158.json" "${envs6[@]}" "$PY" "$E6/s2_lookalike.py" 10 no_gdT158; }
step_nkt_diag() { need_exp nkt || return 0; nkt_envs
  run nkt_gap "$O6/pair_channel_gap.json" "${envs6[@]}" "$PY" "$E6/s4_pair_channel_gap.py"
  run nkt_gdt "$O6/gdt_check.json" "${envs6[@]}" "$PY" "$E6/s6_gdt_check.py"
  run nkt_gap_nogdt "$O6/pair_gap_nogdt.json" "${envs6[@]}" "$PY" "$E6/s8_pair_gap_nogdt.py"
  run nkt_random "$O6/random_pair_control.json" "${envs6[@]}" "$PY" "$E6/s9_random_pair_control.py"; }
step_nkt_k30() { need_exp nkt || return 0; nkt_envs
  run nkt_k30 "$O6/exp6_results_k30.json" "${envs6[@]}" "$PY" "$E6/s2_lookalike.py" 30 main; }
step_nkt_report() { need_exp nkt || return 0; nkt_envs
  # stdout goes to a temp file first so a failed run leaves no marker behind
  run nkt_tables "$O6/tables_k10.md" bash -c "${envs6[*]} '$PY' '$E6/s3_tables.py' k10 > '$O6/tables_k10.md.tmp' && mv '$O6/tables_k10.md.tmp' '$O6/tables_k10.md'"
  run nkt_donor "$O6/per_donor_k10.md" bash -c "${envs6[*]} '$PY' '$E6/s7_per_donor.py' k10 s4_d0.10 > '$O6/per_donor_k10.md.tmp' && mv '$O6/per_donor_k10.md.tmp' '$O6/per_donor_k10.md'"; }
step_nkt() { step_nkt_predict; step_nkt_k10; step_nkt_k10_nogdt; step_nkt_diag; step_nkt_k30; step_nkt_report; }

step_examples() {
  need_exp examples || return 0
  run examples "$WORK/examples_official.json" env WT="$ROOT" V2="$B" OUT_JSON="$WORK/examples_official.json" \
      FIXED_IDS="${FIXED_IDS:-cite_site4_85386,cite_site4_75462}" "$PY" "$EXP_ROOT/deck/examples.py"
}

step_pooling_first() {
  # pooling sensitivity: the same head recipe on token 0 of the same forward pass
  local PF="$WORK/cite_official_first" HF="$MAIN/outputs/cite_phase1_official_first"
  if [ ! -e "$P/z_rna_first.npy" ]; then say "FAIL pooling_first: $P/z_rna_first.npy missing (embedding not finished?)"; return 1; fi
  if [ ! -e "$PF/z_rna.npy" ]; then
    mkdir -p "$PF"; ln -sf "$P/z_rna_first.npy" "$PF/z_rna.npy"
    ln -sf "$L/cite_arrays.npz" "$PF/cite_arrays.npz"; ln -sf "$L/meta.json" "$PF/meta.json"
  fi
  run train_first "$HF/metrics.json" "$PY" "$ROOT/scripts/04_train.py" --processed "$PF" --out "$HF" --seed 0 --device "$DEVICE" \
      ${TRAIN_X[@]+"${TRAIN_X[@]}"}
  run eval_first "$HF/test_eval_meta.json" "$PY" "$ROOT/scripts/05_eval.py" --processed "$PF" --ckpt "$HF/best.pt" --out-dir "$HF" \
      --size-factor train-median --seed 0 --device "$DEVICE"
}

ALL=(medians legacy_check pilot embed mode_a train mode_a eval export hard_proof exactness scope_refine depth nkt examples pooling_first)
STEPS=("$@"); [ ${#STEPS[@]} -eq 0 ] && STEPS=("${ALL[@]}")
say "rerun_official: MAIN=$MAIN PY=$PY CKPT=$CKPT ANM_ROOT=$ANM_ROOT WORK=$WORK EXP_ROOT=${EXP_ROOT:-} DEVICE=$DEVICE steps=${STEPS[*]}"
say "  P95_MANIFEST=$P95_MANIFEST SCHEMA_O0=$SCHEMA_O0"
for s in "${STEPS[@]}"; do
  "step_$s" || { say "FAIL step $s"; exit 1; }
done
say "all requested steps done in $(( $(date +%s) - T0 ))s"

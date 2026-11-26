"""Generates run_c5b_eval.sh and run_c5b_memtrial.sh (VM wrappers) from the watchdog / trap / run_step blocks of the reviewed run_c5b_train_v2.sh (taken verbatim between the markers below) and the sha table of the scripts in this directory.
Resource limits fixed by 实验管理: Pss <= 12 GiB (12582912 kB), GPU <= 4 GiB (4096 MiB), 8 threads, nice 19; refuses to start if MemAvailable < 24 GiB (headroom for the CAN host diagnostics), /data >= 90%, or the output partition has < 5 GiB free.
Run from this directory (results_src/c5b); never overwrites."""
import hashlib, re
from pathlib import Path
T = Path("../../code/run_c5b_train_v2.sh").read_text()
assert hashlib.sha256(T.encode()).hexdigest().startswith("6429fc15709e8fb8"), "template changed"
FIN = "".join(l for l in T[T.index("finish() {"):T.index("OTHER=$(pgrep")].splitlines(True) if "${TAG}" not in l); WATCH = T[T.index("tree_pids() {"):T.index('cd "$SD" || exit 1')]
sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
A2BLOCK_PLACEHOLDER = Path("a2_block.inc").read_text()      # the eval_a2 step (inserted before the inference step of run_c5b_eval.sh)
EVAL = ["c5b_bbbc_rev1.py", "c5b_mh_rev1.py", "c5b_lh_rev1.py", "c5b_lbsame_rev1.py", "c5b_bbbc_s11.py", "c5b_mh_s11.py", "c5b_lh_s11.py"]
S8 = "54a2af8d3704c16b0994"      # s8_bounds.py sha prefix (results_src/mh/orig_iter1/s8_bounds.py)
HEAD = '''#!/usr/bin/env bash
# P5 C5b %(what)s wrapper (VM only). Limits fixed by 实验管理: Pss <= 12 GiB (12582912 kB) of the whole process tree, GPU <= 4 GiB (4096 MiB), 8 threads, nice 19; the evaluation steps run Pool(4) workers with OMP_NUM_THREADS=2 (4 x 2 = 8 threads).
# Refuses to start if MemAvailable < 24 GiB (headroom for the CAN host diagnostics), /data >= 90%%, the output partition has < 5 GiB free, a sha differs from the reviewed table, or the output directory exists. Watchdog / trap / run_step are the blocks of the reviewed run_c5b_train_v2.sh.
set -u
PY=<VM_HOME>/miniconda3/envs/pyproject/bin/python
SD=${C5B_SCRIPTS:-<VM_HOME>/P5_StainRef/results_src_c5b}        # directory holding the scripts below and s8_bounds.py
CK=<VM_HOME>/P5_StainRef/runs/c5b/ckpt
ARCH=<VM_HOME>/P5_StainRef/runs/iter1_20261021/infer            # archived A2 inputs: ids_iter1.npz (read only)
IDS_MH=${IDS_MH:-<VM_HOME>/P5_StainRef/runs/mpidb_hc_20261025T130513Z/infer/ids_mh.npz}
IDS_LH=${IDS_LH:-<VM_HOME>/P5_StainRef/runs/lacuna_hc_20261025T221212Z/infer/ids_lh.npz}
IDS_LBS=${IDS_LBS:-<VM_HOME>/P5_StainRef/runs/lacuna_samepop_paired_20261026T092505Z/data/ids_lb_same.npz}
TS=$(date -u +%%Y%%m%%dT%%H%%M%%SZ)
OUT=<VM_HOME>/P5_StainRef/runs
OUTDIR=$OUT/%(dirname)s_$TS; RUNDIR=$OUTDIR; STEPDIR=$OUTDIR
RSS_LIMIT_KB=${RSS_LIMIT_KB:-12582912}
GPU_LIMIT_MIB=${GPU_LIMIT_MIB:-4096}
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8 P5_TORCH_THREADS=8 NPROC=4 P5_MODEL_PY=<DATA>/medprune/MedPrune/models/malarianet.py
USED=$(df --output=pcent /data | tail -1 | tr -dc '0-9'); [ "${USED:-100}" -ge 90 ] && { echo "/data usage ${USED}%% >= 90%%, not starting"; exit 1; }
AVAIL_KB=$(df -Pk "$OUT" | awk 'NR==2{print $4}'); [ "${AVAIL_KB:-0}" -lt 5242880 ] && { echo "output partition has ${AVAIL_KB:-0} KB free (< 5 GiB), not starting"; exit 1; }
MEMAV_KB=$(awk '/MemAvailable/{print $2}' /proc/meminfo); [ "${MEMAV_KB:-0}" -lt 25165824 ] && { echo "MemAvailable ${MEMAV_KB:-0} kB < 24 GiB, not starting (headroom for the CAN host diagnostics)"; exit 1; }

# 0) the invoked scripts must be exactly the reviewed versions
declare -A SHA=( %(shas)s )
for f in "${!SHA[@]}"; do
  [ "$(sha256sum "$SD/$f" | cut -d' ' -f1)" = "${SHA[$f]}" ] || { echo "sha mismatch: $f"; exit 1; }
done
[ "$(sha256sum "$SD/s8_bounds.py" | cut -c1-20)" = "%(s8)s" ] || { echo "s8_bounds.py sha mismatch"; exit 1; }
[ "$(ls "$CK"/*.pth 2>/dev/null | wc -l)" -eq 10 ] || { echo "expected exactly 10 C5b checkpoints in $CK"; exit 1; }
[ -e "$OUTDIR" ] && { echo "output dir exists: $OUTDIR"; exit 1; }
mkdir -p "$OUTDIR" || exit 1
LOG=$OUTDIR/run.log; WDLOG=$OUTDIR/watchdog.log
log() { echo "$(date -u +%%FT%%TZ) $*" | tee -a "$LOG"; }

'''
OTHER = '''OTHER=$(pgrep -af "p5_main.py|p5_c5b_train.py|p5_calib.py|p6_pathotrunk.py|lh_infer|c4_feat|run_c5b_train|c5b_infer|c5b_memtrial|c5b_[a-z]*_rev1|c5b_[a-z]*_s11" | grep -v pgrep | grep -v -e "%(self)s" || true)
if [ -n "$OTHER" ]; then log "another job of ours is running, refusing:"; echo "$OTHER" | tee -a "$LOG"; exit 1; fi
log "start $TS; /data ${USED}%%; MemAvailable ${MEMAV_KB} kB; GPU: $(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader); limits Pss ${RSS_LIMIT_KB} kB, GPU ${GPU_LIMIT_MIB} MiB; torch/cuda recorded by the infer step"
'''
shas_e = " ".join("[%s]=%s" % (f, sha(f)) for f in ["c5b_infer.py", "check_c5b.py", "c5b_a2_link.py", "c5b_a2_compare.py"] + EVAL)
ev = HEAD % dict(what="inference + evaluation", dirname="c5b_eval", shas=shas_e, s8=S8)
ev += FIN.replace("find . -type f \\( -name \"*.json\" -o -name \"*.npz\" \\) -print0", "find . -type f -print0") if False else FIN
ev += OTHER % dict(self="run_c5b_eval.sh") + WATCH + '''cd "$SD" || exit 1
'''+A2BLOCK_PLACEHOLDER+'''# 1) inference (GPU), 10 checkpoints x 5 sites; logits and ids go to $OUTDIR/infer
run_step infer "$PY" c5b_infer.py --out "$OUTDIR/infer" --ckpt-dir "$CK" --ids-iter1 "$ARCH/ids_iter1.npz" --ids-lb-same "$IDS_LBS" || exit $RC_LAST
echo "$(sha256sum "$CK"/*.pth | awk '{print $1" "$2}' | tr '\\n' ';')" > "$OUTDIR/ckpt_sha.txt"; log "checkpoint shas written to ckpt_sha.txt"
# 2) check: rows, fields, labels and field lists equal to the frozen A2 ids; all logits files present, shaped, finite (prints PASS or the failing keys)
run_step check "$PY" check_c5b.py --new-dir "$OUTDIR/infer" --ids-iter1 "$ARCH/ids_iter1.npz" --ids-mh "$IDS_MH" --ids-lh "$IDS_LH" --ids-lb-same "$IDS_LBS" || { log "ERROR: check failed: $(tr '\\n' ' ' < "$OUTDIR/check.stdout.log")"; exit 1; }
[ "$(tr -d '\\n' < "$OUTDIR/check.stdout.log")" = "PASS" ] || { log "ERROR: check did not print PASS"; exit 1; }
log "check: PASS"
# 3) evaluation (CPU only), one arm label after the other; ARMS = ['a2'] inside the scripts keeps ai = 0 (the frozen A2 streams)
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 C5B_LOGDIR="$OUTDIR/infer"
for LAB in c5b-f c5b-fs; do
  export C5B_ARM=$LAB; E=$OUTDIR/eval_$LAB
  run_step ${LAB}_bbbc_rev1 "$PY" c5b_bbbc_rev1.py "$ARCH" "$ARCH" "$ARCH" "$E/bbbc_rev1" || exit $RC_LAST
  run_step ${LAB}_mh_rev1 "$PY" c5b_mh_rev1.py "$ARCH" "$ARCH" "$OUTDIR/infer" "$E/mh_rev1" || exit $RC_LAST
  run_step ${LAB}_lh_rev1 "$PY" c5b_lh_rev1.py "$ARCH" "$ARCH" "$OUTDIR/infer" "$E/lh_rev1" || exit $RC_LAST
  run_step ${LAB}_lbsame_rev1 "$PY" c5b_lbsame_rev1.py "$ARCH" "$ARCH" "$OUTDIR/infer" "$E/lbsame_rev1" || exit $RC_LAST
  run_step ${LAB}_bbbc_s11 "$PY" c5b_bbbc_s11.py "$ARCH" "$ARCH" "$E/bbbc_s11" || exit $RC_LAST
  run_step ${LAB}_mh_s11 "$PY" c5b_mh_s11.py "$ARCH" "$OUTDIR/infer" "$E/mh_s11" || exit $RC_LAST
  run_step ${LAB}_lh_s11 "$PY" c5b_lh_s11.py "$ARCH" "$OUTDIR/infer" "$E/lh_s11" || exit $RC_LAST
  for f in bbbc_rev1/c5b_bbbc_rev1_results_$LAB.json bbbc_rev1/c5b_bbbc_rev1_perdraw_${LAB}_m22.npz mh_rev1/c5b_mh_rev1_results_$LAB.json mh_rev1/c5b_mh_rev1_perdraw_${LAB}_m22.npz lh_rev1/c5b_lh_rev1_results_$LAB.json lh_rev1/c5b_lh_rev1_perdraw_${LAB}_m22.npz lbsame_rev1/c5b_lbsame_rev1_results_$LAB.json lbsame_rev1/c5b_lbsame_rev1_perdraw_${LAB}_m22.npz bbbc_s11/c5b_bbbc_s11_results_$LAB.json bbbc_s11/c5b_bbbc_s11_perdraw_${LAB}_m22.npz mh_s11/c5b_mh_s11_results_$LAB.json mh_s11/c5b_mh_s11_perdraw_${LAB}_m22.npz lh_s11/c5b_lh_s11_results_$LAB.json lh_s11/c5b_lh_s11_perdraw_${LAB}_m22.npz; do
    [ -e "$E/$f" ] || { log "ERROR: missing output $E/$f"; exit 1; }
  done
done
log "all steps ok"
exit 0
'''
Path("run_c5b_eval.sh").write_text(ev) if not Path("run_c5b_eval.sh").exists() else (_ for _ in ()).throw(SystemExit("refusing to overwrite run_c5b_eval.sh"))
shas_m = " ".join("[%s]=%s" % (f, sha(f)) for f in ["c5b_infer.py", "c5b_memtrial.py"])
mt = HEAD % dict(what="memory trial", dirname="c5b_memtrial", shas=shas_m, s8=S8)
mt = mt.replace('[ "$(sha256sum "$SD/s8_bounds.py" | cut -c1-20)" = "%s" ] || { echo "s8_bounds.py sha mismatch"; exit 1; }\n' % S8, "").replace('[ "$(ls "$CK"/*.pth 2>/dev/null | wc -l)" -eq 10 ] || { echo "expected exactly 10 C5b checkpoints in $CK"; exit 1; }\n', '[ -e "$CK/c5b-f_s41.pth" ] || { echo "checkpoint missing"; exit 1; }\n')
mt += FIN + OTHER % dict(self="run_c5b_memtrial.sh") + WATCH + '''cd "$SD" || exit 1
# memory trial: first 5,000 BBBC-ALL rows, ONE checkpoint (c5b-f_s41), extrapolated to the 85,486 rows; the watchdog Pss and GPU peaks of the step are logged as well
run_step memtrial "$PY" c5b_memtrial.py "$CK/c5b-f_s41.pth" 5000 || exit $RC_LAST
log "memtrial output: $(tr '\\n' '|' < "$OUTDIR/memtrial.stdout.log" | cut -c1-1500)"
exit 0
'''
Path("run_c5b_memtrial.sh").write_text(mt) if not Path("run_c5b_memtrial.sh").exists() else (_ for _ in ()).throw(SystemExit("refusing to overwrite run_c5b_memtrial.sh"))

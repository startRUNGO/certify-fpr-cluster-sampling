# P5a: preregistered / frozen thresholds quoted in the SI, as one small JSON for add-only key registration (replaces make_candA_strata.py 857cabbe, never run).  WRITTEN, NOT RUN, until 实验管理 reviews it.  No data read.
# Each constant is written here and checked against the exact source line of the pinned frozen reader (the line must contain every listed fragment verbatim; exit 3 otherwise, nothing written):
#   candA|strata_cut|0 = 0.5, candA|strata_cut|1 = 0.8   c1_A_read.py            b9141c49  line 72: '("tpr<=0.5", -1, 0.5)', '("0.5<tpr<=0.8", 0.5, 0.8)', '("tpr>0.8", 0.8, 2)'
#   r1prime|rho_bar_max = -0.5, r1prime|min_strata_negative = 3   c1_r1r3_v1_9e2b45a2.py (the verdict reader of r1r3v1, script_sha256 9e2b45a2)  line 67: 'rho_bar <= -0.5 and n_neg >= 3'
#   r3prime|p_wrong_min = 0.6, r3prime|diff_min = 0.2              same file  line 72: 'main["p_wrong"] >= 0.6 and main["diff"] >= 0.2'
#   c1P2|p_max = 0.05, c1P2|other_over_tested_max = 0.5            c1_read.py (card 1 frozen reader, a2241f08)  line 85: 'obs_o < 0.5 * obs_t', 'pv <= 0.05 and obs_o < 0.5 * obs_t'
#   v1.1 (实验管理, section 4 TBDs):
#   m51|eta = 0.5, m51|pi_cap = 0.5     c1card/c1_m51_run.py (v2.2 a787897f, the script of the M5-1 full run FULL3, RUNLOG 2026-11-23T13:53:44Z)  line 30: 'LAMBDA = 0.5; PI_MAX = 0.5'
#   b1u|iccu_when_msb_msw_zero = 0.0    m3/k2prep/m3_replay_lib.py (ed9d1f17, pinned by c1_smoke)  line 32: 'if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0'
#   Paths are now relative to --results_src (v1 took --c1card).
# Output: <out>/p5a_consts.json (flat keys above).  Usage: python3 -B make_p5a_consts.py --results_src <Paper5_StainRef/results_src> --out <new dir>
import argparse, hashlib, json, sys
from pathlib import Path
SRC = {"c1card/c1_m51_run.py": "a787897fbc17d6c1b902c7c9311098c6970c5ff0933a04e29c913800f36fc1b8", "m3/k2prep/m3_replay_lib.py": "ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de",
       "c1card/c1_A_read.py": "b9141c49a51e7e725d06d1523fd4d7cf2faee15831be8be12ed0170fb0ef5361",
       "c1card/c1_r1r3_v1_9e2b45a2.py": "9e2b45a29f63bdfcfb930b4db5495e6445ec6f2ddc606532bd65b5134b0498f0",
       "c1card/c1_read.py": "a2241f08d5211600456ebfb09b88496d26dc534d176dbfa2a656f939a3510d3e"}
CHECKS = [("c1card/c1_m51_run.py", 30, ("LAMBDA = 0.5; PI_MAX = 0.5",)), ("m3/k2prep/m3_replay_lib.py", 32, ("if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0",)),
          ("c1card/c1_A_read.py", 72, ('("tpr<=0.5", -1, 0.5)', '("0.5<tpr<=0.8", 0.5, 0.8)', '("tpr>0.8", 0.8, 2)')),
          ("c1card/c1_r1r3_v1_9e2b45a2.py", 67, ("rho_bar <= -0.5 and n_neg >= 3",)),
          ("c1card/c1_r1r3_v1_9e2b45a2.py", 72, ('main["p_wrong"] >= 0.6 and main["diff"] >= 0.2',)),
          ("c1card/c1_read.py", 85, ("obs_o < 0.5 * obs_t", "pv <= 0.05 and obs_o < 0.5 * obs_t"))]
CONSTS = {"candA|strata_cut|0": 0.5, "candA|strata_cut|1": 0.8, "r1prime|rho_bar_max": -0.5, "r1prime|min_strata_negative": 3,
          "r3prime|p_wrong_min": 0.6, "r3prime|diff_min": 0.2, "c1P2|p_max": 0.05, "c1P2|other_over_tested_max": 0.5,
          "m51|eta": 0.5, "m51|pi_cap": 0.5, "b1u|iccu_when_msb_msw_zero": 0.0}


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--results_src", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(argv); D, out = Path(a.results_src), Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    for f, h in SRC.items():
        if hashlib.sha256((D / f).read_bytes()).hexdigest() != h: print("GATE (exit 3): sha differs", f, file=sys.stderr); return 3
    for f, ln, frags in CHECKS:
        L = (D / f).read_text(encoding="utf-8").split("\n")
        if len(L) < ln or any(fr not in L[ln - 1] for fr in frags): print("GATE (exit 3): line check failed", f, ln, file=sys.stderr); return 3
    out.mkdir(parents=True)
    with open(out / "p5a_consts.json", "x") as fh: json.dump(CONSTS, fh, indent=1)
    print(CONSTS); return 0


if __name__ == "__main__": sys.exit(main())

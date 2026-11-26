# P5a (CertBound): the preregistered M5-1 thresholds for M2 and M3 as a small JSON, so that the manuscript can quote them through results_master keys.  WRITTEN, NOT RUN, until 实验管理 reviews it.  No data read.
# Source of the values: the countersigned prereg 预注册_M5-1_不等概率视野抽样_v3_会签定稿_20261122T213704Z.md (sha 139c75b0, pinned), lines 67 (M2) and 68 (M3);
#   cross-check: the frozen reader results_src/c1card/c1_m51_read.py (sha fe426a87, pinned) uses the same constants (line 54 'a2 / a0 <= 0.80', line 65 'd.mean() >= 0.01', lines 55/66 '>= 2').
# The values are written as constants below and the script only checks (exit 3 on any failure, nothing written) that each appears verbatim in both pinned files:
#   prereg line 67 contains 'A2 的半宽均值 ≤ 0.80 × A0 的半宽均值' and '在至少 2 个开发站上成立';  prereg line 68 contains '≥ +0.01' and '在至少 2 个开发站上成立';
#   reader contains 'a2 / a0 <= 0.80', 'd.mean() >= 0.01', and twice '>= 2'.
# Output (new directory, 'x' mode): m51_gates.json = {"M2|hw_ratio_max": 0.80, "M2|min_sites": 2, "M3|sens_diff_min": 0.01, "M3|min_sites": 2, "n_dev_sites": 3}
# Usage: python3 -B make_m51_gates.py --prereg <...v3_会签定稿_20261122T213704Z.md> --reader <results_src/c1card/c1_m51_read.py> --out <new dir>
import argparse, hashlib, json, sys
from pathlib import Path
PREREG_SHA = "139c75b02fd6e6c18bf0a9d6cad15aa7bd2b784accabf418e212aa5b23f59c7a"
READER_SHA = "fe426a8709940b0ab72b10c4945933befa003f5ad4a6b5397ea3351b62d35c18"
GATES = {"M2|hw_ratio_max": 0.80, "M2|min_sites": 2, "M3|sens_diff_min": 0.01, "M3|min_sites": 2, "n_dev_sites": 3}


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(); [ap.add_argument("--" + k, required=True) for k in ("prereg", "reader", "out")]; a = ap.parse_args(argv); out = Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    if sha(a.prereg) != PREREG_SHA or sha(a.reader) != READER_SHA: print("GATE (exit 3): pinned sha differs", file=sys.stderr); return 3
    L = open(a.prereg, encoding="utf-8").read().split("\n"); r = open(a.reader, encoding="utf-8").read()
    checks = [("prereg l67 M2", L[66].startswith("- **M2 省标注**") and "A2 的半宽均值 ≤ 0.80 × A0 的半宽均值" in L[66] and "在至少 2 个开发站上成立" in L[66]),
              ("prereg l68 M3", L[67].startswith("- **M3 灵敏度**") and "≥ +0.01" in L[67] and "在至少 2 个开发站上成立" in L[67]),
              ("reader M2", "a2 / a0 <= 0.80" in r), ("reader M3", "d.mean() >= 0.01" in r), ("reader sites", r.count(">= 2))") >= 2), ("reader 3 sites", "at >= 2 of 3 sites" in r)]
    bad = [n for n, ok in checks if not ok]
    if bad: print("GATE (exit 3): text check failed", bad, file=sys.stderr); return 3
    out.mkdir(parents=True)
    with open(out / "m51_gates.json", "x") as fh: json.dump(GATES, fh, indent=1)
    print(GATES); return 0


if __name__ == "__main__": sys.exit(main())

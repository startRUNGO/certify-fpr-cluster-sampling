"""Finite-site (hypergeometric) minimum numbers of fields for Proposition prop:minm at level delta/K, pi = 0, r = 1 (q = 1 - eps).  WRITTEN FOR REVIEW, NOT RUN.  Pure arithmetic, no data, no rule.
Definition (checked by 实验管理 with finiteF_minm_numbers.py 3a83383a...): a site of F fields puts K0 = ceil(F q) - 1 fields (the largest integer strictly below F q) on the 'y' value, so the null holds strictly; m fields drawn without replacement are all y-fields with probability
prod_{i<m} (K0 - i)/(F - i); certification needs this probability <= alpha; the smallest such m (m = K0 + 1 if none).  Cross-check inside the script: the i.i.d. value ceil(ln(1/alpha)/ln(1/q)) must equal the registered prop_minm.json 'n_min|delta_over_K' (33).
Usage: python3 -B make_finiteF_minm.py --out <new dir>   (writes <dir>/finiteF_minm.json with flat keys 'n_min_hyp|delta_over_K|F<F>' and 'n_min_iid|delta_over_K')"""
import argparse, json, math, sys
from pathlib import Path
EPS, DELTA, K = 0.10, 0.10, 3          # protocol/eps_main, protocol/delta, protocol/candidates_K (same constants as make_prop_minm.py)
FS = (100, 300, 1000)


def m_iid(q, alpha): return math.ceil(math.log(1 / alpha) / math.log(1 / q))


def m_hyp(F, q, alpha):
    K0 = math.ceil(F * q) - 1
    m = 1
    while m <= K0:
        if sum(math.log((K0 - i) / (F - i)) for i in range(m)) <= math.log(alpha): return m
        m += 1
    return m


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args(argv); out = Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    q = 1 - EPS; alpha = DELTA / K
    reg = json.load(open(Path(__file__).resolve().parent / "prop_minm.json"))["n_min|delta_over_K"]
    iid = m_iid(q, alpha)
    if iid != reg: print("STOP (exit 3): i.i.d. value", iid, "differs from the registered prop_minm.json n_min|delta_over_K", reg, file=sys.stderr); return 3
    res = {"n_min_iid|delta_over_K": iid}
    for F in FS: res[f"n_min_hyp|delta_over_K|F{F}"] = m_hyp(F, q, alpha)
    if not (res["n_min_hyp|delta_over_K|F100"] <= res["n_min_hyp|delta_over_K|F300"] <= res["n_min_hyp|delta_over_K|F1000"] <= iid): print("STOP (exit 3): finite-F values are not monotone in F", res, file=sys.stderr); return 3
    out.mkdir(parents=True)
    with open(out / "finiteF_minm.json", "x") as fh: json.dump(res, fh, indent=1)
    print(res); return 0


if __name__ == "__main__": sys.exit(main())

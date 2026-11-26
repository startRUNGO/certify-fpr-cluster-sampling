# P5a (CertBound): numbers of Proposition minm+ (positives-only finite-sample impossibility), draft v3 756296df, section 3 table + finite-site remark.  WRITTEN, NOT RUN, until 实验管理 reviews it.  Pure arithmetic, no data.
# i.i.d.: m+ = ceil( ln(1/alpha) / ln(1 + d/(eps c_max)) ),  i.e. q+ = eps c_max / (eps c_max + d)  (q+ = probability of the w-triple in the two-point law, w = -d/c_max).
# Finite site (remark): K = largest integer strictly below F q+ (so TPR < eps holds strictly); m fields without replacement are all w-fields with probability prod_{i<m} (K-i)/(F-i);
#   necessary m = smallest m with that probability <= alpha (m = K+1 if none).  Same definition as results_src/minm/make_finiteF_minm.py (4d794251), there with q = 1 - eps.
# Grid (as v3): eps = 0.10; c_max in {90, 60}; alpha in {0.05, 0.10}; d in {0.9, 0.8} for i.i.d.; finite F in {100, 300, 1000} for d = 0.9 only.
# Gate: every computed value must equal the value printed in v3 (EXPECTED below, copied from v3 section 3); any mismatch -> exit 3, nothing written.
# Output (new directory, 'x' mode): minm_plus.json with flat keys  "iid|cmax<c>|alpha<a>|d<d>|{ratio, m_plus, m_plus_real}"  and  "hyp|cmax<c>|alpha<a>|d0.9|F<F>"  plus constants "eps", "shown|c_max" (90), "shown|d" (0.9), "shown|alpha" ([0.05, 0.10]; v1.1, 实验管理: the minm+ alpha gets its own keys, not zplus_F1meta/alphas).
# Usage: python3 -B make_minm_plus.py --out <new dir>   (v1.1 output dir: inputs_20261124/minmplus_v11; v1 e2ea6c3e kept, its output renamed NOT_USED_v1_minm_plus.json)
import argparse, json, math, sys
from pathlib import Path
EPS = 0.10; CMAX = (90, 60); ALPHAS = (0.05, 0.10); DS = (0.9, 0.8); FS = (100, 300, 1000)
EXPECTED_IID = {(90, 0.05, 0.9): 32, (90, 0.10, 0.9): 25, (60, 0.05, 0.9): 22, (60, 0.10, 0.9): 17,
                (90, 0.05, 0.8): 36, (90, 0.10, 0.8): 28, (60, 0.05, 0.8): 24, (60, 0.10, 0.8): 19}
EXPECTED_HYP = {(90, 0.05): (25, 30, 31), (90, 0.10): (20, 23, 24), (60, 0.05): (18, 21, 22), (60, 0.10): (15, 16, 17)}


def m_iid_real(cmax, alpha, d): return math.log(1 / alpha) / math.log(1 + d / (EPS * cmax))


def m_hyp(F, q, alpha):
    K = math.ceil(F * q) - 1
    m = 1
    while m <= K:
        if sum(math.log((K - i) / (F - i)) for i in range(m)) <= math.log(alpha): return m
        m += 1
    return m


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args(argv); out = Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    res = {"eps": EPS, "shown|c_max": 90, "shown|d": 0.9, "shown|alpha": list(ALPHAS)}; bad = []   # the case the main text quotes (v3 table row 1-2 and the finite-site remark)
    for c in CMAX:
        for al in ALPHAS:
            for d in DS:
                r = m_iid_real(c, al, d); m = math.ceil(r); tag = f"iid|cmax{c}|alpha{al:g}|d{d:g}"
                res[tag + "|ratio"] = d / (EPS * c); res[tag + "|m_plus_real"] = r; res[tag + "|m_plus"] = m
                if m != EXPECTED_IID[(c, al, d)]: bad.append((tag, m, EXPECTED_IID[(c, al, d)]))
            q = EPS * c / (EPS * c + 0.9)
            for F, e in zip(FS, EXPECTED_HYP[(c, al)]):
                tag = f"hyp|cmax{c}|alpha{al:g}|d0.9|F{F}"; res[tag] = m_hyp(F, q, al)
                if res[tag] != e: bad.append((tag, res[tag], e))
                if res[tag] > res[f"iid|cmax{c}|alpha{al:g}|d0.9|m_plus"]: bad.append((tag, "finite > iid"))
    if bad: print("GATE (exit 3): differs from v3", bad, file=sys.stderr); return 3
    out.mkdir(parents=True)
    with open(out / "minm_plus.json", "x") as fh: json.dump(res, fh, indent=1, allow_nan=False)
    print(json.dumps(res, indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

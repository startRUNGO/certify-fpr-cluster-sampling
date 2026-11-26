"""Minimum numbers of fields for Proposition prop:minm (two-point bound, pi = 0). Writes prop_minm.json.
Run from Paper5_StainRef/: python3 results_src/minm/make_prop_minm.py"""
import json, math
from pathlib import Path
SRC = Path(__file__).resolve().parents[1]
EPS, DELTA, K = 0.10, 0.10, 3          # protocol/eps_main, protocol/delta, protocol/candidates_K (PROTOCOL_FREEZE_20261022.md)
r = {"bbbc041_all": json.load(open(SRC/"minm/bbbc_rbar_vmread.json"))["value"],
     "mpidb_full": json.load(open(SRC/"hc/mpidb_hc_b2/s7/mh_s7_results.json"))["analytic_fields"]["a2|mh"]["ni_over_nmax_mean"],
     "lacuna_B": json.load(open(SRC/"hc/lacuna_hc_paired_B/s7/lh_s7_results.json"))["analytic_fields"]["a2|lh"]["ni_over_nmax_mean"]}
mstar = lambda a, rr: math.ceil(math.log(1/a)/math.log(1+EPS*rr/(1-EPS)))
assert mstar(DELTA, 1.0) == 22   # = s5_results/analytic|0.1|n_min, already registered
out = {"n_min|delta_over_K": mstar(DELTA/K, 1.0), "r_mean|bbbc041_all": r["bbbc041_all"]}   # MP-IDB / Lacuna-B r are already registered (hc s7 analytic_fields)
for s, rr in r.items():
    out[f"m_star|delta|{s}"] = mstar(DELTA, rr)
    out[f"m_star|delta_over_K|{s}"] = mstar(DELTA/K, rr)
json.dump(out, open(SRC/"minm/prop_minm.json", "w"), indent=1)
print(out)

"""Pure arithmetic, no data read (review F1b, 2026-10-28): does B1u, as S3 calls it (level delta/K = 0.10/3),
certify FPR <= 0.10 from m = 22 identical fields with zero exceedances?  The function b1u below is copied verbatim
from results_src/mh/mh_rev1.py (the REV-1 code behind rev1_results); the script asserts the copy is identical.
Run from Paper5_StainRef/: python3 -I results_src/minm/b1u_identical_fields.py"""
import json, math, re
from pathlib import Path
import numpy as np
from scipy.stats import beta as betad, f as fdist
SRC = Path(__file__).resolve().parents[1]
def b1u(x, n, delta):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X/N
    n0 = (N-(n**2).sum()/N)/(m-1); pi = x/n; msb = (n*(pi-ph)**2).sum()/(m-1); msw = (n*pi*(1-pi)).sum()/(N-m)
    if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0
    else:
        F = msb/msw; flo = fdist.ppf(delta, m-1, N-m); th = (F/flo-1)/n0; iccu = 0.0 if th <= 0 else th/(1+th)
    nbw = (n**2).sum()/N; neff = N/(1+(nbw-1)*iccu); xeff = X*neff/N
    return float(betad.ppf(1-delta, xeff+1, neff-xeff))
# --- the copy must equal the source function line for line
src = (SRC/"mh/mh_rev1.py").read_text().splitlines(); me = Path(__file__).read_text().splitlines()
grab = lambda L: L[next(i for i, l in enumerate(L) if l.startswith("def b1u(")):][:8]
assert grab(src) == grab(me), "b1u copy differs from mh_rev1.py"
M, EPS, ALPHA = 22, 0.10, 0.10/3          # protocol/m, protocol/eps_main, delta/K as called in mh_rev1.py (b1u(x, ni, 0.10/3))
def case(n_per_field):
    n = [n_per_field]*M; x = [0]*M
    ub = b1u(x, n, ALPHA)
    return {"n_per_field": n_per_field, "N": n_per_field*M, "UB": ub, "certifies": ub <= EPS}
# smallest field size at which 22 identical zero-exceedance fields are certified
first = next(k for k in range(1, 1000) if case(k)["certifies"])
out = {"note": "x_i = 0 for all fields gives pi = 0, msw = 0 and msb = 0, so the code's branch sets ICC_U = 0 (no F interval is formed); n_eff = N, x_eff = 0, UB = Beta^{-1}(1-alpha; 1, N) = 1 - alpha^(1/N).",
       "alpha": ALPHA, "eps": EPS, "m": M, "smallest_n_per_field_certified": first,
       "closed_form_check_n_ge_2": [abs(case(k)["UB"] - (1 - ALPHA**(1/(k*M)))) < 1e-12 for k in (2, 5, 50)],
       "n_per_field_1": "N - m = 0, so msw = 0/0 = NaN and the code returns UB = NaN (no certification); reported as computed",
       "grid_n_per_field": [case(k) for k in (1, 2, 3, 5, 10, 50)]}
# site cases: n = round(r * n_max); n_max known locally for MP-IDB and Lacuna-B only (BBBC n_max not archived locally)
rm = json.load(open(SRC.parent/"paper_en/results_master.json"))["values"]
sites = {"mpidb_full": ("mpidb_hc_b2__mh_s7_results/analytic_fields|a2|mh|nmax", "mpidb_hc_b2__mh_s7_results/analytic_fields|a2|mh|ni_over_nmax_mean"),
         "lacuna_B": ("lacuna_hc_paired_B__lh_s7_results/analytic_fields|a2|lh|nmax", "lacuna_hc_paired_B__lh_s7_results/analytic_fields|a2|lh|ni_over_nmax_mean")}
out["sites"] = {}
for s, (kn, kr) in sites.items():
    nmax, r = rm[kn], rm[kr]
    out["sites"][s] = {"nmax": nmax, "r_mean": r, "r1": case(int(nmax)), "r_mean_case": case(int(round(r*nmax)))}
nb = rm["rev1_results/agg|a2|bbbc041_all|22|S3|0.1|n_labeled_neg"]
out["bbbc041_all"] = {"r_mean": rm["prop_minm/r_mean|bbbc041_all"], "nmax": "not archived locally", "mean_labeled_negatives_per_field_at_m22": nb/M, "case_at_that_mean": case(int(round(nb/M)))}
print(json.dumps(out, indent=1))

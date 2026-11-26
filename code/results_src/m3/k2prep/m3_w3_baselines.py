# P5 W3 (gap table draft2 16bcfbb6): three commonly used cluster-aware upper bounds as BASELINES for the certification rules, evaluated on the same replayed REV-1
# draws of the three development sites.  LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 approves.
# NATURE: comparison baselines for the paper and controls for K2; no selection, no judgement line, not confirmatory evidence.
# Each baseline is a single-candidate rule (candidates c0 and c1s), level delta = 0.10, run with the REV-1 fixed-sequence test on the same descending grid as S3/B1u
# (03_theory remark: a bound valid at each fixed grid point gives a valid fixed-sequence procedure).  Per grid point the field-level data are x_i (negatives above the
# threshold) and n_i (negatives) for the m = 22 drawn fields; p_hat = X/N.
#   CB  cluster bootstrap percentile: B = 2000 resamples of the 22 fields with replacement; ratio sum(w x)/sum(w n); UB = (1 - delta) quantile.  The resampling weights
#       are drawn once per (site, seed, draw, eps, candidate) from default_rng([20261031, site index, seed, draw, eps code, candidate code]) (common random numbers
#       along the grid); this stream is separate from, and does not consume, the REV-1 stream.
#   GEE cluster-robust (sandwich) Wald bound on the ratio: V = m/(m-1) * sum((x_i - p_hat n_i)^2) / N^2; UB = p_hat + t_{m-1, 1-delta} sqrt(V), clipped to [0, 1].
#       Plain form as used in practice; when X = 0 it gives UB = 0 (known anti-conservatism at zero counts; kept, not patched, because the baseline must be the
#       common practice).
#   BB  beta-binomial profile-likelihood bound: x_i ~ BetaBin(n_i, mu, rho); profile log-likelihood l_p(mu) = max_rho l(mu, rho), rho in [1e-6, 0.999];
#       UB = the mu > mu_hat with 2 (l_p(mu_hat) - l_p(mu)) = chi2_1(1 - 2 delta) (one-sided level delta); UB = 1 if not reached.  Asymptotic, model-based.
# Also recorded from the same replay for reference: S3, B1u_c0, B1u_c1s.  The replay self-check (m3_replay_lib, 1e-12 vs stored REV-1) must pass before anything is
#   written.  Expected runtime: the BB profile is the slowest part (bounds are evaluated only along the fixed-sequence path, about a dozen grid points per draw);
#   rough estimate on the host: 1-2 hours for 3 sites x 5 seeds x 200 draws x 2 eps x 2 candidates.
# Exceptions (实验管理 2026-10-30): the REV-1 fixed-sequence test turns an exception in the bound into "not certified" (j = 0).  That is kept (the library is not
#   changed), but every baseline bound is wrapped so that each exception is COUNTED by (site, eps, rule) and type, then re-raised.  Per cell the output reports
#   n_bound_exceptions and the three most frequent types; a cell with exceptions > 1% of its draws is flagged "exceptions>1%" and that baseline is reference only.
#   S3 and B1u are not wrapped (their behaviour must equal REV-1).
# Output per rule, site, eps: n, coverage, fpr_mean, fpr_median, cannot_certify, sens, ds_eps, ds_own, n_bound_exceptions, exception_types, flag.  One new JSON in a new directory ("x", allow_nan=False).
# Usage: python -I -B m3_w3_baselines.py --a2link .. --ids_iter1 .. --mhlh .. --rev1_bbbc .. --rev1_mh .. --rev1_lh .. --mh_rev1_py .. --out <dir/json>
import argparse, json, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit, logit
from scipy.stats import betabinom, chi2, t as tdist

sys.path.insert(0, str(Path(__file__).resolve().parent)); import m3_replay_lib as L   # noqa: E402

DELTA = 0.10; B = 2000; BOOT_SEED = 20261031; CAND_CODE = {"c0": 0, "c1s": 1}


def gee_ub(x, n, delta=DELTA):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); p = x.sum()/N
    V = m/(m-1)*((x - p*n)**2).sum()/N**2
    return float(min(1.0, max(0.0, p + tdist.ppf(1-delta, m-1)*np.sqrt(V))))


def cb_ub(x, n, W, delta=DELTA):
    x = np.asarray(x, float); n = np.asarray(n, float); num = (W * x).sum(1); den = (W * n).sum(1); ok = den > 0   # not W @ x: spurious Accelerate warnings
    return float(np.quantile(num[ok]/den[ok], 1-delta))


def bb_loglik(x, n, mu, rho):
    a = mu*(1-rho)/rho; b = (1-mu)*(1-rho)/rho; return float(betabinom.logpmf(x, n, a, b).sum())


def bb_profile(x, n, mu):
    r = minimize_scalar(lambda z: -bb_loglik(x, n, mu, expit(z)), bounds=(logit(1e-6), logit(0.999)), method="bounded", options={"xatol": 1e-6})
    return -float(r.fun)


def bb_ub(x, n, delta=DELTA):
    x = np.asarray(x, float); n = np.asarray(n, float); lo, hi = 1e-9, 1-1e-9
    r = minimize_scalar(lambda z: -bb_profile(x, n, expit(z)), bounds=(logit(lo), logit(hi)), method="bounded", options={"xatol": 1e-7})
    mu_hat = float(expit(r.x)); lmax = -float(r.fun); crit = chi2.ppf(1-2*delta, 1)
    g = lambda mu: 2*(lmax - bb_profile(x, n, mu)) - crit
    if g(hi) < 0: return 1.0
    return float(brentq(g, max(mu_hat, lo), hi, xtol=1e-9))


EXC = defaultdict(Counter)                                                  # (site, eps, rule) -> Counter of exception type names


def counted(name, fn):
    def g(x, n, ctx):
        try: return fn(x, n, ctx)
        except Exception as e:
            EXC[(ctx[0], ctx[3], name)][type(e).__name__] += 1; raise
    return g


def make_rules():
    Wcache = {}

    def weights(ctx, cand):
        site, sd, d_, eps, _ = ctx; key = (site, sd, d_, eps, cand)
        if key not in Wcache:
            rng = np.random.default_rng([BOOT_SEED, L.DEV.index(site), sd, d_, int(round(eps*100)), CAND_CODE[cand]])
            Wcache.clear(); Wcache[key] = rng.multinomial(L.M, np.full(L.M, 1/L.M), size=B).astype(float)
        return Wcache[key]
    rules = {}
    for k in ("c0", "c1s"):
        rules["CB_" + k] = (k, counted("CB_" + k, lambda x, n, ctx, k=k: cb_ub(x, n, weights(ctx, k))))
        rules["GEE_" + k] = (k, counted("GEE_" + k, lambda x, n, ctx: gee_ub(x, n)))
        rules["BB_" + k] = (k, counted("BB_" + k, lambda x, n, ctx: bb_ub(x, n)))
    return rules


def summarise(R, eps, site, rule):
    fpr = np.array([r["fpr"] for r in R]); ex = EXC.get((site, eps, rule), Counter()); nex = int(sum(ex.values()))
    return dict(n=len(R), coverage=float((fpr <= eps).mean()), fpr_mean=float(fpr.mean()), fpr_median=float(np.median(fpr)),
                cannot_certify=float(1-np.mean([r["cert"] for r in R])), sens=float(np.mean([r["sens"] for r in R])),
                ds_eps=float(np.mean([r["dse"] for r in R])), ds_own=float(np.mean([r["dso"] for r in R])),
                n_bound_exceptions=nex, exception_types=[[k, v] for k, v in ex.most_common(3)],
                flag=("exceptions>1%" if nex > 0.01*len(R) else None))


def main():
    ap = argparse.ArgumentParser()
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); out = Path(a.out); assert not out.exists(), "output exists"
    inp = L.Inputs(a.a2link, a.ids_iter1, a.mhlh, a.rev1_bbbc, a.rev1_mh, a.rev1_lh, a.mh_rev1_py)
    EXC.clear(); rules = make_rules(); rows = []
    for site in L.DEV:
        for sd in L.SEEDS:
            r_, _ = L.replay(inp, sd, site, extra=rules); rows += r_
    bad = L.self_check(inp, rows)
    if bad: raise SystemExit("self-check failed, nothing written: " + repr(bad[:5]))
    names = ["S3", "B1u_c0", "B1u_c1s"] + list(rules)
    cells = [dict(site=s, eps=e, rule=nm, **summarise([r for r in rows if r["site"] == s and r["eps"] == e and r["rule"] == nm], e, s, nm))
             for s in L.DEV for e, _ in L.SETTINGS for nm in names]
    rep = dict(script_sha256=L.sha(__file__), lib_sha256=L.sha(L.__file__), inputs_sha256=dict(L.PIN), score_files_sha256=inp.used,
               nature="comparison baselines (gap table W3); no judgement, not confirmatory evidence", delta=DELTA, bootstrap=dict(B=B, seed=BOOT_SEED),
               self_check="passed (S3, B1u_c0, B1u_c1s vs stored REV-1 within 1e-12)", cells=cells)
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False)
    out.parent.mkdir(parents=False, exist_ok=False)
    with open(out, "x") as fh: fh.write(txt)
    print("written", out)


if __name__ == "__main__":
    main()

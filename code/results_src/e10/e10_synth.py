# P5 E10 controlled synthetic simulation with known truth (spec P5_规格_E10_合成受控模拟_20261103_draft1.md; reviewer A2).  LOCAL, no project data: every score is generated here.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# NATURE: a direct test of the validity claims on a population whose generating law is known, over the factors the review named: between-field shift (hence the ICC of the exceedance indicators), cluster size, positive
#   fraction and tail weight.  No rule and no constant is changed.  Real data, M5, IML, PlasmoCount, NIH, any image, any classifier: never read.
# Generator (per site realisation, F = 300 fields): negatives per field n_f ~ max(5, Poisson(r)); positives per field ~ Poisson(pi / (1 - pi) n_f); field location mu_f ~ N(0, tau^2) (shared by both classes);
#   negative score = mu_f + e, positive score = mu_f + SEP + e, e ~ N(0, 1) or Student t(5) rescaled to unit variance.  A source set (tau = 0, about 4000 negative cells = ceil(4000 / r) fields, same pi, SEP and tail) is drawn per realisation only to fit the Platt q that the frozen
#   K1-v bound needs (label-using, as the NIH validation fold is in the real runs).  r3: the frozen KV.platt (BFGS) can lose precision on a very separable source set ("platt fit" assertion); then and only then a Newton fit
#   is used (platt_synth, as the A0-ctrl run r1 did), and the fallback seeds are listed in the output (platt_fallback).
# Truth: the realised FPR / sensitivity of the synthetic site (F = 300 fields, so close to the super-population), evaluated by the frozen library population exactly as for the real sites; the realised ICC of the exceedance
#   indicators at the own-oracle threshold (eps = 0.10) is reported per realisation (one-way ANOVA, b1u formulas).
# (A) Rules S3 and SPLIT6 (and B1u_c1s): frozen m3_split_dev_eval.run_site_seed on the synthetic site as a Fake Inputs object, m = 22 fields, 200 draws per realisation, 5 realisations per cell (seeds 41..45 as realisation ids),
#   eps in {0.10, 0.05}; coverage per realisation, seed-cluster interval over the 5 realisations (ebatch_lib), cannot-certify share, sensitivity cost against the own oracle.  Library conventions.
# (B) The 2-alpha question at a FIXED threshold (no grid, no selection): tau_0 = the c0 threshold whose site FPR is the smallest value >= eps (H0 true at the boundary) and tau_half (FPR about eps / 2); for 2000 draws of 22 fields
#   per realisation the upper bound is computed and the type-I rate P(UB <= eps | FPR >= eps) and the power P(UB <= eps | FPR ~ eps/2) are reported for: B1 (point ICC, s8_bounds.b1), B1u (frozen lib b1u: F limit and
#   Clopper-Pearson limit each at level alpha) at alpha = 0.10 and at alpha = 0.10 / 3 (the per-candidate level of S3), and the DIAGNOSTIC arm B1u-split (same formulas with alpha / 2 for each of the two steps; NOT a rule of the paper).
# Design: one-factor-at-a-time around the base cell (tau 0.5, r 100, pi 0.05, SEP 2.5, normal): tau in {0, 0.25, 0.5, 1.0}; r in {20, 100, 500}; pi in {0.01, 0.05, 0.2, 0.4}; tail in {normal, t5}  (10 cells).  --full runs the
#   factorial (4 x 3 x 4 x 2 = 96 cells).  Streams: default_rng([20261103, 10, cell index, realisation id, purpose]).
# Anchors (exit 3, no output directory): library shas (replay lib, split eval, s8_bounds); generator checks run before anything is written: the realised ICC at tau = 0 is close to 0 and at tau = 1 close to the analytic value
#   Var(p(mu)) / (p (1 - p)) computed by quadrature; at tau = 0 (iid cells) the type-I rate of <EMAIL> at the boundary threshold, 2000 draws of 22 fields, is <= alpha + 3 standard errors (r2: implemented in generator_checks).
# Reading (fixed before results): E10-a: S3 and SPLIT6 coverage >= 0.90 in every cell where the realised ICC <= 0.2 and pi <= 0.2 -> validity on the synthetic grid; any cell below 0.90 is listed with its factor values (the
#   boundary of the approximate guarantee), whatever the factor.  E10-b: B1u type-I at the boundary above alpha (+ 3 SE) -> the two steps cost more than alpha and the size is reported; below alpha -> no inflation measured in
#   these cells (the union bound 2 alpha is a bound, not a measurement).  E10-c: B1 (point ICC) under-covers where the ICC is large (expected, reported).  Nothing is chosen.
# Usage: python e10_synth.py [--full] [--cells 0,3] [--reps 5] [--ndraw 200] [--fixed-draws 2000] --out <new dir>/e10.json
import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[_k] = "8"
import argparse, itertools, json, sys, time
from pathlib import Path
import numpy as np
from scipy.stats import beta as betad, f as fdist, norm
HERE = Path(__file__).resolve().parent; R = HERE.parents[0]
sys.path[:0] = [str(R / "e_batch"), str(R / "m3" / "split"), str(R / "m3" / "k2prep"), str(R / "c5b")]
import ebatch_lib as EB

PINS = {"m3_replay_lib": "ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de", "m3_split_dev_eval": "cbf1e7f4f33dd1c74cbfd0f3ac4d34a397e8a7ef44b82c9365ca117733bf6457", "s8_bounds": "54a2af8d3704c16b0994dcdc584d4f407cd9564cc07d49ad74cdf7c4cad5ae7f"}
F_SITE = 300; SEP = 2.5; BASE = dict(tau=0.5, r=100, pi=0.05, tail="normal"); LEVELS = dict(tau=(0.0, 0.25, 0.5, 1.0), r=(20, 100, 500), pi=(0.01, 0.05, 0.2, 0.4), tail=("normal", "t5")); EPSS = (0.10, 0.05); DELTA = 0.10; M = 22


def cells(full):
    if full: return [dict(tau=t, r=r, pi=p, tail=tl) for t, r, p, tl in itertools.product(*[LEVELS[k] for k in ("tau", "r", "pi", "tail")])]
    out = [dict(BASE)]
    for k in LEVELS: out += [dict(BASE, **{k: v}) for v in LEVELS[k] if v != BASE[k]]
    return out


def noise(rng, n, tail):
    if tail == "normal": return rng.normal(size=n)
    return rng.standard_t(5, size=n) / np.sqrt(5.0 / 3.0)


def make_site(rng, tau, r, pi, tail, F=F_SITE, sep=SEP):
    """Scores s, labels y, field ids g of one synthetic site."""
    nf = np.maximum(5, rng.poisson(r, F)); npos = rng.poisson(pi / (1 - pi) * nf); mu = rng.normal(0, tau, F) if tau > 0 else np.zeros(F); s, y, g = [], [], []
    for f in range(F):
        s.append(mu[f] + noise(rng, nf[f], tail)); y.append(np.zeros(nf[f], int)); g.append(np.full(nf[f], f))
        if npos[f]: s.append(mu[f] + sep + noise(rng, npos[f], tail)); y.append(np.ones(npos[f], int)); g.append(np.full(npos[f], f))
    return np.concatenate(s), np.concatenate(y), np.concatenate(g)


def anova_icc(x, n):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); ph = x.sum() / N; pi = x / n
    n0 = (N - (n ** 2).sum() / N) / (m - 1); msb = (n * (pi - ph) ** 2).sum() / (m - 1); msw = (n * pi * (1 - pi)).sum() / (N - m)
    d = msb + (n0 - 1) * msw; return 0.0 if d <= 0 else float(min(1.0, max(0.0, (msb - msw) / d)))


def realised_icc(s, y, g, eps=0.10):
    neg = y == 0; t = float(np.quantile(s[neg], 1 - eps)); fl = np.unique(g[neg]); n = np.array([(g[neg] == f).sum() for f in fl], float); x = np.array([(s[neg][g[neg] == f] >= t).sum() for f in fl], float); return anova_icc(x, n)


def analytic_icc(tau, eps=0.10):
    """Gaussian model, equal field sizes: ICC of exceedance indicators at the threshold with site FPR = eps: Var(p(mu)) / (p (1 - p))."""
    t = norm.isf(eps) * np.sqrt(1 + tau ** 2) if tau > 0 else norm.isf(eps); z = np.linspace(-8, 8, 4001); w = norm.pdf(z); w /= w.sum(); p = norm.sf(t - tau * z)
    return float(((w * p ** 2).sum() - (w * p).sum() ** 2) / (eps * (1 - eps)))


def _platt_nll(s, y):
    from scipy.special import log_expit
    return lambda p: float(-(y * log_expit(p[0] * s + p[1]) + (1 - y) * log_expit(-(p[0] * s + p[1]))).sum())


def _platt_grad(s, y, p):
    r = 1.0 / (1.0 + np.exp(-(p[0] * s + p[1]))) - y; return np.array([float((r * s).sum()), float(r.sum())])


def platt_newton(s, y, tol=1e-9, max_iter=100):
    """Same two-parameter NLL as KV.platt, Newton-Raphson from (1, 0) with analytic gradient and Hessian; the step is halved only if the NLL would increase (as results_src/m3/a0ctrl r1, platt_newton)."""
    s = np.asarray(s, float); y = np.asarray(y, float); nll = _platt_nll(s, y); p = np.array([1.0, 0.0]); f = nll(p)
    for _ in range(max_iter):
        g = _platt_grad(s, y, p)
        if np.linalg.norm(g) < tol: break
        sg = 1.0 / (1.0 + np.exp(-(p[0] * s + p[1]))); w = sg * (1 - sg)
        H = np.array([[float((w * s * s).sum()), float((w * s).sum())], [float((w * s).sum()), float(w.sum())]]); step = np.linalg.solve(H, g); lam = 1.0
        while nll(p - lam * step) > f and lam > 1e-12: lam /= 2
        p = p - lam * step; f = nll(p)
    g = _platt_grad(s, y, p); assert np.linalg.norm(g) < tol, ("platt newton did not converge", float(np.linalg.norm(g)))
    return p, float(f), float(np.linalg.norm(g))


def platt_synth(KV, s, y):
    """(r3) KV.platt first; only if the frozen fit raises its 'platt fit' assertion (BFGS precision loss on a very separable synthetic source set) the Newton fit above.  Only for the synthetic source sets; the frozen
    library is not changed.  Any other assertion propagates.  Returns ((a, b), info)."""
    s = np.asarray(s, float); y = np.asarray(y, float)
    try: a_, b_ = KV.platt(s, y); return (a_, b_), dict(fallback=False)
    except AssertionError as ex:
        if not (ex.args and isinstance(ex.args[0], tuple) and ex.args[0][0] == "platt fit"): raise
        msg = repr(ex.args)
    pn, fn, gn = platt_newton(s, y)
    from scipy.optimize import minimize
    nll = _platt_nll(s, y); rb = minimize(nll, x0=np.array([1.0, 0.0]), method="BFGS"); assert fn <= float(nll(rb.x)) + 1e-9, ("newton NLL above the BFGS end point", fn, float(nll(rb.x)))
    return (float(pn[0]), float(pn[1])), dict(fallback=True, nll=fn, grad_norm=gn, bfgs_nll_at_end=float(nll(rb.x)), frozen_assertion=msg)


def source_set(rng, r, pi, tail, n_cells=4000, sep=SEP):
    """Source set for the Platt fit: tau = 0, about n_cells NEGATIVE cells (F_src = ceil(n_cells / r) fields; positives add pi / (1 - pi) of that), same pi, separation and tail as the site."""
    return make_site(rng, 0.0, r, pi, tail, F=int(np.ceil(n_cells / r)), sep=sep)


def b1u_split(x, n, delta):
    """DIAGNOSTIC only: the b1u formulas with delta / 2 for each of the two steps (not a rule of the paper)."""
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X / N; n0 = (N - (n ** 2).sum() / N) / (m - 1); pi = x / n
    msb = (n * (pi - ph) ** 2).sum() / (m - 1); msw = (n * pi * (1 - pi)).sum() / (N - m); a = delta / 2
    if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0
    else:
        F = msb / msw; flo = fdist.ppf(a, m - 1, N - m); th = (F / flo - 1) / n0; iccu = 0.0 if th <= 0 else th / (1 + th)
    nbw = (n ** 2).sum() / N; neff = N / (1 + (nbw - 1) * iccu); xeff = X * neff / N; return float(betad.ppf(1 - a, xeff + 1, neff - xeff))


def fixed_threshold_test(L, B8, s, y, g, eps, rng, ndraw, kappa_levels=(1.0, 0.5)):
    """Type-I (kappa = 1: smallest site FPR >= eps) and power (kappa = 0.5) of the upper bounds at one FIXED threshold; fields drawn as in the rules (22 of the fields with negatives)."""
    neg = y == 0; sneg = np.sort(s[neg]); N = len(sneg); fl = np.unique(g[neg]); negf = [s[neg][g[neg] == f] for f in fl]; out = {}
    u = np.unique(sneg); fp = 1 - np.searchsorted(sneg, u, side="left") / N
    for kappa in kappa_levels:
        target = kappa * eps; ok = np.where(fp >= target)[0]; j = ok[-1]; tau0 = float(u[j]); site_fpr = float(fp[j])           # smallest FPR still >= target
        arms = {"B1": [], "<EMAIL>": [], "<EMAIL>/3": [], "<EMAIL>": []}
        for d in range(ndraw):
            F = rng.choice(len(fl), M, replace=False); n = np.array([len(negf[i]) for i in F], float); x = np.array([(negf[i] >= tau0).sum() for i in F], float)
            arms["B1"].append(B8.b1(x, n, DELTA)[0] <= eps); arms["<EMAIL>"].append(L.b1u(x, n, DELTA) <= eps); arms["<EMAIL>/3"].append(L.b1u(x, n, DELTA / 3) <= eps); arms["<EMAIL>"].append(b1u_split(x, n, DELTA) <= eps)
        out[str(kappa)] = dict(site_fpr=site_fpr, certify_rate={k: float(np.mean(v)) for k, v in arms.items()}, n_draws=ndraw)
    return out


def generator_checks(rng, L=None, B8=None, ndraw=2000):
    """Anchors on the generator, run before any output is written.  With L and B8 (the frozen libraries) the iid type-I anchor of the spec is included: at tau = 0 the cells are exchangeable, so the type-I rate of
    <EMAIL> at the boundary threshold (site FPR = the smallest value >= eps) over `ndraw` draws of 22 fields must not exceed alpha + 3 standard errors, else exit 3."""
    res = {}
    for tau in (0.0, 1.0):
        s, y, g = make_site(rng, tau, 100, 0.05, "normal"); res[f"icc_tau{tau}"] = (realised_icc(s, y, g), analytic_icc(tau))
    if res["icc_tau0.0"][0] > 0.03: raise EB.AnchorError(("realised ICC at tau = 0 not near 0", res["icc_tau0.0"]))
    if abs(res["icc_tau1.0"][0] - res["icc_tau1.0"][1]) > 0.06: raise EB.AnchorError(("realised ICC at tau = 1 differs from the analytic value", res["icc_tau1.0"]))
    if L is not None:
        s, y, g = make_site(rng, 0.0, 100, 0.05, "normal"); r = fixed_threshold_test(L, B8, s, y, g, 0.10, rng, ndraw, kappa_levels=(1.0,))["1.0"]
        rate = r["certify_rate"]["<EMAIL>"]; tol = DELTA + 3 * float(np.sqrt(DELTA * (1 - DELTA) / ndraw)); res["iid_type1_B1u_0.10"] = dict(rate=rate, limit=tol, site_fpr=r["site_fpr"], n_draws=ndraw)
        if rate > tol or r["site_fpr"] < 0.10: raise EB.AnchorError(("iid type-I rate of <EMAIL> at tau = 0 above alpha + 3 SE (or boundary FPR below eps)", res["iid_type1_B1u_0.10"]))
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--full", action="store_true"); ap.add_argument("--cells", default=None); ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--ndraw", type=int, default=200); ap.add_argument("--fixed-draws", type=int, default=2000); a = ap.parse_args(); out = Path(a.out); assert not out.parent.exists(), "output directory exists: use a new directory"
    import m3_replay_lib as L, m3_split_dev_eval as SP, m3_k1v_precheck as KV, s8_bounds as B8
    for name, mod in (("m3_replay_lib", L), ("m3_split_dev_eval", SP), ("s8_bounds", B8)):
        if EB.sha(mod.__file__) != PINS[name]: raise EB.AnchorError((name, "differs from its pin"))
    assert SP.L is L; L.M = M; L.NDRAW = a.ndraw; site = "bbbc041_all"; assert site in L.SITECFG and site in L.DEV
    platt_info = {}; chk = generator_checks(np.random.default_rng([20261103, 10, 9999, 0, 0]), L, B8, ndraw=2000); grid = cells(a.full); sel = [int(x) for x in a.cells.split(",")] if a.cells else range(len(grid)); res_cells = {}; t0 = time.time()
    class Synth:
        used = {}
        def __init__(self, sites): self.sites = sites
        def load_site(self, sd, st): return self.sites[sd]
    for ci in sel:
        cell = grid[ci]; sites, platt, icc, rows_by = {}, {}, {}, {}; fixed = {e: {} for e in EPSS}
        for rid in range(a.reps):
            sd = 41 + rid; s, y, g = make_site(np.random.default_rng([20261103, 10, ci, sd, 0]), **cell); sites[sd] = (s, y, g); icc[sd] = realised_icc(s, y, g)
            ss, ys, _ = source_set(np.random.default_rng([20261103, 10, ci, sd, 1]), cell["r"], cell["pi"], cell["tail"]); platt[sd], pinfo = platt_synth(KV, ss, ys); platt_info[(ci, sd)] = pinfo
        SP.E.EXC.clear(); SP.E.DEG.clear(); SP.E.DROP.clear(); SP.CHOSEN.clear(); SP.C2STOP.clear(); inp = Synth(sites)
        for sd in sites: rows_by[sd] = SP.run_site_seed(inp, sd, site, platt[sd])
        summ = {}
        for eps in EPSS:
            for rule in ("S3", "SPLIT6", "B1u_c1s"):
                cov = [float(np.mean([r["fpr"] <= eps for r in rows_by[sd] if r["rule"] == rule and abs(r["eps"] - eps) < 1e-12])) for sd in sites]
                cc = [float(1 - np.mean([r["cert"] for r in rows_by[sd] if r["rule"] == rule and abs(r["eps"] - eps) < 1e-12])) for sd in sites]
                dso = [float(np.mean([r["dso"] for r in rows_by[sd] if r["rule"] == rule and abs(r["eps"] - eps) < 1e-12])) for sd in sites]
                summ[f"{rule}|{eps}"] = dict(coverage=float(np.mean(cov)), coverage_per_realisation=cov, coverage_ci90_seedcluster=EB.seed_cluster_ci(cov), cannot_certify=float(np.mean(cc)), ds_own=float(np.mean(dso)))
            for sd in sites: fixed[eps][str(sd)] = fixed_threshold_test(L, B8, *sites[sd], eps, np.random.default_rng([20261103, 10, ci, sd, 2]), a.fixed_draws)
        res_cells[str(ci)] = dict(factors=cell, realised_icc_per_realisation=[icc[sd] for sd in sites], realised_icc_mean=float(np.mean(list(icc.values()))), library_rules=summ, fixed_threshold=fixed)
    flags = [dict(cell=int(ci), factors=c["factors"], rule=k, coverage=v["coverage"]) for ci, c in res_cells.items() for k, v in c["library_rules"].items() if k.split("|")[0] in ("S3", "SPLIT6") and v["coverage"] < 0.90]
    res = dict(version="e10_r1", spec="P5_规格_E10_合成受控模拟_20261103_draft1", reps=a.reps, ndraw=a.ndraw, fixed_draws=a.fixed_draws, m=M, f_site=F_SITE, sep=SEP, cells=res_cells, generator_checks=chk, reading_E10_a_cells_below_090=flags,
               platt_fallback=[[ci, sd] for (ci, sd), v in platt_info.items() if v["fallback"]], platt_info={f"{ci}|{sd}": v for (ci, sd), v in platt_info.items() if v["fallback"]}, pins=PINS, seconds=time.time() - t0, ci_definition=EB.CI_DEFINITION, note="synthetic; no project data; library conventions for S3 / SPLIT6")
    EB.write_out(out, res)


if __name__ == "__main__":
    try: main()
    except EB.AnchorError as e: print("ANCHOR MISMATCH (exit 3, no output):", e, file=sys.stderr); sys.exit(3)

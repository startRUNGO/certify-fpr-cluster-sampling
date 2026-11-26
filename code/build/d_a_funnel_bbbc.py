# P5a D-a (实验管理 2026-11-25): READ-ONLY per-field export for the funnel figure F-a, development site BBBC-ALL only.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Nature: descriptive diagnostic.  The threshold tau_eps is the site-level threshold at which the pooled FPR of ALL BBBC-ALL negatives is eps = 0.10; it is not the
#   threshold of any rule (the figure caption says so).  No rule is run, nothing is trained, no confirmation site is read.
# Data: arm A2, seed 41 (as Figure 2, figmech_v39), through m3_replay_lib (ed9d1f17): Inputs(...).load_site / population(...), uncentered scores (candidate c0, shift 0).
#   Per field f of the site (every field with at least one negative, the library population): n_f negatives and X_f = #{negatives with score >= tau_eps}.
#   tau_eps: with N negatives sorted in descending order and k = round(eps * N), tau_eps is the k-th largest negative score; the realized pooled FPR
#   fpr_site = #{negatives >= tau_eps} / N is exported (ties make it >= k / N).
#   Funnel: binomial 5 % and 95 % quantiles of X / n under p = fpr_site for every n from 1 to max n_f (two-sided 90 % limits for exchangeable cells);
#   a field is outside when X_f < q05(n_f) or X_f > q95(n_f); the expected outside share under the binomial model (sum over fields of P(X < q05) + P(X > q95), / F) is exported
#   next to the observed share.  Dispersion: Pearson X^2 = sum (X_f - n_f p)^2 / (n_f p (1 - p)) and X^2 / (F - 1); one-way ANOVA ICC point estimate (as in b1u, without the
#   upper confidence limit) and the design effect 1 + (n_bar_w - 1) ICC with n_bar_w = sum n_f^2 / N.
# Anchors (exit 3, nothing written): library pin; figmech_v39_numbers.json pin (a4b23923); N == figmech n_neg_bbbc_all; the pooled FPR at the figmech source threshold
#   (tau_src_logit) recomputed here == figmech fpr_at_tau_src_bbbc_all EXACTLY; n_f of the 22 fields of figmech draw 0 (its 'fields' list) == its 'n_f' list.
# Guard: no argument may contain plasmo / iml / r8 / m5 (case-insensitive), as e12_stress; the library refuses confirmation-site paths too.
# Output: NEW directory (argument --out, must not exist), files in 'x' mode: d_a_funnel_bbbc.json, SHA256SUMS.  allow_nan=False.
# Run (VM, pyproject env, nice 19, one thread): python -I -B d_a_funnel_bbbc.py --out <new dir> --figmech <figmech_v39_numbers.json> + the 7 Inputs arguments of figmech_v39 run.sh
import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[_k] = "1"
import argparse, hashlib, json, re, sys
from pathlib import Path
import numpy as np
from scipy.stats import binom, f as fdist
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE / "m3" / "k2prep")]
import m3_replay_lib as L
L_PIN = "ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de"
FIGMECH_PIN = "a4b23923a65d4284b9580529b70a1845565aa0e5460177e55eb8b85c5442945c"
SITE = "bbbc041_all"; SD = 41; EPS = 0.10
FORBID = re.compile(r"plasmo|iml|r8|m5", re.I)


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def icc_anova(x, n):
    """One-way ANOVA ICC point estimate of the field-level exceedance indicator (the b1u statistics without the F-quantile step); clipped to [0, 1)."""
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); ph = x.sum() / N
    n0 = (N - (n ** 2).sum() / N) / (m - 1); pi = x / n; msb = (n * (pi - ph) ** 2).sum() / (m - 1); msw = (n * pi * (1 - pi)).sum() / (N - m)
    if msw <= 0: return 0.0 if msb <= 0 else 1.0
    th = (msb / msw - 1) / n0
    return 0.0 if th <= 0 else float(min(th / (1 + th), 1.0))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--figmech", required=True)
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); out = Path(a.out)
    if any(FORBID.search(str(v)) for v in vars(a).values()): print("GATE (exit 4): forbidden argument", file=sys.stderr); return 4
    if out.exists(): return fail("output directory exists", out)
    if sha(L.__file__) != L_PIN: return fail("m3_replay_lib sha differs")
    if sha(a.figmech) != FIGMECH_PIN: return fail("figmech_v39_numbers.json sha differs")
    FM = json.load(open(a.figmech))
    try:
        inp = L.Inputs(a.a2link, a.ids_iter1, a.mhlh, a.rev1_bbbc, a.rev1_mh, a.rev1_lh, a.mh_rev1_py)
        pop = L.population(inp, SD, SITE)
    except AssertionError as ex: return fail("library", ex)
    s, y = pop["s"], pop["y"]; fl, fl_neg = pop["fl"], pop["fl_neg"]["c0"]
    neg = np.sort(s[y == 0])[::-1]; N = len(neg)
    # anchors against the Figure 2 export (same arm, seed, site, library)
    if N != FM["n_neg_bbbc_all"]: return fail("negative count", N, "!= figmech", FM["n_neg_bbbc_all"])
    fpr_src = float((s[y == 0] >= FM["tau_src_logit"]).mean())
    if fpr_src != FM["fpr_at_tau_src_bbbc_all"]: return fail("FPR at tau_src", repr(fpr_src), "!= figmech", repr(FM["fpr_at_tau_src_bbbc_all"]))
    byname = {str(f): len(fl_neg[f]) for f in fl}
    if [byname.get(f) for f in FM["fields"]] != FM["n_f"]: return fail("n_f of the figmech draw-0 fields differ")
    # site threshold at pooled FPR eps (diagnostic only)
    k = int(round(EPS * N)); tau = float(neg[k - 1]); fpr_site = float((neg >= tau).mean())
    n_f = np.array([len(fl_neg[f]) for f in fl]); X_f = np.array([len(fl_neg[f]) - np.searchsorted(fl_neg[f], tau, side="left") for f in fl])
    if int(X_f.sum()) != int((neg >= tau).sum()) or int(n_f.sum()) != N: return fail("field sums differ from the pooled counts")
    nn = np.arange(1, int(n_f.max()) + 1)
    q05 = binom.ppf(0.05, nn, fpr_site); q95 = binom.ppf(0.95, nn, fpr_site)
    lo_f, hi_f = binom.ppf(0.05, n_f, fpr_site), binom.ppf(0.95, n_f, fpr_site)
    below, above = X_f < lo_f, X_f > hi_f
    exp_out = float(np.mean(binom.cdf(lo_f - 1, n_f, fpr_site) + binom.sf(hi_f, n_f, fpr_site)))
    chi2 = float((((X_f - n_f * fpr_site) ** 2) / (n_f * fpr_site * (1 - fpr_site))).sum()); F = len(fl)
    icc = icc_anova(X_f, n_f); nbw = float((n_f ** 2).sum() / N); de = 1 + (nbw - 1) * icc
    res = dict(nature="P5a D-a: per-field negatives and exceedances at the site threshold with pooled FPR eps (diagnostic, not a rule threshold); BBBC-ALL, arm A2, seed 41",
               script_sha256=sha(__file__), library_sha256=L_PIN, figmech_sha256=FIGMECH_PIN, input_sha256={str(p): h for p, h in inp.used.items()},
               site=SITE, seed=SD, eps=EPS, n_neg=int(N), n_fields=int(F), k=k, tau_eps_logit=tau, fpr_site=fpr_site,
               anchor=dict(fpr_at_tau_src=fpr_src, figmech_fields_checked=len(FM["fields"])),
               n_out_below=int(below.sum()), n_out_above=int(above.sum()), share_out=float((below | above).mean()), share_out_expected_binomial=exp_out,
               pearson_x2=chi2, dispersion_x2_over_df=chi2 / (F - 1), icc_anova=icc, n_bar_w=nbw, design_effect=de,
               n_f_min=int(n_f.min()), n_f_median=float(np.median(n_f)), n_f_max=int(n_f.max()),
               fields=[str(f) for f in fl], n_f=n_f.tolist(), X_f=X_f.tolist(), frac_f=(X_f / n_f).tolist(),
               funnel_n=nn.tolist(), funnel_q05_frac=(q05 / nn).tolist(), funnel_q95_frac=(q95 / nn).tolist())
    payload = json.dumps(res, indent=1, allow_nan=False)
    out.mkdir(parents=True)
    with open(out / "d_a_funnel_bbbc.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh: fh.write(f"{sha(out / 'd_a_funnel_bbbc.json')}  d_a_funnel_bbbc.json\n")
    print(json.dumps({k_: v for k_, v in res.items() if not isinstance(v, (list, dict))}, indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

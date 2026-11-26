# P5 E-batch shared library (spec P5_规格_E1E2E3E5E6a_开发站_20261103 draft3).  LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Pure numpy/scipy: certificate rules on pooled labeled negatives (Q, MC, CP, N0), exact analytic violation laws, population evaluation, seed-cluster interval, output writer.
# Convention everywhere: a cell is a false positive when s >= tau (s = logit1 - logit0); tau = +inf means "no cell flagged" (FPR 0, sensitivity 0, reported as vacuous).
# Never read: M5, IML, PlasmoCount, NIH test fold, any image, any classifier.  Paths containing iml / plasmo are refused.
import hashlib, json, math, sys
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.stats import beta as betad, rankdata

ROOT = Path(__file__).resolve().parents[1]
LIB_PIN = "ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de"     # m3_replay_lib.py
MASTER_PIN = "b7be87050d2915c22bd97af0a7b52391ce060adba42a0138586ac795c23ec10b"
SITES = ("bbbc041_all", "mpidb_full", "lacuna_B"); SEEDS = [41, 42, 43, 44, 45]; EPSS = (0.10, 0.05); DELTA = 0.10; M_FIELDS = 22
MASTER_SITE = {"bbbc041_all": "bbbc041_all", "mpidb_full": "mpidb_full", "lacuna_B": "lacuna_B"}   # survey_arm/<site>|a2|auroc_mean


class AnchorError(Exception):
    """A recomputed value overlapping a stored value differs beyond its tolerance (exit 3, no output directory is created)."""


class CheckError(Exception):
    """An implementation check against an analytic law failed (exit 5, nothing is read)."""


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def refuse_path(*ps):
    for p in ps: assert "iml" not in str(p).lower() and "plasmo" not in str(p).lower(), ("confirmatory-site path refused", str(p))


def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, np.floating): o = float(o)
    if isinstance(o, np.integer): o = int(o)
    if isinstance(o, np.bool_): o = bool(o)
    return None if isinstance(o, float) and not math.isfinite(o) else o


def write_out(out, res, arrays=None):
    """New directory only; files opened with 'x'; SHA256SUMS lists the result file and, when given, <stem>_taus.npz holding `arrays` (per-draw thresholds, +inf allowed)."""
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=False); files = [out]
    with open(out, "x") as f: json.dump(clean(res), f, indent=1, allow_nan=False)
    if arrays:
        z = out.parent / (out.stem + "_taus.npz"); files.append(z)
        with open(z, "xb") as f: np.savez(f, **arrays)
    with open(out.parent / "SHA256SUMS", "x") as f: f.write("".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(out)


def load_master(path):
    assert sha(path) == MASTER_PIN, "results_master.json differs from its pin"
    return json.load(open(path))["values"]


def auroc(s, y):
    """Mann-Whitney AUROC with average ranks (threshold free)."""
    s = np.asarray(s, float); y = np.asarray(y); r = rankdata(s); n1 = int((y == 1).sum()); n0 = int((y == 0).sum())
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def p32(s):
    """float32 probability sigmoid(float32 s) in float32: the score scale of the bbbc041_all survey (323 saturated ties at A2 seed 41)."""
    return (1.0 / (1.0 + np.exp(-np.asarray(s, np.float32)))).astype(np.float32)


def auroc_bases(s, y):
    """AUROC on both historical score scales: s (float64 logit difference) and p32."""
    return {"s": auroc(s, y), "p32": auroc(p32(s), y)}


def anchor_auroc(master, site, bases, tol=1e-12):
    """bases: one auroc_bases() dict per seed.  The master stores survey_arm/<site>|a2|auroc_mean, and the historical pipelines used different score scales (diagnosis of 实验管理, 2026-11-03):
    bbbc041_all = AUROC on p32 (bit-identical), lacuna_B = AUROC on s (bit-identical), mpidb_full = the two agree.  The mean over seeds of each basis is compared with the stored value;
    a basis passes when the difference is <= tol (1e-12); AnchorError when none does.  Returns the basis that matched ('s', 'p32' or 'both') for the output (anchor_auroc_basis)."""
    key = f"survey_arm/{MASTER_SITE[site]}|a2|auroc_mean"; means = {b: float(np.mean([x[b] for x in bases])) for b in ("s", "p32")}
    if key not in master: raise AnchorError((key, None, means, tol))
    hit = [b for b in ("s", "p32") if abs(float(master[key]) - means[b]) <= tol]
    if not hit: raise AnchorError((key, float(master[key]), means, tol))
    return "both" if len(hit) == 2 else hit[0]


# ---------------------------------------------------------------- rules on a pooled negative sample
def _ceil(x): return int(math.ceil(round(x, 9)))            # guards 18.000000000000004-type float noise


def kth(neg, k):
    """k-th smallest (1-based) of the sorted array; +inf above n."""
    return float(neg[k - 1]) if 1 <= k <= len(neg) else math.inf


def k_q(n, eps): return max(1, _ceil(n * (1 - eps)))        # empirical (1 - eps) quantile: the ceil(n(1-eps))-th smallest
def k_mc(n, eps): return _ceil((n + 1) * (1 - eps))         # marginal split-conformal order statistic; k > n -> vacuous


@lru_cache(maxsize=None)
def x_star(n, eps, delta=DELTA):
    """Largest sample false-positive count x in 0..n-1 whose Clopper-Pearson upper limit (level 1 - delta, one-sided) is <= eps; -1 when even x = 0 fails."""
    xs = np.arange(0, n); ub = betad.ppf(1 - delta, xs + 1, n - xs); ok = np.where(ub <= eps)[0]
    return int(ok.max()) if len(ok) else -1


def tau_q(neg, eps): return kth(neg, k_q(len(neg), eps))
def tau_mc(neg, eps): return kth(neg, k_mc(len(neg), eps))


def tau_cp(neg, eps, delta=DELTA):
    """R1 exact rule on the pooled negatives: threshold with at most x* sample false positives; x* = 0 puts tau just above the sample maximum; no x* -> +inf (cannot certify)."""
    x = x_star(len(neg), eps, delta)
    if x < 0: return math.inf
    return float(np.nextafter(neg[-1], math.inf)) if x == 0 else kth(neg, len(neg) - x + 1)


RULES = {"Q": tau_q, "MC": tau_mc, "CP": tau_cp}


def viol_analytic(rule, n, eps, delta=DELTA):
    """P(true FPR > eps) for iid continuous scores, exact: the true FPR at the j-th order statistic is Beta(n+1-j, j).  None when the rule is vacuous at this n."""
    if rule == "Q": k = k_q(n, eps); return float(betad.sf(eps, n + 1 - k, k))
    if rule == "MC":
        k = k_mc(n, eps); return None if k > n else float(betad.sf(eps, n + 1 - k, k))
    if rule == "CP":
        x = x_star(n, eps, delta)
        if x < 0: return None
        return float(betad.sf(eps, 1, n)) if x == 0 else float(betad.sf(eps, x, n - x + 1))
    raise ValueError(rule)


# ---------------------------------------------------------------- population
def make_pop(s, y, g):
    s = np.asarray(s, float); y = np.asarray(y); g = np.asarray(g)
    fl = [f for f in np.unique(g) if ((g == f) & (y == 0)).any()]
    return dict(sneg=np.sort(s[y == 0]), spos=np.sort(s[y == 1]), fl=fl, negf={f: s[(g == f) & (y == 0)] for f in fl}, N=int((y == 0).sum()), P=int((y == 1).sum()))


def fpr_at(pop, tau): return 0.0 if math.isinf(tau) else float(1 - np.searchsorted(pop["sneg"], tau, side="left") / pop["N"])
def sens_at(pop, tau): return 0.0 if math.isinf(tau) else float(1 - np.searchsorted(pop["spos"], tau, side="left") / pop["P"])
def oracle_own(pop, eps): return sens_at(pop, float(np.quantile(pop["sneg"], 1 - eps)))


def evaluate(pop, tau, eps, orac):
    f = fpr_at(pop, tau); sn = sens_at(pop, tau)
    return dict(vac=bool(math.isinf(tau)), fpr=f, viol=bool(f > eps), sens=sn, dso=sn - orac)


VAC_CONVENTION = ("vacuous draw = tau +inf (no cell flagged): FPR 0 and counted as covered in coverage / fpr_*; coverage_nonvacuous and fpr_mean_nonvacuous exclude those draws. "
                  "NOT the manuscript convention for cannot-certify (evaluated at the fallback threshold): never put the two in one table.")
CI_DEFINITION = "seed-cluster 90% interval: percentile 5-95 of the mean coverage over B=10000 resamples of the five per-seed coverages with replacement (default_rng([20261103, 99]); equals a seed-cluster bootstrap when every seed has the same number of draws)"


def summarise(recs):
    """recs: list of evaluate() dicts of one cell.  coverage = 1 - violation rate (a vacuous draw never violates; its share is reported)."""
    v = np.array([r["viol"] for r in recs], float); vac = np.array([r["vac"] for r in recs]); fpr = np.array([r["fpr"] for r in recs]); dso = np.array([r["dso"] for r in recs])
    nv = ~vac
    return dict(n_draws=len(recs), coverage=float(1 - v.mean()), vacuous_share=float(vac.mean()), ds_own_nonvacuous=float(dso[nv].mean()) if nv.any() else None,
                coverage_nonvacuous=float(1 - v[nv].mean()) if nv.any() else None, fpr_mean_nonvacuous=float(fpr[nv].mean()) if nv.any() else None,
                fpr_mean=float(fpr.mean()), fpr_median=float(np.median(fpr)), fpr_p90=float(np.percentile(fpr, 90)), convention=VAC_CONVENTION)


def pooled_summary(recs_by_seed):
    """recs_by_seed: {seed: [recs]} -> summary over all draws + per-seed coverage + seed-cluster 90% interval."""
    allr = [r for rs in recs_by_seed.values() for r in rs]; out = summarise(allr)
    cov = [summarise(rs)["coverage"] for rs in recs_by_seed.values()]; out["coverage_per_seed"] = cov; out["coverage_ci90_seedcluster"] = seed_cluster_ci(cov)
    return out


def seed_cluster_ci(cov, B=10000):
    """Percentile 90% interval of the mean coverage when seeds are resampled with replacement (CI_DEFINITION; fixed stream)."""
    c = np.asarray(cov, float); rng = np.random.default_rng([20261103, 99]); idx = rng.integers(0, len(c), (B, len(c))); m = c[idx].mean(axis=1)
    return [float(np.percentile(m, 5)), float(np.percentile(m, 95))]


def check_analytic(rows, ndraw_total, slack):
    """rows: list of (rule, n, eps, observed violation rate over ndraw_total pooled draws).  CheckError listing every cell whose gap exceeds 3 MC standard errors + slack, and every CP cell above delta + that tolerance."""
    bad = []
    for rule, n, eps, obs in rows:
        a = viol_analytic(rule, n, eps)
        if a is None: continue
        tol = 3 * math.sqrt(max(a * (1 - a), 1e-12) / ndraw_total) + slack
        if abs(obs - a) > tol: bad.append((rule, n, eps, obs, a, tol))
        if rule == "CP" and obs > DELTA + 3 * math.sqrt(DELTA * (1 - DELTA) / ndraw_total) + slack: bad.append(("CP>delta", n, eps, obs, DELTA))
    if bad: raise CheckError(bad)


def common_args(ap, with_nihval=False):
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "out", "master"): ap.add_argument("--" + k, required=True)
    if with_nihval: ap.add_argument("--nihval_dir", required=True)
    ap.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    return ap


def open_inputs(a):
    sys.path.insert(0, str(ROOT / "m3" / "k2prep")); import m3_replay_lib as L
    assert sha(L.__file__) == LIB_PIN, "m3_replay_lib differs from its pin"
    refuse_path(a.out, a.master)
    return L, L.Inputs(a.a2link, a.ids_iter1, a.mhlh, a.rev1_bbbc, a.rev1_mh, a.rev1_lh, a.mh_rev1_py)

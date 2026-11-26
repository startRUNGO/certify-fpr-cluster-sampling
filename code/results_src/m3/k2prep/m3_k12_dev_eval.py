# P5 mechanism selection on the development sites (problem card draft5 b3318a71, sections 0, 6, 6A, with addendum A): arms S3, c1s (= B1u_c1s), K1+S3, K2+S3 (main, bootstrap-t),
# K1t+S3, K2chi+S3 (controls, report only).  LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 approves.
# NATURE: design-side mechanism selection on development sites; not confirmatory evidence.  Confirmation only on W1 (M5).  IML and PlasmoCount are not read.
# DATA: the three development sites via m3_replay_lib (score files checked against their manifests; pinned inputs).  The REV-1 draw loop is re-run here (copy of
#   m3_replay_lib.replay, which copies d_icc0.py) with the K arms added; m3_replay_lib.self_check on the S3, B1u_c0 and B1u_c1s rows of THIS loop must pass (1e-12 vs
#   stored REV-1) before anything is written, which also verifies that this copy consumes the REV-1 random stream exactly.
# K ARMS (draft5 section 6A; each candidate k at alpha = delta/3, union over c0, c1s, c2s, choice as S3 = largest share of labelled positives above the certified
#   threshold; fixed-sequence test on the S3 grid; "certified at t" iff the upper confidence bound UB_Z(t) <= 0 for Z(t) = sum_F (x_f - eps n_f)):
#   per drawn field i: z_i = x_i - eps n_i, w_i = g_i - eps c_i (g = all cells >= t, c = all cells; no label), e_i = z_i - w_i; F = fields of the site (>= 1 negative);
#   sum_F w_f and sigma_w = population sd of w_f over F are exact (computed from all cells of the site); fpc = sqrt(1 - m/|F|).
#   K1   : Zhat = sum_F w + |F| mean(e); SE = |F| sd(e, ddof 1)/sqrt(m) fpc; bootstrap-t over the m fields, B = 999, T* = (mean* - mean)/(sd*/sqrt(m)),
#          UB = Zhat - q_alpha(T*) SE.   K1t: UB = Zhat + t_{m-1, 1-alpha} SE.
#   K2   : Zhat2 = |F| mean(z); sigma_U = sigma_w + sd_U(e); UB = Zhat2 + z_{1-alpha_m} |F| sigma_U/sqrt(m) fpc; alpha_m = 0.8 alpha, alpha_d = 0.2 alpha;
#          sd_U(e) by bootstrap-t on log s^2 with SE(log s^2) = sqrt((kappa - 1)/m), kappa = m4/s^4 (non-excess), B = 999, log s2_U = log s2 - q_{alpha_d}(T*) SE.
#   K2chi: sd_U(e) = sd(e) sqrt((m - 1)/chi2_{alpha_d, m-1}).
#   Addendum A: sd(e) = 0 -> that grid point is NOT certified (UB = +inf), counted; bootstrap resamples with sd* = 0 are dropped and the share recorded, and if the
#   dropped share exceeds 5% the grid point is NOT certified, counted.  Interface to the library's fixed_sequence (ub(j) <= eps): to_seq(UB) returns 0.0 iff
#   UB <= 0 and 1.0 otherwise (including +inf and NaN); with eps in {0.05, 0.10}, ub(j) <= eps iff UB_Z <= 0.  UB_Z is never converted to the FPR scale.
#   Identity used in the methods text (addendum A3): e_f = -P_f (TPR_f(t) - eps), P_f = positives of field f.
#   Bootstrap weights: default_rng([20261101, site index, seed, draw, eps code, candidate code]) once per (draw, eps, candidate), shared along the grid and by K1/K2
#   (common random numbers); separate from, and not consuming, the REV-1 stream.  Exceptions in a K bound are counted per (site, eps, arm) and re-raised (the
#   fixed-sequence test then treats that candidate as not certified, as REV-1; at most one exception per draw and candidate, so up to 3 per draw); a cell with
#   exceptions > 1% of its draws is flagged and the arm is reference only.
# CRITERIA (draft5 section 6, fixed before any development-site run), for K1 and K2 (main arms) against S3:
#   gate, every one of the 6 cells (3 sites x 2 eps): coverage >= 0.90, seed-cluster bootstrap 90% interval lower limit >= 0.88 (B = 10,000, seed 20261101),
#     cannot_certify <= S3 + 0.02;  main: three-site mean improvement of dSens_own at eps 0.10 >= +0.02, and no cell worse than S3 by more than 0.01;
#   no flagged cell.  Both pass -> first choice is the larger eps-0.10 mean improvement.  Controls K1t, K2chi and the c1s baseline are reported, not judged.
# Output: one new JSON in a new directory ("x", allow_nan=False).  Prints only the output path.
# Usage: python -I -B m3_k12_dev_eval.py --a2link .. --ids_iter1 .. --mhlh .. --rev1_bbbc .. --rev1_mh .. --rev1_lh .. --mh_rev1_py .. --out <dir/json>
import argparse, json, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import chi2, norm, t as tdist

sys.path.insert(0, str(Path(__file__).resolve().parent)); import m3_replay_lib as L   # noqa: E402

DELTA = 0.10; ALPHA = DELTA/3; A_M = 0.8*ALPHA; A_D = 0.2*ALPHA; B = 999; BOOT_SEED = 20261101; CI_B = 10000; CI_SEED = 20261101
K_ARMS = ("K1", "K1t", "K2", "K2chi"); MAIN = ("K1", "K2"); CAND_CODE = {"c0": 0, "c1s": 1, "c2s": 2}
GATE_COV, GATE_LO, GATE_CC, MAIN_GAIN, MAIN_WORSE, EXC_SHARE = 0.90, 0.88, 0.02, 0.02, 0.01, 0.01
EXC = defaultdict(Counter); DEG = defaultdict(Counter); DROP = defaultdict(list); MAX_DROP = 0.05


def to_seq(U):
    """Addendum A2: 0.0 iff UB_Z <= 0, else 1.0 (also for +inf and NaN)."""
    return 0.0 if (U is not None and np.isfinite(U) and U <= 0) else 1.0


class AllCells:
    """All cells of the site's fields F for one candidate: exact sum_F w_f(t) and population sd of w_f(t), and per-field sorted scores for g_i(t)."""
    def __init__(self, pop, k):
        fl = list(pop["fl"]); idx = {f: i for i, f in enumerate(fl)}; g = pop["g"]; sc = pop["sc"][k]
        keep = np.isin(g, np.asarray(fl)); s = sc[keep]; gi = np.array([idx[v] for v in g[keep]]); self.F = len(fl)
        self.c = np.bincount(gi, minlength=self.F).astype(float); self.idx = idx
        o = np.lexsort((s, gi)); bounds = np.cumsum(self.c.astype(int))[:-1]; self.by_field = np.split(s[o], bounds)
        order = np.argsort(-s, kind="stable"); fd = gi[order]
        grp = np.argsort(fd, kind="stable"); fds = fd[grp]; first = np.searchsorted(fds, fds, side="left"); rank = np.empty(len(fd), int); rank[grp] = np.arange(len(fd)) - first
        self.G2 = np.concatenate([[0.0], np.cumsum(2.0*rank + 1)]); self.GC = np.concatenate([[0.0], np.cumsum(self.c[fd])])
        self.asc = np.sort(s); self.n = len(s); self.Sc = float(self.c.sum()); self.Sc2 = float((self.c**2).sum())

    def w_stats(self, t, eps):
        K = self.n - np.searchsorted(self.asc, t, side="left")
        Sw = K - eps*self.Sc; Sw2 = self.G2[K] - 2*eps*self.GC[K] + eps**2*self.Sc2; var = Sw2/self.F - (Sw/self.F)**2
        return float(Sw), float(np.sqrt(max(var, 0.0)))

    def g_lab(self, fields, t):
        return np.array([len(self.by_field[self.idx[f]]) - np.searchsorted(self.by_field[self.idx[f]], t, side="left") for f in fields], float)

    def c_lab(self, fields):
        return np.array([self.c[self.idx[f]] for f in fields])


def degenerate_tol(e):
    """Shared degeneracy tolerance (approved 2026-10-30): s_e <= 1e-9 max(1, max|e_i|) counts as s_e = 0 (floating-point rounding, not real spread)."""
    return 1e-9*max(1.0, float(np.abs(e).max()))


def k1_ub(e, Sw, F, m, W, mode, key):
    eb = e.mean(); se = e.std(ddof=1); fpc = np.sqrt(1 - m/F); Zh = Sw + F*eb; tol = degenerate_tol(e)
    if se <= tol: DEG[key]["se0"] += 1; return np.inf
    SE = F*se/np.sqrt(m)*fpc
    if mode == "t": return Zh + tdist.ppf(1 - ALPHA, m - 1)*SE
    mu = (W*e).sum(1)/m; ex2 = (W*e*e).sum(1)/m; s2 = (ex2 - mu**2)*m/(m - 1); ok = s2 > tol**2; drop = 1 - ok.mean(); DROP[key].append(drop)
    if drop > MAX_DROP: DEG[key]["drop>5%"] += 1; return np.inf
    T = (mu[ok] - eb)/(np.sqrt(s2[ok])/np.sqrt(m)); return Zh - np.quantile(T, ALPHA)*SE


def sd_upper(e, m, W, mode, key):
    se = e.std(ddof=1); tol = degenerate_tol(e)
    if se <= tol: DEG[key]["se0"] += 1; return np.inf
    if mode == "chi": return float(se*np.sqrt((m - 1)/chi2.ppf(A_D, m - 1)))
    mu0 = e.mean(); v0 = e.var(); k0 = ((e - mu0)**4).mean()/v0**2; SE = np.sqrt(max(k0 - 1, 1e-12)/m); ls2 = np.log(se**2)
    mu = (W*e).sum(1)/m; e2 = (W*e**2).sum(1)/m; e3 = (W*e**3).sum(1)/m; e4 = (W*e**4).sum(1)/m
    v = e2 - mu**2; ok = v > tol**2; drop = 1 - ok.mean(); DROP[key].append(drop)
    if drop > MAX_DROP: DEG[key]["drop>5%"] += 1; return np.inf
    m4 = e4 - 4*mu*e3 + 6*mu**2*e2 - 3*mu**4; kap = m4[ok]/v[ok]**2; SEs = np.sqrt(np.maximum(kap - 1, 1e-12)/m)
    T = (np.log(v[ok]*m/(m - 1)) - ls2)/SEs; return float(np.sqrt(np.exp(ls2 - np.quantile(T, A_D)*SE)))


def k2_ub(z, e, sw, F, m, W, mode, key):
    fpc = np.sqrt(1 - m/F); sU = sw + sd_upper(e, m, W, mode, key)
    return F*z.mean() + norm.ppf(1 - A_M)*F*sU/np.sqrt(m)*fpc


def counted(arm, site, eps, fn):
    def g(j):
        try: return fn(j)
        except Exception as ex:
            EXC[(site, eps, arm)][type(ex).__name__] += 1; raise
    return g


def run_site_seed(inp, sd, site):
    pop = L.population(inp, sd, site); fl, fl_neg, fl_posv = pop["fl"], pop["fl_neg"], pop["fl_posv"]; sens_at, fpr_at = pop["sens_at"], pop["fpr_at"]
    AC = {k: AllCells(pop, k) for k in L.CEN}
    key8, sidx = L.SITECFG[site]; ai = L.AI; rows = []; m = L.M
    rng = np.random.default_rng([20261022, key8, ai, sd, m, sidx])
    for d_ in range(L.NDRAW):
        sel_ = rng.choice(len(fl), m, replace=False); F = [fl[i] for i in sel_]; clicks = sum(len(fl_posv["c0"][f]) for f in F); ni = np.array([len(fl_neg["c0"][f]) for f in F])
        for eps, dl in L.SETTINGS:   # kept: consumes the stream exactly as REV-1 (R1-2 subsets)
            r_ = np.random.default_rng(int(rng.integers(1 << 30))); [r_.choice(n_, 2, replace=False) if n_ >= 2 else np.arange(n_) for n_ in ni]
        D = {}
        for k in L.CEN:
            negs = [fl_neg[k][f] for f in F]; pooled = np.concatenate(negs); grid = np.unique(np.quantile(pooled, 1-np.arange(L.G)/L.G))[::-1]
            X = np.stack([ni[i]-np.searchsorted(negs[i], grid, side="left") for i in range(m)]).astype(float)
            posl = np.concatenate([fl_posv[k][f] for f in F]) if clicks > 0 else np.array([]); D[k] = (negs, pooled, grid, X, posl)
        cache = {}

        def ubf(rule, k, j):
            if (rule, k, j) not in cache: cache[(rule, k, j)] = L.b1u(D[k][3][:, j], ni, 0.10 if rule == "U1" else 0.10/3)
            return cache[(rule, k, j)]

        def rec(name, eps, k, t, cert):
            sn = sens_at(k, t)
            rows.append(dict(site=site, seed=sd, draw=d_, eps=eps, rule=name, cand=k, cert=bool(cert), fpr=float(fpr_at(k, t)), sens=float(sn),
                             dse=float(sn-pop["orac_unc"][eps]), dso=float(sn-pop["orac_own"][(k, eps)]), ni=ni.tolist(), fields=[str(f) for f in F]))

        def choose(T, eps, name):
            cert_c = [k for k in L.CEN if T[k][1]]
            if not cert_c: rec(name, eps, "c0", T["c0"][0], False); return
            if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
            else:
                sl = [(D[k][4] >= T[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
            rec(name, eps, ch, T[ch][0], True)
        for eps, dl in L.SETTINGS:
            choose({k: L.fixed_sequence(D[k][2], lambda j, k=k: ubf("U3", k, j), eps) for k in L.CEN}, eps, "S3")
            for k in ("c0", "c1s"):
                t_, c_, _ = L.fixed_sequence(D[k][2], lambda j, k=k: ubf("U1", k, j), eps); rec("B1u_" + k, eps, k, t_, c_)
            Wk = {k: np.random.default_rng([BOOT_SEED, L.DEV.index(site), sd, d_, int(round(eps*100)), CAND_CODE[k]]).multinomial(m, np.full(m, 1/m), size=B).astype(float)
                  for k in L.CEN}
            cl = {k: AC[k].c_lab(F) for k in L.CEN}
            for arm in K_ARMS:
                def ub(j, k, arm=arm, eps=eps):
                    t = D[k][2][j]; z = D[k][3][:, j] - eps*ni; w = AC[k].g_lab(F, t) - eps*cl[k]; e = z - w; Sw, sw = AC[k].w_stats(t, eps); key = (site, eps, arm)
                    if arm in ("K1", "K1t"): U = k1_ub(e, Sw, AC[k].F, m, Wk[k], "t" if arm == "K1t" else "bt", key)
                    else: U = k2_ub(z, e, sw, AC[k].F, m, Wk[k], "chi" if arm == "K2chi" else "bt", key)
                    return to_seq(U)
                choose({k: L.fixed_sequence(D[k][2], counted(arm, site, eps, lambda j, k=k: ub(j, k)), eps) for k in L.CEN}, eps, arm)
    return rows


def seed_ci(R, eps, rng):
    hits = np.array([sum(r["fpr"] <= eps for r in R if r["seed"] == s) for s in L.SEEDS], float); ns = np.array([sum(r["seed"] == s for r in R) for s in L.SEEDS], float)
    idx = rng.integers(0, len(L.SEEDS), size=(CI_B, len(L.SEEDS))); v = hits[idx].sum(1)/ns[idx].sum(1)
    return [float(np.quantile(v, 0.05)), float(np.quantile(v, 0.95))]


def summarise(R, eps, site, arm, rng):
    fpr = np.array([r["fpr"] for r in R]); cert = np.array([r["cert"] for r in R]); ex = EXC.get((site, eps, arm), Counter()); nex = int(sum(ex.values()))
    cls = {}
    for c in ("c0", "c1s", "c2s", "uncertified"):
        S = [r for r in R if (not r["cert"] if c == "uncertified" else (r["cert"] and r["cand"] == c))]
        cls[c] = dict(n=len(S), coverage=(float(np.mean([r["fpr"] <= eps for r in S])) if S else None))
    dr = DROP.get((site, eps, arm), [])
    return dict(n=len(R), coverage=float((fpr <= eps).mean()), coverage_ci90_seedcluster=seed_ci(R, eps, rng), cannot_certify=float(1 - cert.mean()),
                sens=float(np.mean([r["sens"] for r in R])), ds_own=float(np.mean([r["dso"] for r in R])), ds_eps=float(np.mean([r["dse"] for r in R])),
                fpr_mean=float(fpr.mean()), fpr_median=float(np.median(fpr)), fpr_sd_certified=(float(fpr[cert].std()) if cert.any() else None), by_class=cls,
                n_bound_exceptions=nex, exception_types=[[k, v] for k, v in ex.most_common(3)], flag=("exceptions>1%" if nex > EXC_SHARE*len(R) else None),
                n_se0=int(DEG.get((site, eps, arm), Counter()).get("se0", 0)), n_drop_over_5pct=int(DEG.get((site, eps, arm), Counter()).get("drop>5%", 0)),
                boot_drop_share_median=(float(np.median(dr)) if dr else None))


def judge(cells):
    get = lambda arm, s, e: [c for c in cells if c["arm"] == arm and c["site"] == s and abs(c["eps"] - e) < 1e-12][0]
    out = {}
    for arm in MAIN:
        gate, worse, gains = [], [], []
        for s in L.DEV:
            for e, _ in L.SETTINGS:
                a, b = get(arm, s, e), get("S3", s, e)
                gate.append(dict(site=s, eps=e, cov_ok=a["coverage"] >= GATE_COV, lo_ok=a["coverage_ci90_seedcluster"][0] >= GATE_LO,
                                 cc_ok=a["cannot_certify"] <= b["cannot_certify"] + GATE_CC + 1e-12, flagged=a["flag"] is not None))
                worse.append(dict(site=s, eps=e, diff=a["ds_own"] - b["ds_own"]))
                if abs(e - 0.10) < 1e-12: gains.append(a["ds_own"] - b["ds_own"])
        g_ok = all(x["cov_ok"] and x["lo_ok"] and x["cc_ok"] and not x["flagged"] for x in gate); mg = float(np.mean(gains))
        m_ok = mg >= MAIN_GAIN and all(x["diff"] >= -MAIN_WORSE for x in worse)
        out[arm] = dict(gate=gate, gate_pass=g_ok, mean_gain_eps010=mg, per_cell_diff=worse, main_pass=m_ok, passes=bool(g_ok and m_ok))
    passed = [a for a in MAIN if out[a]["passes"]]
    out["choice"] = (max(passed, key=lambda a: out[a]["mean_gain_eps010"]) if passed else None)
    out["if_none"] = "neither main arm passes -> back to root-cause analysis, criteria not relaxed"
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); out = Path(a.out); assert not out.exists(), "output exists"
    inp = L.Inputs(a.a2link, a.ids_iter1, a.mhlh, a.rev1_bbbc, a.rev1_mh, a.rev1_lh, a.mh_rev1_py)
    EXC.clear(); DEG.clear(); DROP.clear(); rows = []
    for site in L.DEV:
        for sd in L.SEEDS: rows += run_site_seed(inp, sd, site)
    bad = L.self_check(inp, rows)
    if bad: raise SystemExit("self-check failed, nothing written: " + repr(bad[:5]))
    rng = np.random.default_rng(CI_SEED); arms = ["S3", "B1u_c1s"] + list(K_ARMS)
    cells = [dict(arm=arm, site=s, eps=e, **summarise([r for r in rows if r["site"] == s and abs(r["eps"] - e) < 1e-12 and r["rule"] == arm], e, s, arm, rng))
             for arm in arms for s in L.DEV for e, _ in L.SETTINGS]
    rep = dict(script_sha256=L.sha(__file__), lib_sha256=L.sha(L.__file__), inputs_sha256=dict(L.PIN), score_files_sha256=inp.used,
               problem_card="draft5 b3318a71 + addendum A", nature="development-site mechanism selection; not confirmatory evidence",
               settings=dict(delta=DELTA, alpha_per_candidate=ALPHA, alpha_m=A_M, alpha_d=A_D, B=B, boot_seed=BOOT_SEED, ci_B=CI_B, ci_seed=CI_SEED, lambda_K1=1.0),
               self_check="passed (S3, B1u_c0, B1u_c1s of this loop vs stored REV-1 within 1e-12)", cells=cells, judgement=judge(cells))
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False)
    out.parent.mkdir(parents=False, exist_ok=False)
    with open(out, "x") as fh: fh.write(txt)
    print("written", out)


if __name__ == "__main__":
    main()

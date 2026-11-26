"""S5 per PREREG_REV_S5_MIDDLEGROUND_20261024.md (sha256 242f7093...). Adapted from s4.py (sha256 91a34154...): rules S3 (control), S5a, S5b only; S5a/S5b change only the n_eff step of the upper bound. CPU only."""
import json, sys, numpy as np
from pathlib import Path
from multiprocessing import Pool
from scipy.special import betainc
from scipy.stats import beta as betad, f as fdist, norm as normd, t as tdist
INF = Path(sys.argv[1]); E7 = Path(sys.argv[2]); OUT = Path(sys.argv[3]); REV1 = Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=False)
ARMS = ["a2", "a0-ctrl"]; SEEDS = [41,42,43,44,45]; NDRAW = 200; SETTINGS = [(0.10, 0.10), (0.05, 0.10)]; G = 200
SITES = ["bbbc041_all", "mpidb_full", "lacuna_B"]; CEN = ["c0", "c1s", "c2s"]; MS_DEV = [10, 22]; RULES = ["S3", "S5a", "S5b"]
def b1u(x, n, delta):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X/N
    n0 = (N-(n**2).sum()/N)/(m-1); pi = x/n; msb = (n*(pi-ph)**2).sum()/(m-1); msw = (n*pi*(1-pi)).sum()/(N-m)
    if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0
    else:
        F = msb/msw; flo = fdist.ppf(delta, m-1, N-m); th = (F/flo-1)/n0; iccu = 0.0 if th <= 0 else th/(1+th)
    nbw = (n**2).sum()/N; neff = N/(1+(nbw-1)*iccu); xeff = X*neff/N
    return float(betad.ppf(1-delta, xeff+1, neff-xeff))
def ub_s5(rule, x, n, delta):
    """S5a: S4b design effect (linearized field-level variance), then n_eff *= (z_{1-delta}/t_{m-1,1-delta})^2 (Korn-Graubard); p_hat in {0,1}: DEFF=1, df factor still applied.
    S5b: as S3 (ICC upper limit via one-sided F interval) but the F quantile is taken at level 0.10 instead of delta; Clopper-Pearson stays at delta. delta = 0.10/3."""
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X/N
    if rule == "S5a":
        if ph <= 0 or ph >= 1: deff = 1.0
        else:
            v = m/(m-1)*((x-ph*n)**2).sum()/N**2; deff = max(1.0, v/(ph*(1-ph)/N))
        neff = N/deff*(normd.ppf(1-delta)/tdist.ppf(1-delta, m-1))**2
        xeff = ph*neff
    else:
        n0 = (N-(n**2).sum()/N)/(m-1); pi = x/n; msb = (n*(pi-ph)**2).sum()/(m-1); msw = (n*pi*(1-pi)).sum()/(N-m)
        if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0
        else:
            F = msb/msw; flo = fdist.ppf(0.10, m-1, N-m); th = (F/flo-1)/n0; iccu = 0.0 if th <= 0 else th/(1+th)
        nbw = (n**2).sum()/N; neff = N/(1+(nbw-1)*iccu); xeff = X*neff/N
    return float(betad.ppf(1-delta, xeff+1, neff-xeff))
def R1(n, eps, delta):
    for k in range(1, int(n)+1):
        if betainc(n-k+1, k, eps) >= 1-delta: return k, True
    return int(n), False
def load_site(arm, sd, site):
    if site == "lacuna_B":
        lid = np.load(E7/"ids_lacuna_B.npz", allow_pickle=True); z = np.load(E7/f"{arm}_s{sd}_lacuna_B_logits.npz")["logits"]; return (z[:, 1]-z[:, 0]).astype(float), lid["y"], lid["field"]
    ids = np.load(INF/"ids_iter1.npz", allow_pickle=True); z = np.load(INF/f"{arm}_s{sd}_scores.npz")[site+"_logits"]; return (z[:, 1]-z[:, 0]).astype(float), ids[site+"_y"], ids[site+"_group"]
def cstat(v, kind, gm):
    if kind == "c0": return 0.0
    w = len(v)/(len(v)+20.0)
    if kind == "c1s": return float(gm+w*(np.median(v)-gm))
    return float(gm+w*(np.median(np.sort(v)[:int(np.ceil(0.7*len(v)))])-gm))
def task(a):
    ai, sd, si = a; arm = ARMS[ai]; site = SITES[si]; key8 = 8; sidx = si
    s, y, g = load_site(arm, sd, site); fl = list(np.unique(g)); gm = float(np.median(s)); cs = {f: {k: cstat(s[g == f], k, gm) for k in CEN} for f in fl}
    sc = {k: s-np.array([cs[f][k] for f in g]) for k in CEN}
    fl_neg = {k: {f: np.sort(sc[k][(g == f) & (y == 0)]) for f in fl} for k in CEN}; fl_posv = {k: {f: sc[k][(g == f) & (y == 1)] for f in fl} for k in CEN}
    fl = [f for f in fl if len(fl_neg["c0"][f]) > 0]
    sneg = {k: np.sort(sc[k][y == 0]) for k in CEN}; spos = {k: np.sort(sc[k][y == 1]) for k in CEN}; N = len(sneg["c0"])
    sens_at = lambda k, t: 1-np.searchsorted(spos[k], t, side="left")/len(spos[k]); fpr_at = lambda k, t: 1-np.searchsorted(sneg[k], t, side="left")/N
    orac_unc = {e: sens_at("c0", np.quantile(sneg["c0"], 1-e)) for e in (0.10, 0.05)}; orac_own = {(k, e): sens_at(k, np.quantile(sneg[k], 1-e)) for k in CEN for e in (0.10, 0.05)}; old05 = sens_at("c0", np.quantile(sneg["c0"], 0.95))
    # Proposition 1 conservativeness at the oracle_0.10 threshold (population level, uncentered)
    tau = np.quantile(sneg["c0"], 0.90); pbar = float(np.mean([ (fl_neg["c0"][f] >= tau).mean() for f in fl])); q = float(np.mean([ (fl_neg["c0"][f].max() >= tau) for f in fl]))
    res = {}
    def rec(name, m, eps, k, t, cert, clicks, nlab):
        e = res.setdefault((arm, site, m, name, eps), {"fpr": [], "sens": [], "dse": [], "dso": [], "dold": [], "cert": [], "clicks": [], "nlab": []}); sn = sens_at(k, t)
        e["fpr"].append(fpr_at(k, t)); e["sens"].append(sn); e["dse"].append(sn-orac_unc[eps]); e["dso"].append(sn-orac_own[(k, eps)]); e["dold"].append(sn-old05); e["cert"].append(cert); e["clicks"].append(clicks); e["nlab"].append(nlab)
    for m in MS_DEV:
        if m > len(fl): continue
        rng = np.random.default_rng([20261022, key8, ai, sd, m, sidx])
        for _ in range(NDRAW):
            sel_ = rng.choice(len(fl), m, replace=False); F = [fl[i] for i in sel_]; clicks = sum(len(fl_posv["c0"][f]) for f in F); ni = np.array([len(fl_neg["c0"][f]) for f in F]); subs = {}
            for eps, dl in SETTINGS:
                r_ = np.random.default_rng(int(rng.integers(1 << 30))); subs[eps] = [r_.choice(n_, 2, replace=False) if n_ >= 2 else np.arange(n_) for n_ in ni]
            D = {}
            for k in CEN:
                negs = [fl_neg[k][f] for f in F]; pooled = np.concatenate(negs); grid = np.unique(np.quantile(pooled, 1-np.arange(G)/G))[::-1]
                X = np.stack([ni[i]-np.searchsorted(negs[i], grid, side="left") for i in range(m)]).astype(float)
                posl = np.concatenate([fl_posv[k][f] for f in F]) if clicks > 0 else np.array([]); D[k] = (negs, pooled, grid, X, posl)
            cache = {}
            def ub(rule, k, j):
                key = (rule, k, j)
                if key not in cache:
                    x = D[k][3][:, j]
                    cache[key] = b1u(x, ni, 0.10/3) if rule == "S3" else ub_s5(rule, x, ni, 0.10/3)
                return cache[key]
            def tau_of(rule, k, eps):
                grid = D[k][2]; j = 0
                try:
                    while j < len(grid) and ub(rule, k, j) <= eps: j += 1
                except Exception: j = 0
                return (grid[j-1], True) if j > 0 else (grid[0], False)
            for eps, dl in SETTINGS:
                nlab_all = len(D["c0"][1])
                for rule in RULES:
                    T = {k: tau_of(rule, k, eps) for k in CEN}
                    cert_c = [k for k in CEN if T[k][1]]
                    if not cert_c: rec(rule, m, eps, "c0", T["c0"][0], False, clicks, nlab_all)
                    else:
                        if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
                        else:
                            sl = [(D[k][4] >= T[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
                        rec(rule, m, eps, ch, T[ch][0], True, clicks, nlab_all)
    print("done", arm, sd, site, flush=True)
    return res, {"pbar": pbar, "q": q}
if __name__ == "__main__":
    tasks = [(ai, sd, si) for ai in range(2) for sd in SEEDS for si in range(len(SITES))]
    with Pool(4) as p: parts = p.map(task, tasks, chunksize=1)
    merged = {}
    for (ai, sd, si), (r, pr) in zip(tasks, parts):
        for k, e in r.items():
            d = merged.setdefault(k, {kk: [] for kk in e})
            for kk in d: d[kk] += e[kk]
    agg = {}
    for (arm, site, m, name, eps), e in merged.items():
        fpr = np.array(e["fpr"])
        agg[f"{arm}|{site}|{m}|{name}|{eps}"] = {"cov": float((fpr <= eps).mean()), "cov_lb010": float((np.maximum(0, fpr-0.10) <= eps).mean()), "cannot_certify": float(1-np.mean(e["cert"])), "ds_eps": float(np.mean(e["dse"])), "ds_own": float(np.mean(e["dso"])), "ds_old05": float(np.mean(e["dold"])),
            "fpr_mean": float(fpr.mean()), "fpr_median": float(np.median(fpr)), "fpr_p90": float(np.quantile(fpr, 0.9)), "sens": float(np.mean(e["sens"])), "clicks": float(np.mean(e["clicks"])), "n_labeled_neg": float(np.mean(e["nlab"])), "n": len(fpr)}
    paired = {}
    for (arm, site, m, name, eps), e in merged.items():
        if name == "S3": continue
        c = merged[(arm, site, m, "S3", eps)]
        paired[f"{arm}|{site}|{m}|{name}-S3|{eps}"] = {"sens_diff": float(np.mean(np.array(e["sens"])-np.array(c["sens"]))), "cov_diff": float((np.array(e["fpr"]) <= eps).mean()-(np.array(c["fpr"]) <= eps).mean())}
    json.dump({"agg": agg, "paired": paired}, open(OUT/"s5_results.json", "w"), indent=1)
    ref = json.load(open(REV1))["agg"]; bad = []
    for k, v in agg.items():
        if "|S3|" in k:
            for kk, vv in v.items():
                if ref[k][kk] != vv: bad.append((k, kk, ref[k][kk], vv))
    nchk = sum(len(v) for k, v in agg.items() if "|S3|" in k)
    json.dump({"s3_values_compared": nchk, "mismatches": len(bad), "first": bad[:5]}, open(OUT/"selfcheck.json", "w"), indent=1)
    print("SELF-CHECK S3 vs rev1:", nchk, "values compared,", len(bad), "mismatches", flush=True)
    if bad: print(bad[:5]); sys.exit(7)
    print("saved all")

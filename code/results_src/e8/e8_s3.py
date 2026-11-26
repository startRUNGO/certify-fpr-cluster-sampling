"""E8 PlasmoCount per PREREG_E8_PLASMOCOUNT_20261023.md (sha256 2ce88dfd...). Frozen rule S3 (s11.py c7d3ad62...), only the loader and rng key changed. Run once."""
import json, sys, numpy as np
from pathlib import Path
from multiprocessing import Pool
from scipy.special import betainc
from scipy.stats import beta as betad, f as fdist
INF = Path(sys.argv[1]); E7 = Path(sys.argv[2]); OUT = Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
ARMS = ["a2", "a0-ctrl"]; SEEDS = [41,42,43,44,45]; NDRAW = 200; MS = [22]; SETTINGS = [(0.10, 0.10), (0.05, 0.10)]; G = 200
SITES = ["plasmo"]; CEN = ["c0", "c1s", "c2s"]
def b1u(x, n, delta):
    """B1u at level delta (CP level delta, ICC_U F-interval one-sided level delta)"""
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X/N
    n0 = (N-(n**2).sum()/N)/(m-1); pi = x/n
    msb = (n*(pi-ph)**2).sum()/(m-1); msw = (n*pi*(1-pi)).sum()/(N-m)
    if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0
    else:
        F = msb/msw; flo = fdist.ppf(delta, m-1, N-m); th = (F/flo-1)/n0; iccu = 0.0 if th <= 0 else th/(1+th)
    nbw = (n**2).sum()/N; neff = N/(1+(nbw-1)*iccu); xeff = X*neff/N
    return float(betad.ppf(1-delta, xeff+1, neff-xeff))
def R1(n, eps, delta):
    for k in range(1, int(n)+1):
        if betainc(n-k+1, k, eps) >= 1-delta: return k, True
    return int(n), False
def load_site(arm, sd, site):
    lid = np.load(E7/"ids_plasmo.npz", allow_pickle=True)
    z = np.load(E7/f"{arm}_s{sd}_plasmo_logits.npz")["logits"]; s = z[:, 1]-z[:, 0]; y = lid["y"]; g = lid["field"]
    return s.astype(float), y, g
def cstat(v, kind, gm):
    if kind == "c0": return 0.0
    w = len(v)/(len(v)+20.0)
    if kind == "c1s": return float(gm+w*(np.median(v)-gm))
    return float(gm+w*(np.median(np.sort(v)[:int(np.ceil(0.7*len(v)))])-gm))
RULES = {"S3": (3, CEN), "S2": (2, ["c0", "c1s"]), "U_c0": (1, ["c0"]), "U_c1s": (1, ["c1s"]), "U_c2s": (1, ["c2s"])}
def task(a):
    ai, sd, si = a; arm = ARMS[ai]; site = SITES[si]
    s, y, g = load_site(arm, sd, site); fl = list(np.unique(g)); gm = float(np.median(s)); cs = {f: {k: cstat(s[g == f], k, gm) for k in CEN} for f in fl}
    sc = {k: s-np.array([cs[f][k] for f in g]) for k in CEN}
    fl_neg = {k: {f: np.sort(sc[k][(g == f) & (y == 0)]) for f in fl} for k in CEN}; fl_posv = {k: {f: sc[k][(g == f) & (y == 1)] for f in fl} for k in CEN}
    fl = [f for f in fl if len(fl_neg["c0"][f]) > 0]
    sneg = {k: np.sort(sc[k][y == 0]) for k in CEN}; spos = {k: np.sort(sc[k][y == 1]) for k in CEN}; N = len(sneg["c0"])
    sens_at = lambda k, t: 1-np.searchsorted(spos[k], t, side="left")/len(spos[k]); fpr_at = lambda k, t: 1-np.searchsorted(sneg[k], t, side="left")/N
    ref05 = {k: sens_at(k, np.quantile(sneg[k], 0.95)) for k in CEN}; refe = {(k, e): sens_at(k, np.quantile(sneg[k], 1-e)) for k in CEN for e in (0.10, 0.05)}
    res = {}
    def rec(name, m, eps, k, t, cert, sel, clicks, nlab):
        e = res.setdefault((arm, site, m, name, eps), {"fpr": [], "sens": [], "dso": [], "dse": [], "dsc0": [], "cert": [], "sel": [], "clicks": [], "nlab": []})
        sn = sens_at(k, t); e["fpr"].append(fpr_at(k, t)); e["sens"].append(sn); e["dso"].append(sn-ref05[k]); e["dse"].append(sn-refe[(k, eps)]); e["dsc0"].append(sn-ref05["c0"]); e["cert"].append(cert); e["sel"].append(sel); e["clicks"].append(clicks); e["nlab"].append(nlab)
    for m in MS:
        rng = np.random.default_rng([20261022, 13, ai, sd, m, si])
        for _ in range(NDRAW):
            sel_ = rng.choice(len(fl), m, replace=False); F = [fl[i] for i in sel_]; clicks = sum(len(fl_posv["c0"][f]) for f in F)
            ni = np.array([len(fl_neg["c0"][f]) for f in F]); subs = {}
            for eps, dl in SETTINGS:
                r_ = np.random.default_rng(int(rng.integers(1 << 30))); subs[eps] = [r_.choice(n_, 2, replace=False) if n_ >= 2 else np.arange(n_) for n_ in ni]
            D = {}
            for k in CEN:
                negs = [fl_neg[k][f] for f in F]; pooled = np.concatenate(negs); grid = np.unique(np.quantile(pooled, 1-np.arange(G)/G))[::-1]
                X = np.stack([ni[i]-np.searchsorted(negs[i], grid, side="left") for i in range(m)]).astype(float)
                posl = np.concatenate([fl_posv[k][f] for f in F]) if clicks > 0 else np.array([]); D[k] = (negs, pooled, grid, X, posl)
            cache = {}
            def tau(k, K, eps):
                grid, X = D[k][2], D[k][3]; j = 0
                while j < len(grid):
                    key = (k, K, j)
                    if key not in cache: cache[key] = b1u(X[:, j], ni, 0.10/K)
                    if cache[key] > eps: break
                    j += 1
                return (grid[j-1], True) if j > 0 else (grid[0], False)
            for eps, dl in SETTINGS:
                T = {(k, K): tau(k, K, eps) for K in (1, 2, 3) for k in CEN if not (K == 2 and k == "c2s")}
                for name, (K, cands) in RULES.items():
                    cert_c = [k for k in cands if T[(k, K)][1]]
                    if not cert_c:
                        k0 = cands[0] if len(cands) == 1 else "c0"; rec(name, m, eps, k0, T[(k0, K)][0], False, -1, clicks, len(D[k0][1])); continue
                    if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
                    else:
                        sens_l = [(D[k][4] >= T[(k, K)][0]).mean() for k in cert_c]; mxs = max(sens_l); ch = cert_c[[i for i, v in enumerate(sens_l) if v == mxs][0]]
                    rec(name, m, eps, ch, T[(ch, K)][0], True, CEN.index(ch), clicks, len(D[ch][1]))
                k2, ok2 = R1(2*m, eps, dl); negs = D["c0"][0]; two = np.sort(np.concatenate([negs[i][subs[eps][i]] for i in range(m)])); rec("R1_2c0", m, eps, "c0", two[min(k2, len(two))-1], ok2, 0, clicks, 2*m)
    print("done", arm, sd, site, flush=True)
    return res
if __name__ == "__main__":
    tasks = [(ai, sd, si) for ai in range(2) for sd in SEEDS for si in range(3)]
    with Pool(8) as p: parts = p.map(task, tasks, chunksize=1)
    merged = {}
    for r in parts:
        for k, e in r.items():
            d = merged.setdefault(k, {kk: [] for kk in e})
            for kk in d: d[kk] += e[kk]
    agg = {}
    for (arm, site, m, name, eps), e in merged.items():
        fpr = np.array(e["fpr"]); sel = np.array(e["sel"])
        agg[f"{arm}|{site}|{m}|{name}|{eps}"] = {"cov": float((fpr <= eps).mean()), "cov_lb010": float((np.maximum(0, fpr-0.10) <= eps).mean()), "fpr_mean": float(fpr.mean()), "sens": float(np.mean(e["sens"])), "ds_own05": float(np.mean(e["dso"])), "ds_ownEps": float(np.mean(e["dse"])), "ds_c0ref05": float(np.mean(e["dsc0"])),
            "cannot_certify": float(1-np.mean(e["cert"])), "sel_freq": {"c0": float((sel == 0).mean()), "c1s": float((sel == 1).mean()), "c2s": float((sel == 2).mean()), "none": float((sel == -1).mean())}, "clicks": float(np.mean(e["clicks"])), "n_labeled_neg": float(np.mean(e["nlab"])), "n": len(fpr)}
    json.dump(agg, open(OUT/"e8_plasmo_results.json", "w"), indent=1); print("saved all")

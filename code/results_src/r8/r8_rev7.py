"""R8 IML per PREREG_R8_IML_20261023.md; adapted from rev7.py (only site list and data loading changed; REV-7 definitions, 06f6b5b5...). Post hoc: S3-G (source-fixed grid) and S3-G-ext (extended fixed grid). CPU only."""
import json, sys, numpy as np
from pathlib import Path
from multiprocessing import Pool
from scipy.special import betainc
from scipy.stats import beta as betad, f as fdist
sys.path.insert(0, str(Path(__file__).parent)); import s8_bounds as B8
INF = Path(sys.argv[1]); E7 = Path(sys.argv[2]); E8 = Path(sys.argv[3]); OUT = Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=True)
ARMS = ["a2", "a0-ctrl"]; SEEDS = [41,42,43,44,45]; NDRAW = 200; SETTINGS = [(0.10, 0.10), (0.05, 0.10)]; G = 200
SITES = ["iml"]; CEN = ["c0", "c1s", "c2s"]; MS_DEV = [10, 22]
def b1u(x, n, delta):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X/N
    n0 = (N-(n**2).sum()/N)/(m-1); pi = x/n; msb = (n*(pi-ph)**2).sum()/(m-1); msw = (n*pi*(1-pi)).sum()/(N-m)
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
    if site == "iml":
        lid = np.load(E8/"ids_iml.npz", allow_pickle=True); z = np.load(E8/f"{arm}_s{sd}_iml_logits.npz")["logits"]; return (z[:, 1]-z[:, 0]).astype(float), lid["y"], lid["field"]
    if site == "lacuna_B":
        lid = np.load(E7/"ids_lacuna_B.npz", allow_pickle=True); z = np.load(E7/f"{arm}_s{sd}_lacuna_B_logits.npz")["logits"]; return (z[:, 1]-z[:, 0]).astype(float), lid["y"], lid["field"]
    ids = np.load(INF/"ids_iter1.npz", allow_pickle=True); z = np.load(INF/f"{arm}_s{sd}_scores.npz")[site+"_logits"]; return (z[:, 1]-z[:, 0]).astype(float), ids[site+"_y"], ids[site+"_group"]
def cstat(v, kind, gm):
    if kind == "c0": return 0.0
    w = len(v)/(len(v)+20.0)
    if kind == "c1s": return float(gm+w*(np.median(v)-gm))
    return float(gm+w*(np.median(np.sort(v)[:int(np.ceil(0.7*len(v)))])-gm))
def task(a):
    ai, sd, si = a; arm = ARMS[ai]; site = SITES[si]; key8 = 14; sidx = 0
    s, y, g = load_site(arm, sd, site); fl = list(np.unique(g)); gm = float(np.median(s)); cs = {f: {k: cstat(s[g == f], k, gm) for k in CEN} for f in fl}
    sc = {k: s-np.array([cs[f][k] for f in g]) for k in CEN}
    fl_neg = {k: {f: np.sort(sc[k][(g == f) & (y == 0)]) for f in fl} for k in CEN}; fl_posv = {k: {f: sc[k][(g == f) & (y == 1)] for f in fl} for k in CEN}
    fl = [f for f in fl if len(fl_neg["c0"][f]) > 0]
    sneg = {k: np.sort(sc[k][y == 0]) for k in CEN}; spos = {k: np.sort(sc[k][y == 1]) for k in CEN}; N = len(sneg["c0"])
    sens_at = lambda k, t: 1-np.searchsorted(spos[k], t, side="left")/len(spos[k]); fpr_at = lambda k, t: 1-np.searchsorted(sneg[k], t, side="left")/N
    orac_unc = {e: sens_at("c0", np.quantile(sneg["c0"], 1-e)) for e in (0.10, 0.05)}; orac_own = {(k, e): sens_at(k, np.quantile(sneg[k], 1-e)) for k in CEN for e in (0.10, 0.05)}; old05 = sens_at("c0", np.quantile(sneg["c0"], 0.95))
    # Proposition 1 conservativeness at the oracle_0.10 threshold (population level, uncentered)
    tau = np.quantile(sneg["c0"], 0.90); pbar = float(np.mean([ (fl_neg["c0"][f] >= tau).mean() for f in fl])); q = float(np.mean([ (fl_neg["c0"][f].max() >= tau) for f in fl]))
    zs = np.load(INF/f"{arm}_s{sd}_scores.npz")["nih_validation_logits"]; ssrc = (zs[:, 1]-zs[:, 0]).astype(float); idsrc = np.load(INF/"ids_iter1.npz", allow_pickle=True); ysrc = idsrc["nih_validation_y"]; gsrc = idsrc["nih_validation_group"]; gms = float(np.median(ssrc))
    csrc = {f: {k: cstat(ssrc[gsrc == f], k, gms) for k in CEN} for f in np.unique(gsrc)}; GRID = {}; GRIDX = {}
    for k in CEN:
        sk = ssrc-np.array([csrc[f][k] for f in gsrc]); sn = sk[ysrc == 0]; gr = np.unique(np.quantile(sn, 1-np.arange(G)/G))[::-1]; iq = float(np.percentile(sn, 75)-np.percentile(sn, 25))
        GRID[k] = gr; GRIDX[k] = np.concatenate([gr[0]+np.arange(40, 0, -1)*0.5*iq, gr])
    res = {}
    def rec(name, m, eps, k, t, cert, clicks, nlab):
        e = res.setdefault((arm, site, m, name, eps), {"fpr": [], "sens": [], "dse": [], "dso": [], "dold": [], "cert": [], "clicks": [], "nlab": []}); sn = sens_at(k, t)
        e["fpr"].append(fpr_at(k, t)); e["sens"].append(sn); e["dse"].append(sn-orac_unc[eps]); e["dso"].append(sn-orac_own[(k, eps)]); e["dold"].append(sn-old05); e["cert"].append(cert); e["clicks"].append(clicks); e["nlab"].append(nlab)
    for m in [22]:
        if m > len(fl): continue
        rng = np.random.default_rng([20261022, key8, ai, sd, m, sidx])
        for _ in range(NDRAW):
            sel_ = rng.choice(len(fl), m, replace=False); F = [fl[i] for i in sel_]; clicks = sum(len(fl_posv["c0"][f]) for f in F); ni = np.array([len(fl_neg["c0"][f]) for f in F]); subs = {}
            for eps, dl in SETTINGS:
                r_ = np.random.default_rng(int(rng.integers(1 << 30))); subs[eps] = [r_.choice(n_, 2, replace=False) if n_ >= 2 else np.arange(n_) for n_ in ni]
            negsC = {k: [fl_neg[k][f] for f in F] for k in CEN}; posC = {k: (np.concatenate([fl_posv[k][f] for f in F]) if clicks > 0 else np.array([])) for k in CEN}; cache = {}
            def tau_of(var, k, eps):
                grid = (GRID if var == "G" else GRIDX)[k]; negs = negsC[k]; j = 0
                try:
                    while j < len(grid):
                        key = (var, k, j)
                        if key not in cache:
                            x = np.array([ni[i]-np.searchsorted(negs[i], grid[j], side="left") for i in range(m)], float); cache[key] = b1u(x, ni, 0.10/3)
                        if cache[key] > eps: break
                        j += 1
                except Exception: j = 0
                return (grid[j-1], True) if j > 0 else (grid[0], False)
            nl = int(ni.sum())
            for eps, dl in SETTINGS:
                for var, nm in (("G", "S3G"), ("X", "S3Gx")):
                    T3 = {k: tau_of(var, k, eps) for k in CEN}; cert_c = [k for k in CEN if T3[k][1]]
                    if not cert_c: rec(nm, m, eps, "c0", T3["c0"][0], False, clicks, nl)
                    else:
                        if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
                        else:
                            sl = [(posC[k] >= T3[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
                        rec(nm, m, eps, ch, T3[ch][0], True, clicks, nl)
    print("done", arm, sd, site, flush=True)
    return res, {"pbar": pbar, "q": q}
if __name__ == "__main__":
    tasks = [(ai, sd, si) for ai in range(2) for sd in SEEDS for si in range(len(SITES))]
    with Pool(4) as p: parts = p.map(task, tasks, chunksize=1)
    merged = {}; prop1 = {}
    for (ai, sd, si), (r, pr) in zip(tasks, parts):
        prop1.setdefault(f"{ARMS[ai]}|{SITES[si]}", []).append(pr)
        for k, e in r.items():
            d = merged.setdefault(k, {kk: [] for kk in e})
            for kk in d: d[kk] += e[kk]
    agg = {}
    for (arm, site, m, name, eps), e in merged.items():
        fpr = np.array(e["fpr"])
        agg[f"{arm}|{site}|{m}|{name}|{eps}"] = {"cov": float((fpr <= eps).mean()), "cov_lb010": float((np.maximum(0, fpr-0.10) <= eps).mean()), "cannot_certify": float(1-np.mean(e["cert"])), "ds_eps": float(np.mean(e["dse"])), "ds_own": float(np.mean(e["dso"])), "ds_old05": float(np.mean(e["dold"])),
            "fpr_mean": float(fpr.mean()), "fpr_median": float(np.median(fpr)), "fpr_p90": float(np.quantile(fpr, 0.9)), "sens": float(np.mean(e["sens"])), "clicks": float(np.mean(e["clicks"])), "n_labeled_neg": float(np.mean(e["nlab"])), "n": len(fpr)}
    p1 = {k: {"pbar": float(np.mean([x["pbar"] for x in v])), "q": float(np.mean([x["q"] for x in v])), "ratio_q_over_pbar": float(np.mean([x["q"]/x["pbar"] for x in v if x["pbar"] > 0]))} for k, v in prop1.items()}
    json.dump({"agg": agg, "prop1": p1}, open(OUT/"rev7_results.json", "w"), indent=1); print("saved all")

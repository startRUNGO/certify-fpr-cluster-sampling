"""REV-4 per PREREG_REV_PARASITEMIA_20261023.md (sha256 df7e2d59...). pi-stratified sensitivity of frozen S3 and injection simulation. CPU only."""
import json, sys, numpy as np
from pathlib import Path
from multiprocessing import Pool
from scipy.stats import beta as betad, f as fdist
INF = Path(sys.argv[1]); E7 = Path(sys.argv[2]); E8 = Path(sys.argv[3]); OUT = Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=True)
ARMS = ["a2", "a0-ctrl"]; SEEDS = [41, 42, 43, 44, 45]; NDRAW = 200; EPS = 0.10; G = 200; M = 22; SITES = ["bbbc041_all", "mpidb_full", "lacuna_B", "plasmo"]; CEN = ["c0", "c1s", "c2s"]; PIS = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50]; NSIM = 40
def b1u(x, n, delta):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X/N
    n0 = (N-(n**2).sum()/N)/(m-1); pi = x/n; msb = (n*(pi-ph)**2).sum()/(m-1); msw = (n*pi*(1-pi)).sum()/(N-m)
    if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0
    else:
        F = msb/msw; flo = fdist.ppf(delta, m-1, N-m); th = (F/flo-1)/n0; iccu = 0.0 if th <= 0 else th/(1+th)
    nbw = (n**2).sum()/N; neff = N/(1+(nbw-1)*iccu); xeff = X*neff/N
    return float(betad.ppf(1-delta, xeff+1, neff-xeff))
def load_site(arm, sd, site):
    if site == "plasmo":
        lid = np.load(E8/"ids_plasmo.npz", allow_pickle=True); z = np.load(E8/f"{arm}_s{sd}_plasmo_logits.npz")["logits"]; return (z[:, 1]-z[:, 0]).astype(float), lid["y"], lid["field"]
    if site == "lacuna_B":
        lid = np.load(E7/"ids_lacuna_B.npz", allow_pickle=True); z = np.load(E7/f"{arm}_s{sd}_lacuna_B_logits.npz")["logits"]; return (z[:, 1]-z[:, 0]).astype(float), lid["y"], lid["field"]
    ids = np.load(INF/"ids_iter1.npz", allow_pickle=True); z = np.load(INF/f"{arm}_s{sd}_scores.npz")[site+"_logits"]; return (z[:, 1]-z[:, 0]).astype(float), ids[site+"_y"], ids[site+"_group"]
def cstat(v, kind, gm, raw=False):
    if kind == "c0": return 0.0
    w = 1.0 if raw else len(v)/(len(v)+20.0)
    if kind == "c1s": return float(gm+w*(np.median(v)-gm))
    return float(gm+w*(np.median(np.sort(v)[:int(np.ceil(0.7*len(v)))])-gm))
def task(a):
    ai, sd, si = a; arm = ARMS[ai]; site = SITES[si]; s, y, g = load_site(arm, sd, site); fl = list(np.unique(g)); gm = float(np.median(s))
    cs = {f: {k: cstat(s[g == f], k, gm) for k in CEN} for f in fl}; sc = {k: s-np.array([cs[f][k] for f in g]) for k in CEN}
    fl_neg = {k: {f: np.sort(sc[k][(g == f) & (y == 0)]) for f in fl} for k in CEN}; fl_posv = {k: {f: sc[k][(g == f) & (y == 1)] for f in fl} for k in CEN}
    fl = [f for f in fl if len(fl_neg["c0"][f]) > 0]; npos = {f: len(fl_posv["c0"][f]) for f in fl}; nneg = {f: len(fl_neg["c0"][f]) for f in fl}
    pif = {f: npos[f]/(npos[f]+nneg[f]) for f in fl if npos[f] > 0}; pv = np.array(list(pif.values())); q33, q67 = np.quantile(pv, [1/3, 2/3])
    desc = {"median": float(np.median(pv)), "p90": float(np.quantile(pv, 0.9)), "max": float(pv.max()), "frac_gt_0.10": float((pv > 0.10).mean()), "frac_gt_0.30": float((pv > 0.30).mean()), "n_fields_with_pos": len(pv)}
    # positives of the site with their field pi
    pos_idx = np.where(y == 1)[0]; pos_pi = np.array([pif.get(g[i], np.nan) for i in pos_idx]); pos_sc = {k: sc[k][pos_idx] for k in CEN}
    bins = {"t1": pos_pi <= q33, "t2": (pos_pi > q33) & (pos_pi <= q67), "t3": pos_pi > q67, "gt10": pos_pi > 0.10, "gt30": pos_pi > 0.30, "all": np.ones(len(pos_pi), bool)}
    rng = np.random.default_rng([20261022, 15, ai, sd, si])
    simcand = [f for f in fl if nneg[f] >= 20]; simF = list(rng.choice(len(simcand), min(NSIM, len(simcand)), replace=False)); simF = [simcand[i] for i in simF]; fl_cert = [f for f in fl if f not in set(simF)]
    # injection set-up (fixed per task)
    sim = []
    for f in simF:
        negraw = s[(g == f) & (y == 0)]; other = s[(y == 1) & (g != f)]; n = len(negraw); per = {}
        for pi in PIS:
            a_ = int(round(pi/(1-pi)*n)) if pi > 0 else 0; inj = rng.choice(other, a_, replace=True) if a_ > 0 else np.array([]); cells = np.concatenate([negraw, inj])
            c = {k: cstat(cells, k, gm) for k in CEN}; craw = {k: cstat(cells, k, gm, raw=True) for k in CEN}; lev = {k: float((negraw < craw[k]).mean()) for k in ("c1s", "c2s")}
            per[pi] = {"inj": inj, "c": c, "lev": lev}
        sim.append(per)
    acc = {}; acc_sim = {}
    def add(d, key, v, c=1):
        e = d.setdefault(key, [0.0, 0]); e[0] += v; e[1] += c
    for _ in range(NDRAW):
        selc = rng.choice(len(fl_cert), M, replace=False); F = [fl_cert[i] for i in selc]; ni = np.array([nneg[f] for f in F]); D = {}
        for k in CEN:
            negs = [fl_neg[k][f] for f in F]; pooled = np.concatenate(negs); grid = np.unique(np.quantile(pooled, 1-np.arange(G)/G))[::-1]
            X = np.stack([ni[i]-np.searchsorted(negs[i], grid, side="left") for i in range(M)]).astype(float); posl = np.concatenate([fl_posv[k][f] for f in F]); D[k] = (grid, X, posl)
        cache = {}
        def tau_of(rule, k):
            grid, X, _ = D[k]; j = 0
            try:
                while j < len(grid):
                    key = (rule, k, j)
                    if key not in cache: cache[key] = b1u(X[:, j], ni, 0.10 if rule == "U1" else 0.10/3)
                    if cache[key] > EPS: break
                    j += 1
            except Exception: j = 0
            return (grid[j-1], True) if j > 0 else (grid[0], False)
        T1 = {k: tau_of("U1", k) for k in CEN}; T3 = {k: tau_of("U3", k) for k in CEN}; cert_c = [k for k in CEN if T3[k][1]]
        if cert_c:
            if len(cert_c) == 1: ch = cert_c[0]
            else:
                sl = [(D[k][2] >= T3[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
        else: ch = None
        rules = {"S3": (ch, T3[ch][0]) if ch else None}
        for k in CEN: rules["B1u_"+k] = (k, T1[k][0]) if T1[k][1] else None
        for rn, rv in rules.items():
            if rv is None: continue
            k, t = rv
            for bn, mk in bins.items():
                if mk.sum() > 0: add(acc, (rn, bn), float((pos_sc[k][mk] >= t).mean()))
            for fi, per in enumerate(sim):
                for pi in PIS[1:]:
                    inj = per[pi]["inj"]; add(acc_sim, (rn, pi), float((inj-per[pi]["c"][k] >= t).mean()))
    out = {"desc": desc, "bins": {bn: int(mk.sum()) for bn, mk in bins.items()}, "sens": {f"{rn}|{bn}": [e[0], e[1]] for (rn, bn), e in acc.items()}, "sim_sens": {f"{rn}|{pi}": [e[0], e[1]] for (rn, pi), e in acc_sim.items()},
           "shift": {f"{k}|{pi}": float(np.mean([per[pi]["c"][k]-per[0.0]["c"][k] for per in sim])) for k in ("c1s", "c2s") for pi in PIS}, "level": {f"{k}|{pi}": float(np.mean([per[pi]["lev"][k] for per in sim])) for k in ("c1s", "c2s") for pi in PIS}}
    print("done", arm, sd, site, flush=True); return out
if __name__ == "__main__":
    tasks = [(ai, sd, si) for ai in range(2) for sd in SEEDS for si in range(len(SITES))]
    with Pool(4) as p: parts = p.map(task, tasks, chunksize=1)
    agg = {}
    for (ai, sd, si), o in zip(tasks, parts):
        key = f"{ARMS[ai]}|{SITES[si]}"; e = agg.setdefault(key, {"desc": [], "bins": [], "sens": {}, "sim_sens": {}, "shift": {}, "level": {}})
        e["desc"].append(o["desc"]); e["bins"].append(o["bins"])
        for dn in ("sens", "sim_sens"):
            for k, (sm, c) in o[dn].items():
                z = e[dn].setdefault(k, [0.0, 0]); z[0] += sm; z[1] += c
        for dn in ("shift", "level"):
            for k, v in o[dn].items(): e[dn].setdefault(k, []).append(v)
    final = {}
    for key, e in agg.items():
        final[key] = {"desc": {k: float(np.mean([d[k] for d in e["desc"]])) for k in e["desc"][0]}, "bins_positives": {k: float(np.mean([b[k] for b in e["bins"]])) for k in e["bins"][0]},
                      "sens": {k: v[0]/v[1] for k, v in e["sens"].items()}, "sim_sens": {k: v[0]/v[1] for k, v in e["sim_sens"].items()}, "shift": {k: float(np.mean(v)) for k, v in e["shift"].items()}, "level": {k: float(np.mean(v)) for k, v in e["level"].items()}}
        cost = {f"{pi}": final[key]["sim_sens"][f"B1u_c0|{pi}"]-final[key]["sim_sens"][f"S3|{pi}"] for pi in PIS[1:] if f"B1u_c0|{pi}" in final[key]["sim_sens"]}; final[key]["centering_cost"] = cost
    trig = {s: float(final[f"a2|{s}"]["centering_cost"]["0.2"]) for s in SITES[:3]}; json.dump({"agg": final, "trigger_a2_cost_at_pi0.20": trig, "triggered": bool(any(v > 0.05 for v in trig.values()))}, open(OUT/"rev4_results.json", "w"), indent=1)
    print("trigger costs (A2, pi=0.20):", trig, "triggered:", any(v > 0.05 for v in trig.values())); print("saved all")

"""REV-5 per PREREG_REV_UNLABELED_20261023.md (sha256 de5e3960...). UPA (source prevalence / oracle site prevalence), prevalence-constrained quantile; realized FPR and sensitivity. CPU only."""
import json, sys, numpy as np
from pathlib import Path
INF = Path(sys.argv[1]); E7 = Path(sys.argv[2]); E8 = Path(sys.argv[3]); OUT = Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=True)
SITES = ["bbbc041_all", "mpidb_full", "lacuna_B", "plasmo", "nih_locked_test"]; SEEDS = [41, 42, 43, 44, 45]; ARMS = ["a2", "a0-ctrl"]; EPS = [0.10, 0.05]
ids = np.load(INF/"ids_iter1.npz", allow_pickle=True); lid = np.load(E7/"ids_lacuna_B.npz", allow_pickle=True); pid = np.load(E8/"ids_plasmo.npz", allow_pickle=True)
def scores(stem, site):
    if site == "lacuna_B": z = np.load(E7/f"{stem}_lacuna_B_logits.npz")["logits"]; y = lid["y"]
    elif site == "plasmo": z = np.load(E8/f"{stem}_plasmo_logits.npz")["logits"]; y = pid["y"]
    else: z = np.load(INF/f"{stem}_scores.npz")[site+"_logits"]; y = ids[site+"_y"]
    return (z[:, 1]-z[:, 0]).astype(float), y
def wcdf(vals, w, x):
    o = np.argsort(vals); v = vals[o]; c = np.cumsum(w[o]); c = c/c[-1]; i = np.searchsorted(v, x, side="right"); return 0.0 if i == 0 else float(c[i-1])
rows = []
for arm in ARMS:
    for sd in SEEDS:
        stem = f"{arm}_s{sd}"; thr = json.load(open(INF/f"{stem}_eval.json"))["source_threshold"]["value"]
        zs = np.load(INF/f"{stem}_scores.npz")["nih_validation_logits"]; ss = zs[:, 1]-zs[:, 0]; ys = ids["nih_validation_y"]
        for site in SITES:
            s, y = scores(stem, site); pi = float((y == 1).mean()); fp = lambda t: float((s[y == 0] >= t).mean()); se = lambda t: float((s[y == 1] >= t).mean())
            res = {"source_threshold": (fp(thr), se(thr))}
            # UPA, source prevalence: reference = all NIH validation cells
            u0 = wcdf(ss, np.ones(len(ss)), thr); t = np.quantile(s, u0); res["upa_source_prev"] = (fp(t), se(t))
            # UPA, oracle site prevalence: reference = reweighted mixture of source neg/pos
            w = np.where(ys == 1, pi/(ys == 1).sum(), (1-pi)/(ys == 0).sum()); u1 = wcdf(ss, w, thr); t = np.quantile(s, u1); res["upa_oracle_prev"] = (fp(t), se(t))
            for e in EPS:
                lvl = 1-pi-(1-pi)*e; t = np.quantile(s, lvl); res[f"prev_quantile_{e}"] = (fp(t), se(t))
            for k, (f_, s_) in res.items(): rows.append({"arm": arm, "seed": sd, "site": site, "method": k, "fpr": f_, "sens": s_, "pi_site": pi})
out = {"rows": rows, "agg": {}}
meths = sorted({r["method"] for r in rows})
for site in SITES:
    for arm in ARMS:
        for m in meths:
            R = [r for r in rows if r["site"] == site and r["arm"] == arm and r["method"] == m]; f = np.array([r["fpr"] for r in R]); s_ = np.array([r["sens"] for r in R])
            out["agg"][f"{arm}|{site}|{m}"] = {"fpr_mean": float(f.mean()), "fpr_sd": float(f.std(ddof=1)), "sens_mean": float(s_.mean()), "sens_sd": float(s_.std(ddof=1)), "frac_fpr_le_0.10": float((f <= 0.10).mean()), "frac_fpr_le_0.05": float((f <= 0.05).mean()), "pi_site": float(R[0]["pi_site"])}
json.dump(out, open(OUT/"rev5_results.json", "w"), indent=1)
for site in SITES:
    print(site, {m: (round(out["agg"][f"a2|{site}|{m}"]["fpr_mean"], 3), round(out["agg"][f"a2|{site}|{m}"]["sens_mean"], 3)) for m in meths})

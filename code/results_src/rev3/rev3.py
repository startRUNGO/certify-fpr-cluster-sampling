"""REV-3 per PREREG_REV_FRONTEND_20261023.md (sha256 fcfb7f09...). Front-end design vs source-threshold FPR inflation; score-scale decomposition. CPU only, saved logits."""
import json, sys, glob, os, itertools, numpy as np
from pathlib import Path
INF = Path(sys.argv[1]); E7 = Path(sys.argv[2]); E8 = Path(sys.argv[3]); OUT = Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=True)
SITES = ["bbbc041_all", "mpidb_full", "lacuna_B", "plasmo", "nih_locked_test"]; SEEDS = [41, 42, 43, 44, 45]
ARMS = ["b0", "a0-ctrl", "c0", "c1", "a2", "a3", "a4"]
ids = np.load(INF/"ids_iter1.npz", allow_pickle=True); lid = np.load(E7/"ids_lacuna_B.npz", allow_pickle=True); pid = np.load(E8/"ids_plasmo.npz", allow_pickle=True)
def scores(stem, site):
    if site == "lacuna_B": z = np.load(E7/f"{stem}_lacuna_B_logits.npz")["logits"]; y = lid["y"]
    elif site == "plasmo": z = np.load(E8/f"{stem}_plasmo_logits.npz")["logits"]; y = pid["y"]
    else: z = np.load(INF/f"{stem}_scores.npz")[site+"_logits"]; y = ids[site+"_y"]
    return (z[:, 1]-z[:, 0]).astype(float), y
q25 = lambda x: np.percentile(x, 25); q75 = lambda x: np.percentile(x, 75)
rows = []
for arm in ARMS:
    for sd in SEEDS:
        stem = f"{arm}_s{sd}"; thr = json.load(open(INF/f"{stem}_eval.json"))["source_threshold"]["value"]
        zs = np.load(INF/f"{stem}_scores.npz")["nih_validation_logits"]; ss = zs[:, 1]-zs[:, 0]; ys = ids["nih_validation_y"]; sn = ss[ys == 0]; ms = float(np.median(sn)); iqs = float(q75(sn)-q25(sn))
        for site in SITES:
            s, y = scores(stem, site); tn = s[y == 0]; mt = float(np.median(tn)); iqt = float(q75(tn)-q25(tn))
            fpr = float((tn >= thr).mean()); shift = mt-ms; fl = float(((tn-shift) >= thr).mean()); fls = float(((((tn-shift)-ms)*(iqs/iqt)+ms) >= thr).mean())
            rows.append({"arm": arm, "seed": sd, "site": site, "fpr": fpr, "shift_over_iqr_s": shift/iqs, "iqr_ratio": iqt/iqs, "fpr_loc_removed": fl, "fpr_loc_scale_removed": fls})
def perm_p(d):
    d = np.array(d, float); obs = abs(d.sum()); n = len(d); c = 0
    for signs in itertools.product([1, -1], repeat=n): c += abs((np.array(signs)*d).sum()) >= obs-1e-12
    return c/2**n
out = {"rows": rows, "by_arm": {}, "paired": {}}
for site in SITES:
    for arm in ARMS:
        R = [r for r in rows if r["site"] == site and r["arm"] == arm]
        out["by_arm"][f"{site}|{arm}"] = {k: {"mean": float(np.mean([r[k] for r in R])), "sd": float(np.std([r[k] for r in R], ddof=1))} for k in ("fpr", "shift_over_iqr_s", "iqr_ratio", "fpr_loc_removed", "fpr_loc_scale_removed")}
    def per_seed(group): return [np.mean([r["fpr"] for r in rows if r["site"] == site and r["arm"] in group and r["seed"] == sd]) for sd in SEEDS]
    comps = {"a2_minus_a0ctrl": (["a2"], ["a0-ctrl"]), "frontend_group (a2,a3,a4 minus b0,a0-ctrl,c0,c1)": (["a2", "a3", "a4"], ["b0", "a0-ctrl", "c0", "c1"]), "c1_minus_a0ctrl": (["c1"], ["a0-ctrl"]), "a2_minus_c1": (["a2"], ["c1"]), "pseudosite_group (c1,a2,a3,a4 minus b0,a0-ctrl)": (["c1", "a2", "a3", "a4"], ["b0", "a0-ctrl"])}
    for nm, (g1, g2) in comps.items():
        d = np.array(per_seed(g1))-np.array(per_seed(g2)); out["paired"][f"{site}|{nm}"] = {"mean_diff": float(d.mean()), "sd_diff": float(d.std(ddof=1)), "n_positive": int((d > 0).sum()), "p_signflip": float(perm_p(d)), "diffs": [float(x) for x in d]}
json.dump(out, open(OUT/"rev3_results.json", "w"), indent=1)
for site in SITES:
    print(site, {a: round(out["by_arm"][f"{site}|{a}"]["fpr"]["mean"], 3) for a in ARMS}, "| a2-a0:", round(out["paired"][f"{site}|a2_minus_a0ctrl"]["mean_diff"], 3), out["paired"][f"{site}|a2_minus_a0ctrl"]["n_positive"], out["paired"][f"{site}|a2_minus_a0ctrl"]["p_signflip"])

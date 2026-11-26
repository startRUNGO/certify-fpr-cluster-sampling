# REV-3 MH r1 (2026-11-16, WRITTEN, NOT RUN until 实验管理 approves): rev3_thrfix.py restricted to the rebuilt MP-IDB (site key mpidb_full) and the two arms of supplement table S8 (a2, a0-ctrl).  Only the MP-IDB scores and labels change (MH/<arm>_s<seed>_mh_logits.npz, MH/ids_mh.npz);
# the source threshold (INF/<stem>_eval.json) and the NIH validation fold scores (INF/<stem>_scores.npz nih_validation_logits) are the same frozen files as in rev3, because the classifiers are the same.  Output name rev3_mh_results.json.  Usage: python3 rev3_thrfix_mh.py <INF dir> <MH dir with a2 / a0-ctrl mh logits and ids_mh.npz> <out dir>
# THRFIX r2 (2026-11-08): ties at the threshold count as positive: the converted threshold is lowered by 1e-6 (> the float32 resolution of a logit near the source threshold, far below real score gaps), because the original compared float32 p >= t and a score that equals the Youden score up to float32 rounding must stay positive (a3_s43 on the NIH locked fold: 14 cells; thrfix r1 gave 92/1989 instead of the sealed 99/1996).
# THRFIX r1 (2026-11-03): the stored source threshold is a PROBABILITY-scale Youden threshold (softmax[:, 1]); the score s = logit1 - logit0 is now compared with logit(t) = ln(t/(1-t)). This is the only change to the original script.
"""REV-3 per PREREG_REV_FRONTEND_20261023.md (sha256 fcfb7f09...). Front-end design vs source-threshold FPR inflation; score-scale decomposition. CPU only, saved logits."""
import json, sys, glob, os, itertools, numpy as np
from pathlib import Path
INF = Path(sys.argv[1]); MH = Path(sys.argv[2]); OUT = Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
SITES = ["mpidb_full"]; SEEDS = [41, 42, 43, 44, 45]
ARMS = ["a0-ctrl", "a2"]
ids = np.load(INF/"ids_iter1.npz", allow_pickle=True); mid = np.load(MH/"ids_mh.npz", allow_pickle=True)   # ids: only the NIH validation fold labels are used; the rebuilt MP-IDB labels come from ids_mh.npz
def scores(stem, site):
    assert site == "mpidb_full"; z = np.load(MH/f"{stem}_mh_logits.npz")["logits"]; y = mid["y"]   # rebuilt MP-IDB (host-cell set), same stem naming as the s5/rev5 hc scripts
    return (z[:, 1]-z[:, 0]).astype(float), y
q25 = lambda x: np.percentile(x, 25); q75 = lambda x: np.percentile(x, 75)
rows = []
for arm in ARMS:
    for sd in SEEDS:
        stem = f"{arm}_s{sd}"; thr = json.load(open(INF/f"{stem}_eval.json"))["source_threshold"]["value"]; thr = float(np.log(thr / (1.0 - thr))) - 1e-6
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
    comps = {"a2_minus_a0ctrl": (["a2"], ["a0-ctrl"])}
    for nm, (g1, g2) in comps.items():
        d = np.array(per_seed(g1))-np.array(per_seed(g2)); out["paired"][f"{site}|{nm}"] = {"mean_diff": float(d.mean()), "sd_diff": float(d.std(ddof=1)), "n_positive": int((d > 0).sum()), "p_signflip": float(perm_p(d)), "diffs": [float(x) for x in d]}
json.dump(out, open(OUT/"rev3_mh_results.json", "w"), indent=1)
for site in SITES:
    print(site, {a: round(out["by_arm"][f"{site}|{a}"]["fpr"]["mean"], 3) for a in ARMS}, "| a2-a0:", round(out["paired"][f"{site}|a2_minus_a0ctrl"]["mean_diff"], 3), out["paired"][f"{site}|a2_minus_a0ctrl"]["n_positive"], out["paired"][f"{site}|a2_minus_a0ctrl"]["p_signflip"])

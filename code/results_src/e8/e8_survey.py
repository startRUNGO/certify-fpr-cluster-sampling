"""E8 PlasmoCount source-threshold FPR survey per PREREG_E8_PLASMOCOUNT_20261023.md (same definition as e7_analysis (a), sha256 3831fea9...)."""
import json, sys, glob, os, numpy as np
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
E8 = Path(sys.argv[1]); INF = Path(sys.argv[2]); OUT = Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
y = np.load(E8/"ids_plasmo.npz", allow_pickle=True)["y"]; rows = []
for f in sorted(glob.glob(str(INF/"*_eval.json"))):
    stem = os.path.basename(f)[:-10]; ev = json.load(open(f)); thr = ev["source_threshold"]["value"]; z = np.load(E8/f"{stem}_plasmo_logits.npz")["logits"]; s = z[:, 1]-z[:, 0]
    fpr = float((s[y == 0] >= thr).mean()); sens = float((s[y == 1] >= thr).mean()); sneg = np.sort(s[y == 0]); ref05 = float((s[y == 1] >= np.quantile(sneg, 0.95)).mean())
    rows.append({"ckpt": stem, "arm": stem.rsplit("_s", 1)[0], "auroc": float(roc_auc_score(y, s)), "fpr": fpr, "fpr_lb010": max(0.0, fpr-0.10), "sens": sens, "ba": .5*(sens+1-fpr), "sens_ref05": ref05})
rho, p = spearmanr([r["auroc"] for r in rows], [r["fpr"] for r in rows]); fp = np.array([r["fpr"] for r in rows]); lb = np.array([r["fpr_lb010"] for r in rows])
surv = {"spearman_rho": float(rho), "p": float(p), "n": len(rows), "fpr_min_med_max": [float(fp.min()), float(np.median(fp)), float(fp.max())], "frac_fpr_gt_0.10_eta0": float((fp > 0.10).mean()), "frac_fpr_gt_0.10_eta0.10_lb": float((lb > 0.10).mean()), "frac_fpr_gt_0.20_eta0": float((fp > 0.20).mean())}
json.dump({"rows": rows, "summary": surv}, open(OUT/"e8_plasmo_survey.json", "w"), indent=1); print(json.dumps(surv, indent=1))

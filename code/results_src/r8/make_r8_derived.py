"""Derived numbers for R8 (IML), S4 and S5. Reads only saved result files; no VM, no scores. Run from Paper5_StainRef/: python3 results_src/r8/make_r8_derived.py
Output: results_src/r8/r8_derived.json (flat dict id -> value)."""
import json, glob, numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; SRC = ROOT/"results_src"
out = {}
# AUROC of the 35 checkpoints on IML (meta files copied read-only from the VM run ~/P5_StainRef/runs/r8_iml/infer)
rows = [json.load(open(f)) for f in sorted(glob.glob(str(SRC/"r8/results/infer_meta/*_iml_meta.json")))]
for arm in sorted({r["ckpt"].rsplit("_s", 1)[0] for r in rows}):
    a = np.array([r["auroc"] for r in rows if r["ckpt"].rsplit("_s", 1)[0] == arm])
    for nm, v in (("n", len(a)), ("auroc_mean", a.mean()), ("auroc_min", a.min()), ("auroc_max", a.max())): out[f"iml|{arm}|{nm}"] = float(v) if nm != "n" else int(v)
a = np.array([r["auroc"] for r in rows])
for nm, v in (("n", len(a)), ("auroc_mean", a.mean()), ("auroc_min", a.min()), ("auroc_max", a.max())): out[f"iml|all|{nm}"] = float(v) if nm != "n" else int(v)
# oracle sensitivity on IML = S3 sens - ds_eps (ds_eps = sens - oracle_eps sensitivity)
R8 = json.load(open(SRC/"r8/results/r8_rev1_results.json"))["agg"]
for arm in ("a2", "a0-ctrl"):
    for eps in ("0.1", "0.05"):
        k = f"{arm}|iml|22|S3|{eps}"; out[f"iml|{arm}|oracle_sens_{eps}"] = R8[k]["sens"]-R8[k]["ds_eps"]; out[f"iml|{arm}|s3_sens_{eps}"] = R8[k]["sens"]
# development sites: S3 vs S4a/S4b (S4 run) and S5a/S5b (S5 run); minimum coverage and mean paired sensitivity difference over the three sites
S4 = json.load(open(SRC/"rev_s4/s4_results.json")); S5 = json.load(open(SRC/"rev_s5/s5_results.json")); SITES = ["bbbc041_all", "mpidb_full", "lacuna_B"]
for tag, R, rules in (("s4", S4, ["S3", "S4a", "S4b"]), ("s5", S5, ["S3", "S5a", "S5b"])):
    for arm in ("a2", "a0-ctrl"):
        for eps in ("0.1", "0.05"):
            for m in (22, 10):
                for r in rules:
                    cov = [R["agg"][f"{arm}|{s}|{m}|{r}|{eps}"]["cov"] for s in SITES]; i = f"dev|{arm}|{m}|{r}|{eps}|"
                    out[i+"cov_min"] = min(cov); out[i+"cov_argmin_site_index"] = int(np.argmin(cov))
                    if r != "S3":
                        d = [R["paired"][f"{arm}|{s}|{m}|{r}-S3|{eps}"]["sens_diff"] for s in SITES]; out[i+"sens_diff_mean"] = float(np.mean(d)); out[i+"sens_diff_min"] = min(d); out[i+"sens_diff_max"] = max(d)
                        c = [R["paired"][f"{arm}|{s}|{m}|{r}-S3|{eps}"]["cov_diff"] for s in SITES]; out[i+"cov_diff_mean"] = float(np.mean(c))
out["dev_site_order"] = "bbbc041_all, mpidb_full, lacuna_B (cov_argmin_site_index indexes this order)"
json.dump(out, open(SRC/"r8/r8_derived.json", "w"), indent=1); print(len(out), "derived values")

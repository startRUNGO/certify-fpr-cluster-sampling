"""R8 IML: self-check of rev1b vs rev1 and 90% intervals (seed level and draw level). Run on the Mac with Python >= 3.9; needs numpy and scipy. Usage: r8_intervals.py RESULTS_DIR"""
import json, sys, numpy as np
from pathlib import Path
from scipy.stats import t as tdist, beta
D = Path(sys.argv[1]); a = json.load(open(D/"r8_rev1_results.json")); b = json.load(open(D/"r8_rev1b_results.json")); P = json.load(open(D/"r8_rev1b_perdraw.json"))
bad = 0; n = 0
for sec in ("agg", "prop1"):
    assert a[sec].keys() == b[sec].keys()
    for k in a[sec]:
        for kk in a[sec][k]:
            n += 1
            if a[sec][k][kk] != b[sec][k][kk]: bad += 1
print(f"SELF-CHECK: compared {n} values, mismatches {bad}")
def wilson(k, N, z=1.6448536269514722):
    p = k/N; d = 1+z*z/N; c = (p+z*z/(2*N))/d; h = z*np.sqrt(p*(1-p)/N+z*z/(4*N*N))/d; return c-h, c+h
rows = []
for key in sorted(P):
    arm, site, m, name, eps = key.split("|"); eps = float(eps)
    if name not in ("S3", "R1_2"): continue
    e = P[key]; fpr = np.array(e["fpr"]); sens = np.array(e["sens"]); dse = np.array(e["dse"]); ns = len(fpr)//5
    out = {"key": key, "n": len(fpr)}
    for nm, v in (("sens", sens), ("dSens", dse)):
        sm = v.reshape(5, ns).mean(1); mean = float(v.mean()); se = sm.std(ddof=1)/np.sqrt(5); tq = tdist.ppf(0.95, 4)
        out[nm] = {"mean": mean, "seed_min": float(sm.min()), "seed_max": float(sm.max()), "seed_t90": [mean-tq*se, mean+tq*se], "draw_q05_q95": [float(np.quantile(v, .05)), float(np.quantile(v, .95))]}
    ok = int((fpr <= eps).sum()); out["cov"] = {"mean": ok/len(fpr), "wilson90_draw": list(wilson(ok, len(fpr))), "per_seed": [float((fpr.reshape(5, ns)[i] <= eps).mean()) for i in range(5)]}
    rows.append(out)
json.dump({"selfcheck_mismatches": bad, "compared": n, "rows": rows}, open(D/"r8_intervals.json", "w"), indent=1)
for r in rows:
    print(r["key"], "cov %.3f W90[%.3f,%.3f]" % (r["cov"]["mean"], *r["cov"]["wilson90_draw"]))
    for nm in ("sens", "dSens"):
        x = r[nm]; print("   %-6s mean %.3f | seed min-max [%.3f,%.3f] t90 [%.3f,%.3f] | draw 5-95%% [%.3f,%.3f]" % (nm, x["mean"], x["seed_min"], x["seed_max"], *x["seed_t90"], *x["draw_q05_q95"]))

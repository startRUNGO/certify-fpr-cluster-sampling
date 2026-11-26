# P5a D-e (实验管理 2026-11-25, decision 2 (a)): two-sided 90 % seed-cluster intervals for the confirmation-site forest plot F-e, computed LOCALLY from stored per-seed values.
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  post hoc descriptive interval; the preregistered criteria used the point values reported in Section 5.6.
# No rule is re-run, no image or score file is read, nothing is trained.  Only IML Malaria has stored per-seed values: g3_intervals.json (cd9543cf, run g3_intervals_20261030T135230Z),
#   cells iml|a2|{S3, R1_2}|{0.1, 0.05}|{cov, dsens} with per_seed (5 training seeds).  Interval: ebatch_lib.seed_cluster_ci (359b6bce) of the per-seed values.
#   PlasmoCount (all values) and the Delta Sens_eps of M5: per-draw values were not retained; no interval (the figure plots points only and says so).
#   M5 coverage keeps its stored 90 % seed-cluster interval (master m5/cells_p|SPLIT6|0.1|coverage_ci90_seedcluster|0/1); nothing is computed for it here.
# Anchors (exit 3, nothing written): input pins (g3, ebatch_lib, master 21496cc6); for every cell the g3 point equals the master value of its table_id within 1e-12
#   (r8_rev1_results/agg|a2|iml|22|<rule>|<eps>|cov / ds_eps), and the mean of per_seed equals the point within 1e-12; 5 seeds per cell.
# Output: NEW directory --out ('x' mode): d_e_confirm_intervals.json, SHA256SUMS.
# Run (local): python3 -E -B d_e_confirm_intervals.py --g3 <g3_intervals.json> --ebatch_lib <results_src/e_batch/ebatch_lib.py> --master <results_master.json> --out <new dir>
import argparse, hashlib, importlib.util, json, sys
from pathlib import Path
import numpy as np
PINS = {"g3": "cd9543cfc0cfe91b583b28431191f14858cda2b712f8f39b49b4539e33b58b75", "ebatch_lib": "359b6bce2dfdbf50a2acb8b909d631ab40d445447275915e2318e70981ef007f",
        "master": "21496cc647faa98ef17ee0261b7da65c766f5dbfe34aae64290396bf2e35124c"}
CELLS = [(r, e, k) for r in ("S3", "R1_2") for e in ("0.1", "0.05") for k in ("cov", "dsens")]
MASTER_FIELD = {"cov": "cov", "dsens": "ds_eps"}; TOL = 1e-12
NOTE = "post hoc descriptive interval; the preregistered criteria used the point values reported in Section 5.6"


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def main():
    ap = argparse.ArgumentParser()
    for k in ("g3", "ebatch_lib", "master", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); out = Path(a.out)
    if out.exists(): return fail("output directory exists", out)
    for k in PINS:
        if sha(getattr(a, k)) != PINS[k]: return fail("pin differs", k)
    spec = importlib.util.spec_from_file_location("ebatch_lib", a.ebatch_lib); EB = importlib.util.module_from_spec(spec); spec.loader.exec_module(EB)
    g3 = json.load(open(a.g3))["cells"]; V = json.load(open(a.master))["values"]; res = []; bad = []
    for r, e, k in CELLS:
        c = g3.get(f"iml|a2|{r}|{e}|{k}")
        if c is None: bad.append(("missing g3 cell", r, e, k)); continue
        mk = f"r8_rev1_results/agg|a2|iml|22|{r}|{e}|{MASTER_FIELD[k]}"; ps = c["per_seed"]
        if c["table_id"] != mk: bad.append(("table_id", r, e, k, c["table_id"]))
        if len(ps) != 5: bad.append(("seeds", r, e, k, len(ps)))
        if mk not in V or abs(c["point"] - V[mk]) > TOL: bad.append(("point vs master", r, e, k, c["point"], V.get(mk)))
        if abs(float(np.mean(ps)) - c["point"]) > TOL: bad.append(("mean(per_seed) vs point", r, e, k))
        res.append(dict(site="iml", rule=r, eps=float(e), stat=k, point=c["point"], per_seed=ps, ci90_seedcluster=EB.seed_cluster_ci(ps), master_key=mk))
    if bad: return fail("anchors", bad)
    rep = dict(nature="P5a D-e: " + NOTE, script_sha256=sha(__file__), pins=PINS, cells=res,
               no_interval=["PlasmoCount: all values (per-draw values were not retained; no interval)", "M5: Delta Sens_eps (per-draw values were not retained; no interval)"],
               stored_interval=["M5 coverage: m5/cells_p|SPLIT6|0.1|coverage_ci90_seedcluster|0/1 (stored 90% seed-cluster interval)"])
    payload = json.dumps(rep, indent=1, allow_nan=False); out.mkdir(parents=True)
    with open(out / "d_e_confirm_intervals.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh: fh.write(f"{sha(out / 'd_e_confirm_intervals.json')}  d_e_confirm_intervals.json\n")
    print(json.dumps([(x["rule"], x["eps"], x["stat"], x["point"], x["ci90_seedcluster"]) for x in res])); return 0


if __name__ == "__main__": sys.exit(main())

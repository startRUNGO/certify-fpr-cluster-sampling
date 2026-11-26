"""Synthetic smoke of the C5b evaluation scripts (no real data of any site is read). Builds a small synthetic world (60 fields per site, two arms x five seeds of random logits) in the layout of the frozen A2 files and in the C5b layout
(LOGDIR/<LAB>_s<seed>_<site>_logits.npz), runs the REVIEWED original scripts (unmodified except the intended deltas) and the seven c5b_* scripts with LAB = "a2" (A2's logits posing as C5b) and requires:
  * every a2 entry of the result JSON (agg and prop1) of the c5b script equals the original's, value for value;
  * the per-draw arrays (m = 22) equal the original's where the original writes them (lh, lb_same, lh_s11); for the others (BBBC, mh) the per-draw cov is consistent with the aggregated cov;
  * the output files carry the arm label; then the same scripts run with LAB = c5b-f on different logits and the expected output files exist (the list the wrapper checks).
--ndraw N replaces NDRAW = 200 by N in the ORIGINAL and the C5B copies alike (smoke speed only; the full 200 is used when N is 200). Usage: python3 c5b_smoke_eval.py --work <new dir> [--ndraw 20]"""
import argparse, json, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--work", required=True); ap.add_argument("--ndraw", type=int, default=20); a = ap.parse_args()
W = Path(a.work); assert not W.exists(); W.mkdir(parents=True)
HERE = Path(__file__).resolve().parent; RS = HERE.parent; OPT = RS.parent / "_opt_20261021"
SEEDS = [41, 42, 43, 44, 45]; ARMS = ["a2", "a0-ctrl"]
rng = np.random.default_rng(7)


def make_site(n_fields=60):
    y, field = [], []
    for f in range(n_fields):
        n = int(rng.integers(14, 27)); npos = int(rng.integers(0, 4)) if f % 5 else 0
        y += [1] * npos + [0] * (n - npos); field += ["F%03d" % f] * n
    return np.array(y), np.array(field)


SITES = {k: make_site() for k in ("bbbc041_all", "mpidb_full", "lacuna_B", "plasmo", "mh", "lh", "lb_same")}


def logits(y, shift):
    s = np.where(y == 1, rng.normal(1.0 + shift, 1.2, len(y)), rng.normal(-1.5 + shift, 1.0, len(y)))
    return np.stack([-s / 2, s / 2], 1).astype(np.float32)


A2 = {(arm, sd, k): logits(SITES[k][0], 0.1 * i) for i, arm in enumerate(ARMS) for sd in SEEDS for k in SITES}
CB = {(sd, k): logits(SITES[k][0], 0.3) for sd in SEEDS for k in SITES}
# ---- original layout: INF (ids_iter1, scores), E (everything else; passed as both E7 and E8)
INF = W / "orig_inf"; E = W / "orig_e"; INF.mkdir(); E.mkdir()
blob = {}
for k in ("bbbc041_all", "mpidb_full"): blob[k + "_y"], blob[k + "_group"] = SITES[k]
np.savez(INF / "ids_iter1.npz", **blob)
for arm in ARMS:
    for sd in SEEDS:
        np.savez(INF / f"{arm}_s{sd}_scores.npz", **{k + "_logits": A2[(arm, sd, k)] for k in ("bbbc041_all", "mpidb_full")})
        for k in ("lacuna_B", "plasmo", "mh", "lh", "lb_same"): np.savez(E / f"{arm}_s{sd}_{k}_logits.npz", logits=A2[(arm, sd, k)])
for k in ("lacuna_B", "plasmo", "mh", "lh", "lb_same"): np.savez(E / f"ids_{k}.npz", y=SITES[k][0], field=SITES[k][1])
# ---- C5b layout: LOGDIR/<LAB>_s<seed>_<site>_logits.npz
for lab in ("a2", "c5b-f", "c5b-fs"):
    d = W / ("uni_" + lab); d.mkdir()
    for sd in SEEDS:
        for k in ("bbbc041_all", "mh", "lh", "lb_same"):
            np.savez(d / f"{lab}_s{sd}_{k}_logits.npz", logits=(A2[("a2", sd, k)] if lab == "a2" else CB[(sd, k)] + (0.05 if lab == "c5b-fs" else 0)))
# ---- script copies
REF = W / "ref"; NEW = W / "new"; REF.mkdir(); NEW.mkdir()
shutil.copy(RS / "mh/orig_iter1/s8_bounds.py", REF / "s8_bounds.py"); shutil.copy(RS / "mh/orig_iter1/s8_bounds.py", NEW / "s8_bounds.py")
def patch(t): return t.replace("NDRAW = 200", "NDRAW = %d" % a.ndraw)
SRC = {"bbbc_rev1": RS / "rev1/rev1.py", "mh_rev1": RS / "mh/mh_rev1.py", "lh_rev1": RS / "paired/lh_paired_rev1.py", "lbsame_rev1": RS / "samepop/lb_same_paired_rev1.py",
       "bbbc_s11": OPT / "s11.py", "mh_s11": RS / "mh/mh_s11.py", "lh_s11": RS / "paired/lh_paired_s11.py"}
for k, p in SRC.items():
    t = p.read_text()
    if k == "lh_rev1":      # the intended delta: the OLD/NEW negative-score site check is not applicable (different models)
        i = t.index("    # ---- PAIRED site check"); j = t.index("    tasks = [(ai, sd, si)"); t = t[:i] + t[j:]
    assert "NDRAW = 200" in t; (REF / ("ref_" + k + ".py")).write_text(patch(t))
    t2 = (HERE / ("c5b_" + k + ".py")).read_text(); assert "NDRAW = 200" in t2; (NEW / ("c5b_" + k + ".py")).write_text(patch(t2))
ORIG_ARGS = {"bbbc_rev1": [INF, E, E], "mh_rev1": [INF, E, E], "lh_rev1": [INF, E, E], "lbsame_rev1": [INF, E, E], "bbbc_s11": [INF, E], "mh_s11": [INF, E], "lh_s11": [INF, E]}
SITE = {"bbbc_rev1": "bbbc041_all", "mh_rev1": "mh", "lh_rev1": "lh", "lbsame_rev1": "lb_same", "bbbc_s11": "bbbc041_all", "mh_s11": "mh", "lh_s11": "lh"}
ORIG_JSON = {"bbbc_rev1": "rev1_results.json", "mh_rev1": "mh_rev1_results.json", "lh_rev1": "lh_rev1_results.json", "lbsame_rev1": "lb_same_rev1_results.json", "bbbc_s11": "s11_results.json", "mh_s11": "mh_s11_results.json", "lh_s11": "lh_s11_results.json"}
NEW_JSON = {"bbbc_rev1": "c5b_bbbc_rev1_results_%s.json", "mh_rev1": "c5b_mh_rev1_results_%s.json", "lh_rev1": "c5b_lh_rev1_results_%s.json", "lbsame_rev1": "c5b_lbsame_rev1_results_%s.json", "bbbc_s11": "c5b_bbbc_s11_results_%s.json", "mh_s11": "c5b_mh_s11_results_%s.json", "lh_s11": "c5b_lh_s11_results_%s.json"}
ORIG_PD = {"lh_rev1", "lbsame_rev1", "lh_s11"}      # originals that write per-draw files
fails = []


def run(script, args, env=None):
    e = dict(os.environ); e.update(env or {}); e.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    r = subprocess.run([sys.executable, str(script)] + [str(x) for x in args], capture_output=True, text=True, env=e, cwd=str(script.parent))
    if r.returncode != 0: fails.append("%s exited %d: %s" % (script.name, r.returncode, r.stderr[-400:])); return False
    return True


def aggof(d): return d["agg"] if "agg" in d else d


for k in SRC:
    o = W / "out_ref" / k; n = W / "out_new_a2" / k
    ok1 = run(REF / ("ref_" + k + ".py"), ORIG_ARGS[k] + ([o] if k.endswith("s11") else [o]))
    ok2 = run(NEW / ("c5b_" + k + ".py"), ORIG_ARGS[k] + [n], {"C5B_ARM": "a2", "C5B_LOGDIR": str(W / "uni_a2")})
    if not (ok1 and ok2): continue
    jo = json.load(open(o / ORIG_JSON[k])); jn = json.load(open(n / (NEW_JSON[k] % "a2")))
    ao, an = aggof(jo), aggof(jn); site = SITE[k]
    ko = {x for x in ao if x.startswith("a2|%s|" % site)}; kn = {x for x in an}
    if ko != kn: fails.append("%s: key sets differ (%d vs %d)" % (k, len(ko), len(kn)))
    bad = [x for x in ko & kn if json.dumps(ao[x], sort_keys=True) != json.dumps(an[x], sort_keys=True)]
    if bad: fails.append("%s: %d agg entries differ, e.g. %s" % (k, len(bad), bad[0]))
    if "prop1" in jo:
        po = {x: v for x, v in jo["prop1"].items() if x.startswith("a2|%s" % site)}
        if json.dumps(po, sort_keys=True) != json.dumps(jn["prop1"], sort_keys=True): fails.append(k + ": prop1 differs")
    pdn = n / ("c5b_%s_perdraw_a2_m22.npz" % k); assert pdn.exists(), pdn
    zn = np.load(pdn, allow_pickle=True)
    if k in ORIG_PD:
        stem = "ref_" + k; zo = np.load(o / ("%s_perdraw_a2_m22.npz" % stem), allow_pickle=True)
        if set(zo.files) != set(zn.files) or any(not (zo[f].dtype == zn[f].dtype and np.array_equal(zo[f], zn[f])) for f in zo.files): fails.append(k + ": per-draw arrays differ from the original")
    # consistency of the per-draw output with the aggregated cov (both rev1 and s11 families)
    for f in zn.files:
        if f.endswith("|fpr"):
            name, eps = f.split("|")[0], float(f.split("|")[1]); key = "a2|%s|22|%s|%s" % (site, name, eps)
            if key in an and abs(float((zn[f] <= eps).mean()) - an[key]["cov"]) > 1e-12: fails.append("%s: per-draw cov inconsistent for %s" % (k, f)); break
    print("a2-as-c5b vs original:", k, "keys compared", len(ko & kn), "per-draw arrays", len(zn.files), flush=True)
# ---- other labels: file names and the wrapper's expected-output list
for lab, keys in (("c5b-f", list(SRC)), ("c5b-fs", ["lh_rev1", "bbbc_s11"])):
    for k in keys:
        n = W / ("out_new_" + lab) / k
        if run(NEW / ("c5b_" + k + ".py"), ORIG_ARGS[k] + [n], {"C5B_ARM": lab, "C5B_LOGDIR": str(W / ("uni_" + lab))}):
            for f in (NEW_JSON[k] % lab, "c5b_%s_perdraw_%s_m22.npz" % (k, lab)):
                if not (n / f).exists(): fails.append("missing output " + f)
tot = sum(f.stat().st_size for d in (W / "out_new_c5b-f").rglob("*") for f in [d] if f.is_file())
print("output bytes of one label (7 scripts, synthetic world, ndraw %d): %d" % (a.ndraw, tot))
print("RESULT:", "ALL PASSED" if not fails else "FAILED"); [print("  ", x) for x in fails]; sys.exit(1 if fails else 0)

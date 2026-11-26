"""Compares the eval_a2 results (the c5b scripts run on the frozen A2 logits) with the EXISTING A2 results (main-session instruction of 2026-10-27): for every pair in the pairs file, all entries of the new result JSON whose key starts with the pair's prefix
must equal the entries of the existing JSON (key sets equal, every numeric leaf within 1e-9, non-numeric leaves equal); per-draw npz pairs must be equal arrays. Prints PASS, or only the failing keys (pair label : key), and exits 3. No value is printed.
Pairs file (JSON list): [{"label": ..., "new": path, "ref": path, "prefix": "a2|mh|", "prop1_prefix": "a2|mh" or null}, {"label": ..., "new": npz, "ref": npz, "perdraw": true}]
Usage: python3 c5b_a2_compare.py --pairs pairs.json [--tol 1e-9]"""
import argparse, json, math, sys
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--pairs", required=True); ap.add_argument("--tol", type=float, default=1e-9); a = ap.parse_args(); bad = []


def leaf_equal(x, y):
    if isinstance(x, dict) or isinstance(y, dict):
        return isinstance(x, dict) and isinstance(y, dict) and set(x) == set(y) and all(leaf_equal(x[k], y[k]) for k in x)
    if isinstance(x, (list, tuple)) or isinstance(y, (list, tuple)):
        return isinstance(x, (list, tuple)) and isinstance(y, (list, tuple)) and len(x) == len(y) and all(leaf_equal(p, q) for p, q in zip(x, y))
    if isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool) and not isinstance(y, bool):
        return (math.isnan(x) and math.isnan(y)) or math.isclose(x, y, rel_tol=0.0, abs_tol=a.tol)
    return x == y


for p in json.load(open(a.pairs)):
    if p.get("perdraw"):
        zn, zr = np.load(p["new"], allow_pickle=True), np.load(p["ref"], allow_pickle=True)
        if set(zn.files) != set(zr.files): bad.append("%s:key_set" % p["label"]); continue
        for k in sorted(zr.files):
            if not (zn[k].dtype == zr[k].dtype and np.array_equal(zn[k], zr[k])): bad.append("%s:%s" % (p["label"], k))
        continue
    jn, jr = json.load(open(p["new"])), json.load(open(p["ref"]))
    an, ar = (jn["agg"] if "agg" in jn else jn), (jr["agg"] if "agg" in jr else jr)
    kn = {k for k in an if k.startswith(p["prefix"])}; kr = {k for k in ar if k.startswith(p["prefix"])}
    if kn != kr or not kr: bad.append("%s:key_set" % p["label"])
    for k in sorted(kn & kr):
        if not leaf_equal(an[k], ar[k]): bad.append("%s:%s" % (p["label"], k))
    if p.get("prop1_prefix"):
        pn = {k: v for k, v in jn.get("prop1", {}).items() if k.startswith(p["prop1_prefix"])}; pr = {k: v for k, v in jr.get("prop1", {}).items() if k.startswith(p["prop1_prefix"])}
        if not leaf_equal(pn, pr) or not pr: bad.append("%s:prop1" % p["label"])
if bad: print("\n".join(sorted(set(bad))[:60])); sys.exit(3)
print("PASS")

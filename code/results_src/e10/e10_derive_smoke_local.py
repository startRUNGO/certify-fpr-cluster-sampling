# Synthetic smoke for e10_derive.py (no project data).  Run: python3 -B e10_derive_smoke_local.py
import sys, json, math, tempfile
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import e10_derive as D, ebatch_lib as EB
fails = []
def ok(n, c, e=""):
    print(("PASS " if c else "FAIL ") + n + (" " + str(e) if e else ""), flush=True)
    if not c: fails.append(n)
def fake(b1_rates, b1u_rates, reps=5, nd=2000):
    cells = {}
    for ci, (b1, b1u) in enumerate(zip(b1_rates, b1u_rates)):
        ft = {"0.1": {}, "0.05": {}}
        for e in ft:
            for sd in range(41, 41 + reps):
                jitter = (sd - 43) * 0.002                                  # equal-weight mean over the realisations removes this exactly
                ft[e][str(sd)] = {k: {"site_fpr": 0.1, "n_draws": nd, "certify_rate": {"B1": b1 + jitter, "<EMAIL>": b1u + jitter, "<EMAIL>/3": b1u / 3 + jitter, "<EMAIL>": b1u / 2 + jitter}} for k in ("1.0", "0.5")}
        cells[str(ci)] = {"fixed_threshold": ft}
    return {"version": "e10_r1", "spec": "P5_规格_E10_合成受控模拟_20261103_draft3", "reps": reps, "cells": cells}
b1 = [0.108, 0.073, 0.100, 0.126, 0.103, 0.122, 0.109, 0.118, 0.115, 0.122]; b1u = [0.064, 0.046, 0.055, 0.088, 0.058, 0.073, 0.061, 0.07, 0.07, 0.078]
r = D.derive(fake(b1, b1u)); se = math.sqrt(0.1 * 0.9 / 10000)
ok("pooled rate equals the equal-weight mean over the realisations (jitter cancels)", all(abs(r["cells"][str(i)]["fixed_threshold_pooled"]["0.1"]["1.0"]["certify_rate"]["B1"] - b1[i]) < 1e-12 for i in range(10)))
ok("all four arms and both kappas are present for every cell", all(set(r["cells"][str(i)]["fixed_threshold_pooled"]["0.1"][k]["certify_rate"]) == set(D.ARMS) for i in range(10) for k in D.KAPPAS))
ok("pooled n, SE and limit: n = 10000, SE = sqrt(0.1 * 0.9 / 10000), limit = 0.1 + 3 SE", r["pooled_n"]["0.1"] == 10000 and abs(r["pooled_se"]["0.1"] - se) < 1e-15 and abs(r["pooled_limit"]["0.1"] - (0.1 + 3 * se)) < 1e-15 and abs(r["pooled_limit"]["0.1"] - 0.109) < 1e-12)
ok("count of B1 above the limit counts the five cells above 0.109 and not the cell at the limit", r["count_B1_above_limit"]["0.1"]["1.0"] == 5, r["count_B1_above_limit"])
ok("count of B1u above alpha is 0 here and counts an excess when present", r["count_B1u_above_alpha"]["0.1"]["1.0"] == 0 and D.derive(fake(b1, [0.064] * 9 + [0.101]))["count_B1u_above_alpha"]["0.1"]["1.0"] == 1)
ok("a rate exactly at alpha is not counted as above", D.derive(fake(b1, [0.10] * 10))["count_B1u_above_alpha"]["0.1"]["1.0"] == 0)
bad = fake(b1, b1u); bad["cells"]["3"]["fixed_threshold"]["0.1"]["42"]["1.0"]["n_draws"] = 1000
try: D.derive(bad); r_ = False
except AssertionError: r_ = True
ok("unequal numbers of draws between realisations are refused (the equal-weight mean would not be the pooled rate)", r_)
bad = fake(b1, b1u); del bad["cells"]["0"]["fixed_threshold"]["0.1"]["41"]["1.0"]["certify_rate"]["B1"]
try: D.derive(bad); r_ = False
except AssertionError: r_ = True
ok("a missing arm is refused", r_)
try: D.derive(fake(b1, b1u, reps=4)); r_ = False
except AssertionError: r_ = True
ok("a pooled n different from 10000 is refused", r_)
with tempfile.TemporaryDirectory() as td:
    td = Path(td); src = td / "src"; src.mkdir(); p = src / "e10.json"; json.dump(fake(b1, b1u), open(p, "w")); (src / "SHA256SUMS").write_text(f"{EB.sha(p)}  e10.json\n")
    sys.argv = ["e10_derive.py", "--e10", str(p), "--out", str(td / "new" / "e10_derived.json")]; D.main(); d = json.load(open(td / "new" / "e10_derived.json"))
    ok("end to end: derived json written with the source sha, keys as announced, SHA256SUMS present", d["source_sha256"] == EB.sha(p) and "cells" in d and d["count_B1_above_limit"]["0.1"]["1.0"] == 5 and (td / "new" / "SHA256SUMS").exists())
    flat = []
    def walk(o, pre):
        if isinstance(o, dict):
            for k, v in o.items(): walk(v, pre + [str(k)])
        else: flat.append("|".join(pre))
    walk(d, []); ok("the flattened keys have the announced form", "cells|3|fixed_threshold_pooled|0.1|1.0|certify_rate|B1" in flat and "pooled_se|0.1" in flat and "pooled_limit|0.1" in flat and "count_B1_above_limit|0.1|1.0" in flat and "count_B1u_above_alpha|0.1|1.0" in flat)
    (src / "SHA256SUMS").write_text("0" * 64 + "  e10.json\n"); sys.argv[-1] = str(td / "new2" / "e10_derived.json")
    try: D.main(); r_ = False
    except AssertionError: r_ = True
    ok("an e10.json that differs from its SHA256SUMS is refused", r_)
print("ALL PASS" if not fails else "FAILED: " + "; ".join(fails))

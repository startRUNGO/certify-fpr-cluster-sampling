# Local smoke test for e11_levelsplit.py.  Synthetic data only (no project file is read); the stored-REV-1 comparison (L.self_check) cannot run on synthetic data, so the end-to-end test replaces it by a function that returns [] and
# says so; everything else (draw loop, rules, A3, A4, summary, output files) is the real code.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Checks: (1) gate: a path containing plasmo / iml / r8 / m5 (any case) exits 4 in a subprocess BEFORE anything is read, and an already-loaded m5_* / e8_* / r8_* / e7_* module exits 4;
#         (2) fidelity: on synthetic inputs the S3, B1u_c0 and B1u_c1s rows of the copy equal m3_replay_lib.replay row by row (cert, cand, fpr, sens, dse, dso);
#         (3) A3 counts > 0 checked points and 0 violations, and a deliberately broken level (U3c above U3) is detected; A4 holds, and a doctored row is detected;
#         (4) band / near_boundary thresholds; paired differences; (5) end-to-end main(): output files exist, SHA256SUMS matches, summary has 3 sites x 2 eps x 4 rules, reading_flags present, an existing output directory is refused.
import json, subprocess, sys, tempfile, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; R = HERE.parent
sys.path[:0] = [str(HERE), str(R / "m3" / "k2prep")]
import m3_replay_lib as L
import e11_levelsplit as E
fails = []
def check(name, cond, info=""):
    print(("PASS " if cond else "FAIL ") + name + (" " + str(info) if (info and not cond) else ""))
    if not cond: fails.append(name)

# (1) gate
args = ["--a2link", "/x/a", "--ids_iter1", "/x/b", "--mhlh", "/x/c", "--rev1_bbbc", "/x/d", "--rev1_mh", "/x/e", "--rev1_lh", "/x/f", "--mh_rev1_py", "/x/g", "--out", "/x/out"]
for bad in ("<DATA>/plasmo_infer/x", "<DATA>/IML/x", "/runs/r8_iml/infer", "/runs/m5_infer/x"):
    a = list(args); a[1] = bad; rc = subprocess.run([sys.executable, "-I", "-B", str(HERE / "e11_levelsplit.py")] + a, capture_output=True).returncode
    check(f"gate exit 4 for {bad}", rc == 4, rc)
rc = subprocess.run([sys.executable, "-I", "-B", str(HERE / "e11_levelsplit.py")] + args, capture_output=True); check("clean paths pass the gate (fail later for another reason, not 4)", rc.returncode != 4, rc.returncode)
import types
sys.modules["m5_fake_module"] = types.ModuleType("m5_fake_module")
try: E.gate([]); check("loaded m5_* module exits 4", False)
except SystemExit as e: check("loaded m5_* module exits 4", e.code == 4, e.code)
del sys.modules["m5_fake_module"]

# synthetic site: 40 fields, ~30 negatives each, field shift, some positives
rng = np.random.default_rng(1); s, y, g = [], [], []
for f in range(40):
    mu = rng.normal(0, 0.4); n = 25 + int(rng.integers(0, 15)); s.append(mu + rng.normal(size=n)); y.append(np.zeros(n, int)); g.append(np.full(n, f))
    npos = int(rng.integers(0, 4)); s.append(mu + 2.5 + rng.normal(size=npos)); y.append(np.ones(npos, int)); g.append(np.full(npos, f))
S, Y, G = np.concatenate(s), np.concatenate(y), np.concatenate(g)
class Fake:
    used = {"synthetic": "none"}
    def load_site(self, sd, site):
        r = np.random.default_rng([7, sd]); return S + r.normal(0, 0.01, len(S)), Y, G
L.NDRAW = 12; L.SEEDS = [41, 42]
inp = Fake()

# (2) fidelity of the copy
a3 = dict(checked=0, violations=0, exceptions=0); mine = E.replay_e11(L, inp, 41, "bbbc041_all", a3); ref, _ = L.replay(inp, 41, "bbbc041_all")
keyf = lambda r: (r["draw"], r["eps"], r["rule"])
rm = {keyf(r): r for r in mine if r["rule"] in ("S3", "B1u_c0", "B1u_c1s")}; rr = {keyf(r): r for r in ref}
same = set(rm) == set(rr) and all(all(rm[k][f] == rr[k][f] for f in ("cert", "cand", "fpr", "sens", "dse", "dso")) for k in rm)
check("copy equals m3_replay_lib.replay for S3, B1u_c0, B1u_c1s", same, (len(rm), len(rr)))
check("rows drop ni and fields", all("ni" not in r and "fields" not in r for r in mine))
check("new rules present", {r["rule"] for r in mine} >= {"S3c", "B1u_c1s_c"})

# (3) A3 and A4
check("A3 checked > 0 and no violation", a3["checked"] > 0 and a3["violations"] == 0, a3)
old = E.LEVEL["U3c"]; E.LEVEL["U3c"] = 0.5; a3b = dict(checked=0, violations=0, exceptions=0); E.replay_e11(L, inp, 41, "bbbc041_all", a3b); E.LEVEL["U3c"] = old
check("A3 detects a broken level (U3c above U3)", a3b["violations"] > 0, a3b)
check("A4 holds on the replay", E.check_a4(mine) == [])
doc = [dict(r) for r in mine]; i = next((i for i, r in enumerate(doc) if r["rule"] == "S3" and not r["cert"]), None)
if i is None: doc = [dict(r) for r in mine]; i = next(i for i, r in enumerate(doc) if r["rule"] == "S3"); doc[i]["cert"] = False
for r in doc:
    if r["rule"] == "S3c" and (r["draw"], r["eps"]) == (doc[i]["draw"], doc[i]["eps"]): r["cert"] = True
check("A4 detects a doctored row (S3c certified while S3 is not)", len(E.check_a4(doc)) >= 1)
check("S3c never more certified than S3 (cannot_certify)", all(np.mean([r["cert"] for r in mine if r["rule"] == "S3c" and r["eps"] == e]) <= np.mean([r["cert"] for r in mine if r["rule"] == "S3" and r["eps"] == e]) + 1e-12 for e in (0.10, 0.05)))

# (4) bands and paired
check("band edges", E.band(0.90) == "hold_1-delta" and E.band(0.8999) == "below_1-delta_within_1-2delta" and E.band(0.80) == "below_1-delta_within_1-2delta" and E.band(0.7999) == "fail")
check("near boundary", E.near(0.905) and E.near(0.795) and not E.near(0.95) and not E.near(0.85))
rows2 = [dict(r, seed=41) for r in mine]
p = E.paired(rows2, "bbbc041_all", 0.10, "S3c", "S3"); check("paired has the keys", {"delta_coverage", "delta_ds_own", "delta_ds_own_per_seed", "seeds_with_negative_delta_ds_own"} <= set(p))
check("paired delta_ds_own equals the mean of the differences", abs(p["delta_ds_own"] - np.mean([a["dso"] - b["dso"] for a in rows2 for b in rows2 if a["rule"] == "S3c" and b["rule"] == "S3" and a["draw"] == b["draw"] and a["eps"] == 0.10 and b["eps"] == 0.10])) < 1e-12)

# (5) end-to-end with the test seam
tmp = Path(tempfile.mkdtemp()); outd = tmp / "out"
E.main(args[:-1] + [str(outd)], _hooks={"inputs": lambda L_, a_: Fake(), "self_check": lambda i_, r_: []})
check("output files exist", all((outd / n).exists() for n in ("e11_summary.json", "e11_rows.json", "SHA256SUMS")))
sums = {l.split()[1]: l.split()[0] for l in (outd / "SHA256SUMS").read_text().splitlines()}
check("SHA256SUMS matches", all(hashlib.sha256((outd / n).read_bytes()).hexdigest() == h for n, h in sums.items()))
rep = json.loads((outd / "e11_summary.json").read_text()); check("3 sites x 2 eps x 4 rules", len(rep["cells"]) == 24, len(rep["cells"]))
check("reading_flags and A3 recorded", "reading_flags" in rep and rep["anchors"]["A3"]["checked"] > 0)
check("every cell has an interval and a band", all("coverage_ci90_5seeds" in c and c["band"] in ("hold_1-delta", "below_1-delta_within_1-2delta", "fail") for c in rep["cells"].values()))
try: E.main(args[:-1] + [str(outd)], _hooks={"inputs": lambda L_, a_: Fake(), "self_check": lambda i_, r_: []}); check("existing output directory is refused", False)
except SystemExit as e: check("existing output directory is refused", True)
print("ALL PASS" if not fails else f"FAILED: {fails}"); sys.exit(0 if not fails else 1)

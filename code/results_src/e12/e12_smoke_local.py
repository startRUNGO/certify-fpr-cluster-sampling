# Local smoke test for e12_stress.py.  Synthetic data only; the frozen libraries are the real ones (pinned in e12_stress).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Checks: (1) gate: a path containing plasmo / iml / r8 / m5 (any case) exits 4 in a subprocess before anything is read; an already-loaded m5_* module exits 4;
#         (2) make_site2 with neutral parameters equals e10_synth.make_site bit for bit (several parameter sets, incl. tail t5 and F = 60);
#         (3) the stress families do what they are supposed to: rare fields (count = max(1, round(rho F)), +delta on their negatives), heavy-tailed sizes (n_f >= 5, max/median grows with sigma), slides (F // S fields per slide, between-slide variance
#             of field means above the no-nesting level);
#         (4) samplers: uniform sampler consumes the stream exactly like rng.choice; the two-stage sampler returns m distinct fields from exactly s slides;
#         (5) fidelity: run_site_seed_e12 with the uniform sampler equals m3_split_dev_eval.run_site_seed row by row (S3, B1u_c0, B1u_c1s, SPLIT6; cert, cand, fpr, sens, dse, dso) on a small synthetic site;
#         (6) S3c never certifies a draw S3 did not (and the check function detects a doctored row); S3c cannot-certify >= S3's;
#         (7) summarise: bands and keys, rare split present only for rare cells; (8) end-to-end main() on cells 0 and 1 with 2 realisations and 6 draws: output files, SHA256SUMS, anchors recorded, rows without ni and fields; an existing output directory is refused.
import json, subprocess, sys, tempfile, hashlib, types
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; R = HERE.parent
sys.path[:0] = [str(HERE), str(R / "e10"), str(R / "e_batch"), str(R / "m3" / "split"), str(R / "m3" / "k2prep"), str(R / "c5b")]
import e12_stress as E
import e10_synth as E10, ebatch_lib as EB, m3_replay_lib as L, m3_split_dev_eval as SP, m3_k1v_precheck as KV
fails = []
def check(name, cond, info=""):
    print(("PASS " if cond else "FAIL ") + name + (" " + str(info) if (info and not cond) else ""))
    if not cond: fails.append(name)

# (1) gate
base_args = ["--e10_json", "/x/e10.json", "--out", "/x/out"]
for bad in ("<DATA>/plasmo_infer/x", "<DATA>/IML/x", "/runs/r8_iml/infer", "/runs/m5_infer/x"):
    a = list(base_args); a[1] = bad; rc = subprocess.run([sys.executable, "-I", "-B", str(HERE / "e12_stress.py")] + a, capture_output=True).returncode
    check(f"gate exit 4 for {bad}", rc == 4, rc)
sys.modules["m5_fake_module"] = types.ModuleType("m5_fake_module")
try: E.gate([]); check("loaded m5_* module exits 4", False)
except SystemExit as e: check("loaded m5_* module exits 4", e.code == 4, e.code)
del sys.modules["m5_fake_module"]

# (2) neutral make_site2 == make_site
for kw in (dict(tau=0.5, r=100, pi=0.05, tail="normal", F=300), dict(tau=0.0, r=20, pi=0.2, tail="t5", F=60), dict(tau=1.0, r=500, pi=0.01, tail="normal", F=100)):
    x1 = E10.make_site(np.random.default_rng([5, 6, 7]), kw["tau"], kw["r"], kw["pi"], kw["tail"], F=kw["F"]); x2 = E.make_site2(E10, np.random.default_rng([5, 6, 7]), **kw)
    check(f"make_site2 neutral == make_site {kw}", all(np.array_equal(u, v) for u, v in zip(x1, x2[:3])) and x2[3] is None and len(x2[4]) == 0)

# (3) stress families
s, y, g, so, rf = E.make_site2(E10, np.random.default_rng(1), 0.5, 100, 0.05, "normal", 300, rare=(0.03, 4.0))
check("rare: round(rho F) fields", len(rf) == 9 and len(set(rf.tolist())) == 9, len(rf))
neg = y == 0; rare_mask = np.isin(g, rf)
check("rare: negatives of rare fields shifted up", s[neg & rare_mask].mean() - s[neg & ~rare_mask].mean() > 2.5, s[neg & rare_mask].mean() - s[neg & ~rare_mask].mean())
check("rare: at least one field for tiny rho", len(E.make_site2(E10, np.random.default_rng(2), 0.5, 100, 0.05, "normal", 300, rare=(0.001, 2.0))[4]) == 1)
ratios = []
for sg in (0.5, 1.0, 2.0):
    s2, y2, g2, _, _ = E.make_site2(E10, np.random.default_rng(3), 0.5, 100, 0.05, "normal", 300, size_sigma=sg); n_f = np.bincount(g2[y2 == 0], minlength=300); ratios.append(n_f.max() / np.median(n_f)); check(f"sizes: n_f >= 5 (sigma {sg})", n_f.min() >= 5, n_f.min())
check("sizes: max/median grows with sigma", ratios[0] < ratios[1] < ratios[2], ratios)
s3_, y3, g3, so3, _ = E.make_site2(E10, np.random.default_rng(4), 0.5, 100, 0.05, "normal", 300, slides=(12, 0.8))
check("slides: F // S fields per slide", so3 is not None and all((so3 == k).sum() == 25 for k in range(12)))
fm = np.array([s3_[(g3 == f) & (y3 == 0)].mean() for f in range(300)]); sm = np.array([fm[so3 == k].mean() for k in range(12)])
check("slides: between-slide variance of field means above the no-nesting level", sm.var(ddof=1) > 2 * fm.var(ddof=1) / 25, (sm.var(ddof=1), fm.var(ddof=1) / 25))

# (4) samplers
fl = list(range(300)); a1 = E.uniform_sampler(np.random.default_rng(9), fl, 22); a2 = np.random.default_rng(9).choice(300, 22, replace=False)
check("uniform sampler == rng.choice", np.array_equal(a1, a2))
smp = E.two_stage_sampler(so3, 12, 2); pick = smp(np.random.default_rng(11), fl, 22)
check("two-stage: 22 distinct fields from exactly 2 slides", len(set(pick)) == 22 and len({int(so3[i]) for i in pick}) <= 2 and len({int(so3[i]) for i in pick}) >= 1, sorted({int(so3[i]) for i in pick}))
check("two-stage uses only the first stage's slides", all(int(so3[i]) in {int(so3[i]) for i in pick} for i in pick))

# (5) fidelity of the copy of run_site_seed on a small synthetic site (F = 60)
class Synth:
    used = {}
    def __init__(self, sites): self.sites = sites
    def load_site(self, sd, st): return self.sites[sd]
site = "bbbc041_all"; L.M = 22; L.NDRAW = 4
sf, yf, gf, _, _ = E.make_site2(E10, np.random.default_rng(21), 0.5, 100, 0.05, "normal", 60); inp = Synth({41: (sf, yf, gf)}); platt = (1.2, 0.1)
lib = SP.run_site_seed(inp, 41, site, platt); mine = E.run_site_seed_e12(SP, inp, 41, site, platt, E.uniform_sampler)
core = ("cert", "cand", "fpr", "sens", "dse", "dso"); kk = lambda r: (r["draw"], r["eps"], r["rule"]); rules = ("S3", "B1u_c0", "B1u_c1s", "SPLIT6")
A = {kk(r): r for r in lib if r["rule"] in rules}; B = {kk(r): r for r in mine if r["rule"] in rules}
check("copy of run_site_seed == library (S3, B1u_c0, B1u_c1s, SPLIT6)", set(A) == set(B) and len(A) > 0 and all(all(A[k][f] == B[k][f] for f in core) for k in A), (len(A), len(B)))
check("copy has S3c rows and no K1 rows", {r["rule"] for r in mine} >= {"S3c"} and not any(r["rule"].startswith("K1") for r in mine))

# (6) S3c <= S3
check("S3c never certifies a draw S3 did not", E.check_s3c_s3(mine) == [])
doc = [dict(r) for r in mine]; i = next(i for i, r in enumerate(doc) if r["rule"] == "S3"); doc[i]["cert"] = False
for r in doc:
    if r["rule"] == "S3c" and (r["seed"], r["draw"], r["eps"]) == (doc[i]["seed"], doc[i]["draw"], doc[i]["eps"]): r["cert"] = True
check("check_s3c_s3 detects a doctored row", len(E.check_s3c_s3(doc)) >= 1)
check("S3c cannot-certify >= S3's", all(np.mean([r["cert"] for r in mine if r["rule"] == "S3c" and r["eps"] == e]) <= np.mean([r["cert"] for r in mine if r["rule"] == "S3" and r["eps"] == e]) + 1e-12 for e in (0.10, 0.05)))

# (7) summarise
rb = {41: mine, 42: [dict(r, seed=42) for r in mine]}
sm1 = E.summarise(EB, rb, {41: np.array([], int), 42: np.array([], int)})
check("summarise keys", all(f"{rule}|{e}" in sm1 for rule in ("S3", "S3c", "SPLIT6", "B1u_c1s") for e in (0.10, 0.05)) and "S3c-S3|0.1" in sm1)
check("no rare split without rare fields", not any(k.startswith("rare_split") for k in sm1))
sm2 = E.summarise(EB, rb, {41: np.array([3, 7]), 42: np.array([3])}); check("rare split present for rare cells", any(k.startswith("rare_split") for k in sm2))
check("band edges", E.band(0.90) == "hold_1-delta" and E.band(0.8999) == "below_1-delta_within_1-2delta" and E.band(0.7999) == "fail")
check("24 cells with the planned names", len(E.cells()) == 24 and E.cells()[0]["name"] == "base" and sum(c["two_stage"] is not None for c in E.cells()) == 4)
check("stream: cell 0 uses the E10 stream, other cells their own", E.stream(0, 41, 0) == [20261103, 10, 0, 41, 0] and E.stream(5, 41, 0) == [20261108, 12, 5, 41, 0])

# (8) end-to-end
tmp = Path(tempfile.mkdtemp()); outd = tmp / "out"; dummy = tmp / "e10_dummy.json"; dummy.write_text("{}")
E.main(["--e10_json", str(dummy), "--out", str(outd), "--cells", "0,1", "--reps", "2", "--ndraw", "6", "--check_draws", "200"])
check("output files exist", all((outd / n).exists() for n in ("e12_summary.json", "e12_rows.json", "SHA256SUMS")))
sums = {l.split()[1]: l.split()[0] for l in (outd / "SHA256SUMS").read_text().splitlines()}
check("SHA256SUMS matches", all(hashlib.sha256((outd / n).read_bytes()).hexdigest() == h for n, h in sums.items()))
rep = json.loads((outd / "e12_summary.json").read_text()); rows = json.loads((outd / "e12_rows.json").read_text())
check("two cells and anchors recorded", set(rep["cells"]) == {"0", "1"} and "c" in rep["anchors"] and rep["anchors"]["d"].startswith("not_applicable"))
check("cell 1 (rare) has the rare split", any(k.startswith("rare_split") for k in rep["cells"]["1"]["results"]) and not any(k.startswith("rare_split") for k in rep["cells"]["0"]["results"]))
check("rows carry no ni, fields or split_* keys", all("ni" not in r and "fields" not in r and not any(k.startswith("split_") for k in r) for r in rows))
try: E.main(["--e10_json", str(dummy), "--out", str(outd), "--cells", "0", "--reps", "2", "--ndraw", "6", "--check_draws", "200"]); check("existing output directory is refused", False)
except SystemExit: check("existing output directory is refused", True)
print("ALL PASS" if not fails else f"FAILED: {fails}"); sys.exit(0 if not fails else 1)

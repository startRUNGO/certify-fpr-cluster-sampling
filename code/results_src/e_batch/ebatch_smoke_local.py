# Synthetic smoke test for the E-batch scripts (no project data, no GPU, CPU seconds).  Run:  python3 -I -B ebatch_smoke_local.py   -> prints ALL PASS or raises.
import json, math, sys, tempfile
from pathlib import Path
import numpy as np
from scipy.stats import binom
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import ebatch_lib as E, e1_baselines as E1, e2n_resample as E2, e3_unit as E3, e5_clinical as E5

ok = lambda name: print("pass:", name)


def synth(seed=0, nf=120, per=(40, 80), shift_sd=0.0, npos_per=10):
    r = np.random.default_rng(seed); s, y, g = [], [], []
    for f in range(nf):
        k = int(r.integers(*per)); off = r.normal(0, shift_sd); s.append(r.normal(off, 1, k)); y.append(np.zeros(k, int)); g.append(np.full(k, f))
        s.append(r.normal(off + 2.5, 1, npos_per)); y.append(np.ones(npos_per, int)); g.append(np.full(npos_per, f))
    return np.concatenate(s), np.concatenate(y), np.concatenate(g)


# 1. rule definitions against independent implementations
for n in (5, 10, 22, 40, 100, 200):
    for eps in (0.10, 0.05):
        xs = [x for x in range(n) if binom.cdf(x, n, eps) <= E.DELTA]          # CP upper limit <= eps  <=>  P_eps(X <= x) <= delta
        assert E.x_star(n, eps) == (max(xs) if xs else -1), ("x*", n, eps)
assert E.k_mc(19, 0.10) == 18 and E.k_mc(20, 0.10) == 19 and E.k_mc(22, 0.05) == 22 and E.k_q(20, 0.10) == 18 and E.k_q(22, 0.10) == 20; ok("order statistics and float guard")
assert E.x_star(22, 0.10) == 0 and E.x_star(22, 0.05) == -1 and E.x_star(40, 0.05) == -1 and E.x_star(40, 0.10) >= 1; ok("n_min facts of the spec (m=22 only zero-failure at 0.10; 0.05 vacuous at 22 and 40)")
neg = np.sort(np.random.default_rng(1).normal(size=22)); t = E.tau_cp(neg, 0.10); assert t > neg[-1] and (neg >= t).sum() == 0
assert math.isinf(E.tau_cp(neg, 0.05)) and math.isinf(E.tau_mc(neg[:5], 0.10)) and not math.isinf(E.tau_q(neg[:5], 0.10)); ok("tau cases")

# 2. exact analytic law against simulation (iid, large N)
s, y, g = synth(3, nf=300, per=(60, 100)); pop = E.make_pop(s, y, g); assert pop["N"] > 15000
rows = []
for n in (20, 50, 100):
    recs = E2.e2n_cell(pop, 0, 41, 0, n, "iid", 3000)
    for e in E.EPSS:
        for r in E.RULES: rows.append((r, n, e, float(np.mean([x["viol"] for x in recs[e][r]]))))
E.check_analytic(rows, 3000, 2.0 / pop["N"]); ok("iid violation = exact Beta law for Q, MC, CP (n=20,50,100)")
badrows = [(r, n, e, o + (0.05 if (r, n, e) == ("MC", 50, 0.10) else 0)) for r, n, e, o in rows]
try: E.check_analytic(badrows, 3000, 2.0 / pop["N"]); raise SystemExit("negative control did not fire")
except E.CheckError: ok("analytic check fires on a wrong rule")
assert abs(E.viol_analytic("CP", 22, 0.10) - 0.9 ** 22) < 1e-12 and E.viol_analytic("CP", 22, 0.05) is None; ok("CP at n=22 is 0.9^22 <= delta")

# 3. field mode truncates to n, keeps clusters; clustered fields break CP, iid does not
r = np.random.default_rng(5); sizes = []
for i in range(50):
    v = E2.draw_field(pop, r, 37); sizes.append(len(v)); assert np.all(np.diff(v) >= 0)
assert set(sizes) == {37}; assert len(set(E2.draw_iid(pop, np.random.default_rng(7), 200).tolist())) == 200; ok("field draw has exactly n cells, iid draw has no repeats")
sc, yc, gc = synth(4, nf=150, per=(8, 14), shift_sd=1.2); pc = E.make_pop(sc, yc, gc)
vi = np.mean([x["viol"] for x in E2.e2n_cell(pc, 1, 41, 2, 50, "iid", 600)[0.10]["CP"]]); vf = np.mean([x["viol"] for x in E2.e2n_cell(pc, 1, 41, 2, 50, "field", 600)[0.10]["CP"]])
assert vi <= E.DELTA + 0.03 and vf > vi + 0.05, (vi, vf); ok("cluster structure inflates CP violation (field %.3f vs iid %.3f)" % (vf, vi))
assert E2.e2n_cell(pc, 1, 41, 4, 10 ** 7, "iid", 5) is None; ok("n larger than the site is infeasible, not extrapolated")

# 4. E1 cell: reproducible stream, N0 fixed, coverage fields, vacuous bookkeeping
a1, n1, tq1 = E1.e1_cell(pop, 0, 41, 0, 0.10, 0.5, ndraw=40); a2, _, _ = E1.e1_cell(pop, 0, 41, 0, 0.10, 0.5, ndraw=40)
assert [x["fpr"] for x in a1["MC"]] == [x["fpr"] for x in a2["MC"]] and len({x["fpr"] for x in a1["N0"]}) == 1; ok("E1 stream reproducible, N0 draw-independent")
sm = E.summarise(a1["MC"]); assert 0 <= sm["coverage"] <= 1 and sm["vacuous_share"] == 0 and sm["fpr_p90"] >= sm["fpr_median"]
small = E.make_pop(*synth(6, nf=30, per=(3, 4))); b1, nn, _ = E1.e1_cell(small, 0, 41, 1, 0.05, 0.0, ndraw=20)
assert all(x["vac"] for x in b1["MC"]) or min(nn) >= 20, "MC must be vacuous when (n+1)(1-eps) > n"; ok("MC vacuous when the pooled negatives are few")
try: E1.e1_cell(E.make_pop(*synth(8, nf=10)), 0, 41, 0, 0.1, 0.0, ndraw=2); raise SystemExit("m check missing")
except AssertionError: ok("fewer fields than m refused")

# 5. M1 anchor (r4) on synthetic NIH-validation scores: (a) n, (b) AUROC, (c) Youden threshold on the probability scale
r = np.random.default_rng(11); yv = np.r_[np.zeros(800, int), np.ones(400, int)]; sv = np.r_[r.normal(-1, 1, 800), r.normal(1.5, 1, 400)]
tre, Jv = E1.youden_prob(yv, 1 / (1 + np.exp(-sv))); tstored = float(np.float32(tre))             # the stored value went through float32
def write(td, n=None, auc=None, val=None):
    p = Path(td) / "ev.json"; json.dump({"ckpt": "a2_s41", "config": "a2ref", "source_threshold": {"from": "nih_validation", "rule": "Youden J", "value": tstored if val is None else val},
                                         "nih_validation": {"auroc": E.auroc(sv, yv) if auc is None else auc, "n": len(sv) if n is None else n}}, open(p, "w")); return p
with tempfile.TemporaryDirectory() as td:
    tau, tv, info = E1.nih_val_anchor(write(td), sv, yv); assert abs(tau - math.log(tv / (1 - tv))) < 1e-15 and info["n"] == len(sv); ok("M1 anchor passes on consistent values (stored threshold through float32)")
    for name, kw in (("n", dict(n=len(sv) + 1)), ("auroc", dict(auc=E.auroc(sv, yv) + 1e-6)), ("threshold", dict(val=tstored + 5e-4))):
        try: E1.nih_val_anchor(write(td, **kw), sv, yv); raise SystemExit("anchor did not fire: " + name)
        except E.AnchorError: ok("M1 anchor fires when " + name + " is off")
# the Youden helper against a brute-force J maximum, and the first-maximiser rule on a tie
yy = np.array([0, 0, 1, 1, 0, 1]); pp = np.array([0.1, 0.4, 0.35, 0.8, 0.6, 0.9]); tb, Jb = E1.youden_prob(yy, pp)
best = max(((np.mean(pp[yy == 1] >= c) - np.mean(pp[yy == 0] >= c)), c) for c in pp); assert abs(Jb - best[0]) < 1e-12 and tb == best[1]; ok("Youden helper equals the brute-force maximum of TPR - FPR")
tie = E1.youden_prob(np.array([1, 1, 0, 0]), np.array([0.9, 0.8, 0.3, 0.2])); assert tie[0] == 0.8 and tie[1] == 1.0; ok("Youden first maximiser on a flat J stretch is the highest threshold")

# 5c. r5: AUROC anchor on the score scale of the stored value
rg = np.random.default_rng(21); bases = []
for sd in range(5):
    yv5 = np.r_[np.zeros(300, int), np.ones(200, int)]; sv5 = np.r_[rg.normal(-2, 2, 300), rg.normal(3, 2, 200)]; sv5[:6] = 30.0 + np.arange(6); sv5[300:306] = 30.0 + np.arange(6) + 0.5      # saturated float32 ties across classes
    bases.append(E.auroc_bases(sv5, yv5))
assert all(b["s"] != b["p32"] for b in bases); ms, mp = float(np.mean([b["s"] for b in bases])), float(np.mean([b["p32"] for b in bases]))
site = E.SITES[0]; key = f"survey_arm/{site}|a2|auroc_mean"
assert E.anchor_auroc({key: mp}, site, bases) == "p32"; assert E.anchor_auroc({key: ms}, site, bases) == "s"; ok("AUROC anchor passes on the p32 scale and on the s scale, and reports which one matched")
assert E.anchor_auroc({key: 0.7}, site, [dict(s=0.7, p32=0.7)] * 5) == "both"; ok("AUROC anchor reports both when the scales agree")
for bad_ in (ms + 1e-9, mp - 1e-6, 0.5):
    try: E.anchor_auroc({key: bad_}, site, bases); raise SystemExit("anchor did not fire")
    except E.AnchorError: pass
ok("AUROC anchor fires when neither scale matches (1e-9 off is not accepted)")
try: E.anchor_auroc({}, site, bases); raise SystemExit("missing key not refused")
except E.AnchorError: ok("AUROC anchor fires on a missing master key")
p_ = E.p32(np.array([-100.0, 0.0, 17.0, 30.0])); assert p_.dtype == np.float32 and p_[2] == 1.0 and p_[3] == 1.0 and p_[0] == 0.0 and abs(p_[1] - 0.5) < 1e-7; ok("p32 is float32 and saturates where the survey scores saturated")

# 6. E3: one cell per field, m cells; spec facts; infeasible handling
rr = E3.e3_cell(pop, 0, 41, 0, 0.10, 22, ndraw=300); sm = E.summarise(rr); assert sm["vacuous_share"] == 0 or sm["vacuous_share"] < 1
assert all(x["vac"] for x in E3.e3_cell(pop, 0, 41, 1, 0.05, 22, ndraw=20)) and all(x["vac"] for x in E3.e3_cell(pop, 0, 41, 1, 0.05, 40, ndraw=20)); ok("E3: eps 0.05 vacuous at m=22 and 40")
assert E3.e3_cell(E.make_pop(*synth(9, nf=30)), 0, 41, 0, 0.10, 40, ndraw=2) is None; ok("E3: site with fewer than m fields is infeasible")
assert sm["coverage"] >= 0.88, sm; ok("E3 R1-1 coverage on exchangeable synthetic data %.3f (>= 0.88, delta 0.10 + MC noise)" % sm["coverage"])

# 5b. F2/F3/F7/F1: non-vacuous convention, per-draw thresholds, interval definition, pins
sm2 = E.summarise(b1["MC"]); assert sm2["coverage_nonvacuous"] is None or sm2["vacuous_share"] < 1; assert E.summarise(b1["MC"])["convention"] == E.VAC_CONVENTION
mix = [dict(vac=True, fpr=0.0, viol=False, sens=0.0, dso=-0.5), dict(vac=False, fpr=0.2, viol=True, sens=0.8, dso=0.0), dict(vac=False, fpr=0.05, viol=False, sens=0.9, dso=0.1)]
sx = E.summarise(mix); assert abs(sx["coverage"] - 2 / 3) < 1e-12 and abs(sx["coverage_nonvacuous"] - 0.5) < 1e-12 and abs(sx["fpr_mean_nonvacuous"] - 0.125) < 1e-12 and abs(sx["fpr_mean"] - 0.25 / 3) < 1e-12; ok("vacuous and non-vacuous conventions differ as stated")
assert E.summarise([mix[0]] * 3)["coverage_nonvacuous"] is None and E1.cov_nv(E.summarise([mix[0]] * 3)) == 1.0; ok("all-vacuous cell has no non-vacuous coverage")
assert len(tq1["Q"]) == len(tq1["MC"]) == 40 and all(isinstance(v, float) for v in tq1["Q"]); ok("E1 returns one threshold per draw")
ci = E.seed_cluster_ci([0.9, 0.92, 0.94, 0.88, 0.91]); assert ci == E.seed_cluster_ci([0.9, 0.92, 0.94, 0.88, 0.91]) and ci[0] < np.mean([0.9, 0.92, 0.94, 0.88, 0.91]) < ci[1] and "B=10000" in E.CI_DEFINITION; ok("seed-cluster interval deterministic, B=10000")
assert len(E1.K_PIN) == 64 and E1.K_PIN.startswith("0e6aa9f8") and E1.SPLIT6_KEY.format(site="mpidb_full", eps=0.1) == "ssrcps/cells|SPLIT6|mpidb_full|0.1|coverage" and E1.SPLIT6_KEY.format(site="lacuna_B", eps=0.05).endswith("|0.05|coverage"); ok("K pin and SPLIT6 key format")
fake = {s: {e: {r: dict(coverage=0.5, coverage_nonvacuous=0.5, vacuous_share=0.0) for r in ("N0", "Q", "MC")} for e in E.EPSS} for s in E.SITES}
sc = {s: {e: {"S3": 0.95, "R1_all": 0.60, "SPLIT6": 0.95} for e in E.EPSS} for s in E.SITES}
rd = E1.reading(fake, sc); assert rd["E1_a"] is True and rd["R1_all_stored"]["mpidb_full"]["0.1"] == 0.60; ok("E1-a still fires when R1_all is 0.60 (not a condition); R1_all reported unchanged")
sc["bbbc041_all"][0.10]["SPLIT6"] = None; sc["bbbc041_all"][0.05]["S3"] = 0.80; rd = E1.reading(fake, sc)
assert rd["E1_a"] is False; ok("E1-a does not fire when S3 or SPLIT6 is missing or below 0.90 on a site where Q or MC is below")
sc2 = {s: {e: {"S3": 0.95, "R1_all": 0.95, "SPLIT6": 0.95} for e in E.EPSS} for s in E.SITES}; ok_cov = {s: {e: {r: dict(coverage=0.95, coverage_nonvacuous=0.95, vacuous_share=0.0) for r in ("N0", "Q", "MC")} for e in E.EPSS} for s in E.SITES}
rd = E1.reading(ok_cov, sc2); assert rd["E1_a"] is False and rd["E1_b"] is True; ok("E1-b when Q and MC cover at every site")
with tempfile.TemporaryDirectory() as td:
    o = Path(td) / "w" / "e.json"; E.write_out(o, {"v": 1}, {"a__41__0.1__tau_Q": np.array([0.5, np.inf])}); z = np.load(o.parent / "e_taus.npz")
    assert np.isinf(z["a__41__0.1__tau_Q"][1]) and len((o.parent / "SHA256SUMS").read_text().strip().splitlines()) == 2; ok("per-draw thresholds written with 'x' and listed in SHA256SUMS")

# 7. E5 arithmetic, writer, path refusal, no locked-fold access in sources
c = E5.convert(0.02, 50.0); assert abs(c["fp_per_100_fields"] - 100.0) < 1e-12 and abs(c["p_any_fp"]["500"] - (1 - 0.98 ** 500)) < 1e-12; ok("E5 conversion")
with tempfile.TemporaryDirectory() as td:
    o = Path(td) / "new" / "x.json"; E.write_out(o, {"a": float("inf"), "b": np.float64(1.5)}); assert json.load(open(o)) == {"a": None, "b": 1.5} and (o.parent / "SHA256SUMS").exists()
    try: E.write_out(o, {}); raise SystemExit("overwrite allowed")
    except FileExistsError: ok("writer: new directory only, inf -> null")
try: E.refuse_path("/x/plasmo/y"); raise SystemExit("not refused")
except AssertionError: ok("confirmatory path refused")
for f in ("ebatch_lib.py", "e1_baselines.py", "e2n_resample.py", "e3_unit.py", "e5_clinical.py"):
    src = (HERE / f).read_text()
    for bad in ("np.load(", "torch", "PIL", "cv2", "cuda"): assert bad not in src, (f, bad)
assert "locked_test" not in (HERE / "e1_baselines.py").read_text().replace("locked-test fold", "") ; ok("sources: no np.load, no torch/PIL/cv2/cuda, no locked-test key")
print("ALL PASS")

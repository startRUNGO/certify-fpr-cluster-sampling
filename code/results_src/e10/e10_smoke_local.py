# Synthetic smoke for e10_synth.py (everything is synthetic anyway; this checks the generator, the analytic ICC, the fixed-threshold test and a tiny end-to-end run).  LOCAL; CPU.  Run: python3 -B e10_smoke_local.py
import sys, json, tempfile
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import e10_synth as S, ebatch_lib as EB
import m3_replay_lib as L, s8_bounds as B8
fails = []
def ok(n, c, e=""):
    print(("PASS " if c else "FAIL ") + n + (" " + str(e) if e else ""), flush=True)
    if not c: fails.append(n)
g_ = S.cells(False); ok("one-factor-at-a-time design has 10 cells including the base; the factorial has 96", len(g_) == 10 and g_[0] == S.BASE and len(S.cells(True)) == 96 and len({tuple(sorted(c.items())) for c in S.cells(True)}) == 96)
rng = np.random.default_rng(1); s, y, g = S.make_site(rng, 0.5, 100, 0.05, "normal"); pos = y.mean()
ok("generator: 300 fields, positive share near pi, negative scores near N(0, 1 + tau^2)", len(np.unique(g)) == 300 and abs(pos - 0.05) < 0.01 and abs(s[y == 0].std() - np.sqrt(1 + 0.25)) < 0.05, (round(pos, 4), round(s[y == 0].std(), 3)))
st, yt, _ = S.make_site(np.random.default_rng(2), 0.0, 200, 0.05, "t5"); sn, yn, _ = S.make_site(np.random.default_rng(2), 0.0, 200, 0.05, "normal"); st, sn = st[yt == 0], sn[yn == 0]
ok("t5 noise is rescaled to unit variance but has heavier tails than the normal (negatives)", abs(st.std() - 1) < 0.1 and (np.abs(st) > 3.5).mean() > 2 * (np.abs(sn) > 3.5).mean(), (round(st.std(), 3), round(float((np.abs(st) > 3.5).mean()), 5), round(float((np.abs(sn) > 3.5).mean()), 5)))
ic = [S.realised_icc(*S.make_site(np.random.default_rng([3, i]), 1.0, 100, 0.05, "normal")) for i in range(6)]; ok("realised ICC at tau = 1 matches the analytic value (quadrature)", abs(np.mean(ic) - S.analytic_icc(1.0)) < 0.04, (round(float(np.mean(ic)), 3), round(S.analytic_icc(1.0), 3)))
ok("analytic ICC increases with tau and is 0 at tau = 0", S.analytic_icc(0.0) < 1e-6 and S.analytic_icc(0.25) < S.analytic_icc(0.5) < S.analytic_icc(1.0))
for tau in (0.0, 1.0):
    s, y, g = S.make_site(np.random.default_rng([4, int(tau * 10)]), tau, 100, 0.05, "normal")
    r = S.fixed_threshold_test(L, B8, s, y, g, 0.10, np.random.default_rng(5), 800)
    t1 = r["1.0"]["certify_rate"]; pw = r["0.5"]["certify_rate"]
    ok(f"fixed-threshold test at tau = {tau}: boundary FPR >= eps, type-I of <EMAIL> <= alpha + 3 SE when the cells are iid, power above type-I", r["1.0"]["site_fpr"] >= 0.10 and (tau > 0 or t1["<EMAIL>"] <= 0.10 + 3 * np.sqrt(0.09 / 800)) and pw["<EMAIL>"] >= t1["<EMAIL>"], (round(r["1.0"]["site_fpr"], 4), {k: round(v, 3) for k, v in t1.items()}))
    if tau == 1.0: ok("under strong clustering the point-ICC bound B1 has a larger type-I rate than the ICC-upper-limit bound B1u (the reason for B1u)", t1["B1"] > t1["<EMAIL>"], {k: round(v, 3) for k, v in t1.items()})
ok("the split-level diagnostic (alpha / 2 each) is never less conservative than the frozen alpha-each form on random inputs", all(S.b1u_split(x, n, 0.10) >= L.b1u(x, n, 0.10) - 1e-12 for x, n in [(np.random.default_rng(i).binomial(40, 0.1, 22).astype(float), np.full(22, 40.0)) for i in range(20)]))
# r2: the iid type-I anchor of generator_checks
rr = S.generator_checks(np.random.default_rng([20261103, 10, 9999, 0, 0]), L, B8, ndraw=1000); it = rr["iid_type1_B1u_0.10"]
ok("generator anchor: iid type-I of <EMAIL> at the boundary is within alpha + 3 SE and the boundary FPR is >= eps", it["rate"] <= it["limit"] and it["site_fpr"] >= 0.10 and abs(it["limit"] - (0.10 + 3 * (0.09 / 1000) ** 0.5)) < 1e-12, {k: round(v, 4) for k, v in it.items() if k != "n_draws"})
orig_b1u = L.b1u
L.b1u = lambda x, n, delta: 0.0                       # a bound that always certifies: type-I rate 1
try: S.generator_checks(np.random.default_rng([20261103, 10, 9999, 0, 0]), L, B8, ndraw=200); r_ = False
except EB.AnchorError: r_ = True
L.b1u = orig_b1u; ok("generator anchor fires (exit 3 path) when the bound certifies at the boundary too often", r_)
# r3: Platt with Newton fallback, source set size
sys.path.insert(0, str(HERE.parent / "m3" / "k1vb"))
import m3_k1v_precheck as KV
rg = np.random.default_rng(8); ys_ = (rg.random(3000) < 0.3).astype(int); ss_ = rg.normal(0, 1, 3000) + 1.5 * ys_
(a1, b1), i1 = S.platt_synth(KV, ss_, ys_); (a2, b2) = KV.platt(ss_, ys_); ok("platt_synth uses the frozen fit when it succeeds (no fallback)", (not i1["fallback"]) and abs(a1 - a2) < 1e-12 and abs(b1 - b2) < 1e-12)
def boom(s, y): raise AssertionError(("platt fit", "Desired error not necessarily achieved due to precision loss."))
orig_platt = KV.platt; KV.platt = boom; (a3, b3), i3 = S.platt_synth(KV, ss_, ys_); KV.platt = orig_platt
g3 = S._platt_grad(ss_.astype(float), ys_.astype(float), np.array([a3, b3]))
ok("on the frozen assertion the Newton fit is used, reaches a vanishing gradient and agrees with the BFGS optimum", i3["fallback"] and np.linalg.norm(g3) < 1e-8 and abs(a3 - a2) < 1e-3 and abs(b3 - b2) < 1e-3 and i3["nll"] <= i3["bfgs_nll_at_end"] + 1e-9, (round(a3, 4), round(a2, 4)))
def other(s, y): raise AssertionError(("something else", 1))
KV.platt = other
try: S.platt_synth(KV, ss_, ys_); r_ = False
except AssertionError: r_ = True
KV.platt = orig_platt; ok("a different assertion is not swallowed by the fallback", r_)
sep_s = np.r_[rg.normal(0, 1, 200), rg.normal(6, 0.3, 200)]; sep_y = np.r_[np.zeros(200, int), np.ones(200, int)]; sep_y[0] = 1; sep_y[-1] = 0     # nearly separable: still finite
(a4, b4), i4 = S.platt_synth(KV, sep_s, sep_y); ok("a nearly separable source set gives finite Platt parameters", np.isfinite([a4, b4]).all())
for r_, pi_ in ((20, 0.01), (100, 0.05), (500, 0.4)):
    ss2, ys2, gs2 = S.source_set(np.random.default_rng(9), r_, pi_, "normal"); ok(f"source set at r = {r_}, pi = {pi_}: about 4000 negative cells and enough positives for a Platt fit", 3500 < int((ys2 == 0).sum()) < 4700 and ys2.sum() >= 20 and len(np.unique(gs2)) == int(np.ceil(4000 / r_)), (int((ys2 == 0).sum()), int(ys2.sum())))
with tempfile.TemporaryDirectory() as td:
    out = Path(td) / "new" / "e10.json"; sys.argv = ["e10_synth.py", "--out", str(out), "--cells", "0", "--reps", "2", "--ndraw", "3", "--fixed-draws", "60"]; S.main(); res = json.load(open(out))
    c = res["cells"]["0"]; ok("tiny end-to-end run writes the cell with library rules, fixed-threshold results and generator checks", set(c["library_rules"]) >= {"S3|0.1", "SPLIT6|0.1", "S3|0.05"} and "0.1" in c["fixed_threshold"] and res["generator_checks"]["icc_tau0.0"][0] < 0.03 and (out.parent / "SHA256SUMS").exists())
print("ALL PASS" if not fails else "FAILED: " + "; ".join(fails))

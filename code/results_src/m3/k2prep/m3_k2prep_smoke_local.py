# Local smoke for m3_replay_lib.py, m3_r2_aux.py, m3_r3_bound.py and m3_w3_baselines.py.  SYNTHETIC score files, manifests and stored "REV-1" aggregates in an
# empty scratch dir.  The only real file read is results_src/mh/mh_rev1.py (code; for the verbatim b1u check).  WRITTEN, NOT RUN, until 实验管理 approves.
# To keep it fast the loaded library is set to 2 seeds x 6 draws (module attributes changed in this process only); synthetic stored aggregates are built from the
# same replay, so the 1e-12 self-check is exercised for real.
# (1) library: b1u copy check passes against mh_rev1.py and fails against an altered copy; fixed_sequence steps down while ub <= eps and returns not certified when
#     the first grid point fails or ub raises; replay is deterministic (two calls equal); self_check is empty on matching aggregates and flags a 1e-9 change.
# (2) R2: main on synthetic sites with field shifts writes cells for 3 sites x 2 candidates x 2 eps; p_f and a_f of one field equal an independent count;
#     the K2 bound ratio (sd_a + sd_d)^2 / Var(p) >= 1 - 1e-12 in every cell (triangle inequality); with field-shifted synthetic data rho(p, a) > 0.5 for c0.
# (3) R3: main writes; contracted coverage equals base coverage in every cell; per cell mean |sens'(1) - sens| <= 1/n_pos (tolerance = 1/n_pos of the site);
#     one cell's improvement (K2+c1s, lacuna_B, eps 0.10, kappa 0.5) equals an independent recomputation; the reading per arm equals the mean of the three
#     eps-0.10 kappa-0.5 improvements and the verdict follows the stop rule (< 0.03 stop, else continue, criterion +0.02 unchanged).
# (4) W3: (W * x).sum(1) equals W @ x element-wise on synthetic weights; gee_ub equals a hand computation; cb_ub with identity-like weights (each row one field repeated) brackets p_hat; bb_ub > mu_hat and shrinks when n grows;
#     (bb: same field proportions with 4x the fields gives a smaller bound)
#     main writes 3 sites x 2 eps x 9 rules and the reference rows S3/B1u equal the synthetic stored aggregates.
# (4b) W3 exception counting: BB forced to raise -> every BB draw counted once (ValueError), flag set, not certified; CB, GEE: 0 exceptions.
# (5) refusals, nothing written: a path containing "iml"; a pinned input differs; a score file that differs from its manifest; output exists; a stored aggregate
#     changed by 1e-6 -> self-check failure (R3 and W3 exit, nothing written).
# Usage: python -I -B m3_k2prep_smoke_local.py <empty scratch dir>
import hashlib, importlib.util, json, shutil, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; T = Path(sys.argv[1]); assert T.is_dir() and not any(T.iterdir())
MH_REV1 = HERE.parents[1] / "mh" / "mh_rev1.py"
sys.path.insert(0, str(HERE))
import m3_replay_lib as L   # noqa: E402
L.SEEDS = [41, 42]; L.NDRAW = 6


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py"); M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M); return M


R2, R3, W3 = load("m3_r2_aux"), load("m3_r3_bound"), load("m3_w3_baselines")
ok = True


def report(tag, cond, extra=""):
    global ok
    ok &= bool(cond); print(tag, "ok" if cond else "FAIL", extra, flush=True)


def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


# ---------- synthetic development sites: 40 fields x 60 cells, field shift mu_f, 5% positives
rng = np.random.default_rng(11); lk = T / "a2_link"; lk.mkdir(); mhlh = T / "frozen"; mhlh.mkdir(); man = ["source\tsource_sha256\tmode\tfile_sha256"]; sums = []
for site in L.DEV:
    nf, nc = 40, 60; g = np.repeat(np.arange(nf), nc); y = (rng.random(nf*nc) < 0.05).astype(int); mu = rng.normal(0, 1.0, nf)[g]
    for sd in L.SEEDS:
        s = mu + rng.normal(0, 1, nf*nc) + 2.5*y; z = np.c_[np.zeros(len(s)), s]
        if site == "bbbc041_all":
            f = lk / f"a2_s{sd}_bbbc041_all_logits.npz"; np.savez(f, logits=z); man.append(f"/x/iter1/a2_s{sd}_scores.npz[bbbc041_all_logits]\tabc\twritten\t{h(f)}")
        else:
            f = mhlh / f"a2_s{sd}_{L.SHORT[site]}_logits.npz"; np.savez(f, logits=z); sums.append(f"{h(f)}  ./{f.name}")
    if site == "bbbc041_all": np.savez(T / "ids_iter1.npz", bbbc041_all_y=y, bbbc041_all_group=g, nih_locked_test_y=np.ones(3, int))
    else: fi = mhlh / f"ids_{L.SHORT[site]}.npz"; np.savez(fi, y=y, field=g); sums.append(f"{h(fi)}  ./{fi.name}")
(lk / "MANIFEST.tsv").write_text("\n".join(man) + "\n"); (mhlh / "SHA256SUMS").write_text("\n".join(sums) + "\n")
stored = {site: T / f"rev1_{site}.json" for site in L.DEV}
for site in L.DEV: json.dump({"agg": {}}, open(stored[site], "w"))
L.PIN.update(manifest=h(lk / "MANIFEST.tsv"), sums=h(mhlh / "SHA256SUMS"), ids_iter1=h(T / "ids_iter1.npz"), mh_rev1_py=h(MH_REV1),
             rev1_bbbc=h(stored["bbbc041_all"]), rev1_mh=h(stored["mpidb_full"]), rev1_lh=h(stored["lacuna_B"]))
inp = L.Inputs(lk, T / "ids_iter1.npz", mhlh, stored["bbbc041_all"], stored["mpidb_full"], stored["lacuna_B"], MH_REV1)
rows = [r for site in L.DEV for sd in L.SEEDS for r in L.replay(inp, sd, site)[0]]
for site in L.DEV:
    agg = {}
    for eps, _ in L.SETTINGS:
        for rule in ("S3", "B1u_c0", "B1u_c1s"):
            R = [r for r in rows if r["site"] == site and r["eps"] == eps and r["rule"] == rule]; fpr = np.array([r["fpr"] for r in R])
            agg[L.STORED_KEY[site].format(rule=rule, eps=eps)] = dict(cov=float((fpr <= eps).mean()), fpr_mean=float(fpr.mean()), fpr_median=float(np.median(fpr)),
                                                                     cannot_certify=float(1-np.mean([r["cert"] for r in R])), ds_eps=float(np.mean([r["dse"] for r in R])),
                                                                     ds_own=float(np.mean([r["dso"] for r in R])), n=len(R))
    json.dump({"agg": agg}, open(stored[site], "w"))
L.PIN.update(rev1_bbbc=h(stored["bbbc041_all"]), rev1_mh=h(stored["mpidb_full"]), rev1_lh=h(stored["lacuna_B"]))
inp = L.Inputs(lk, T / "ids_iter1.npz", mhlh, stored["bbbc041_all"], stored["mpidb_full"], stored["lacuna_B"], MH_REV1)
ARGS = ["--a2link", str(lk), "--ids_iter1", str(T / "ids_iter1.npz"), "--mhlh", str(mhlh), "--rev1_bbbc", str(stored["bbbc041_all"]),
        "--rev1_mh", str(stored["mpidb_full"]), "--rev1_lh", str(stored["lacuna_B"]), "--mh_rev1_py", str(MH_REV1)]

# ---------- (1) library
bad_copy = T / "mh_rev1_altered.py"; bad_copy.write_text(MH_REV1.read_text().replace("ph = X/N", "ph = X/(N+0)", 1))
try: L.check_b1u_copy(bad_copy); e_copy = None
except AssertionError as e: e_copy = repr(e)
seq = L.fixed_sequence([5, 4, 3, 2, 1], lambda j: [0.01, 0.02, 0.2, 0.3, 0.4][j], 0.1)
seq0 = L.fixed_sequence([5, 4], lambda j: 0.5, 0.1); seqx = L.fixed_sequence([5, 4], lambda j: 1/0, 0.1)
rowsB = [r for site in L.DEV for sd in L.SEEDS for r in L.replay(inp, sd, site)[0]]
same = json.dumps(rowsB, sort_keys=True) == json.dumps(rows, sort_keys=True)
agg = json.load(open(stored["lacuna_B"])); k0 = next(iter(agg["agg"])); agg["agg"][k0]["cov"] += 1e-9; pert = T / "rev1_pert.json"; json.dump(agg, open(pert, "w"))
keep = inp.rev1["lacuna_B"]; inp.rev1["lacuna_B"] = pert; flagged = L.self_check(inp, rows); inp.rev1["lacuna_B"] = keep
report("(1) library", L.check_b1u_copy(MH_REV1) is None and e_copy and "differs" in e_copy and seq == (4, True, 1) and seq0 == (5, False, -1) and seqx == (5, False, -1)
       and same and L.self_check(inp, rows) == [] and len(flagged) == 1, (seq, seq0, seqx, len(flagged)))


def run(M, out, extra_args=None):
    sys.argv = ["x"] + (extra_args or ARGS) + ["--out", str(out)]
    try: M.main(); return None
    except (AssertionError, SystemExit) as e: return repr(e)[:180]


# ---------- (2) R2
e2 = run(R2, T / "r2" / "r2.json"); Q = json.loads((T / "r2" / "r2.json").read_text()) if e2 is None else {}
pop = L.population(inp, 41, "mpidb_full"); f0 = pop["fl"][0]; t0 = float(np.quantile(pop["sneg"]["c0"], 0.9))
negs = pop["sc"]["c0"][(pop["g"] == f0) & (pop["y"] == 0)]; alls = pop["sc"]["c0"][pop["g"] == f0]
fs = R2.field_stats(pop, "c0", t0)
P0 = (negs >= t0).mean(); A0 = (alls >= t0).mean()
Pv = [(pop["sc"]["c0"][(pop["g"] == f) & (pop["y"] == 0)] >= t0).mean() for f in pop["fl"]]
cells = Q.get("cells", [])
report("(2) R2", e2 is None and len(cells) == 12 and abs(fs["sd_p"] - float(np.std(Pv))) < 1e-12 and abs(P0 - ((len(pop["fl_neg"]["c0"][f0]) - np.searchsorted(pop["fl_neg"]["c0"][f0], t0, side="left"))/len(pop["fl_neg"]["c0"][f0]))) < 1e-12
       and all(c["mean_over_seeds"]["k2_bound_over_var_p"] is None or c["mean_over_seeds"]["k2_bound_over_var_p"] >= 1 - 1e-12 for c in cells)
       and all(c["mean_over_seeds"]["rho_pa"] > 0.5 for c in cells if c["cand"] == "c0"), (e2, round(A0, 3)))

# ---------- (3) R3
e3 = run(R3, T / "r3" / "r3.json"); Q3 = json.loads((T / "r3" / "r3.json").read_text()) if e3 is None else {}
c3 = Q3.get("cells", [])
cov_same = all(abs(c["by_kappa"][k]["coverage_contracted"] - c["by_kappa"][k]["coverage_base"]) < 1e-12 for c in c3 for k in c["by_kappa"])
rd_ok = True
for arm in R3.ARMS:
    m = float(np.mean([c["by_kappa"]["0.5"]["improvement"] for c in c3 if c["arm"] == arm and abs(c["eps"] - 0.1) < 1e-12])) if c3 else None
    rd = Q3.get("reading", {}).get(arm, {})
    rd_ok &= m is not None and abs(rd.get("mean_improvement_eps010_kappa05", 9) - m) < 1e-12 and rd.get("verdict", "").startswith("stop" if m < 0.03 else "continue")
tol_ok = all(c["discretisation_mean_abs_kappa1"] <= c["discretisation_tol"] and abs(c["discretisation_tol"] - 1.0/c["n_pos"]) < 1e-15 for c in c3)
# independent recomputation of one cell's improvement (arm K2+c1s, lacuna_B, eps 0.10, kappa 0.5) without R3.contracted
ind_ok = False
if c3:
    eps0 = 0.10; Rc = [r for r in rows if r["site"] == "lacuna_B" and abs(r["eps"] - eps0) < 1e-12 and r["rule"] == "B1u_c1s"]; pp = {sd: L.population(inp, sd, "lacuna_B") for sd in L.SEEDS}
    def sprime(r, kap):
        if not r["cert"]: return r["sens"]
        P = pp[r["seed"]]; neg = P["sneg"]["c1s"]; pos = P["spos"]["c1s"]; t = np.quantile(neg, 1 - (eps0 - kap*(eps0 - r["fpr"])))
        return float((pos >= t).mean())
    impI = np.mean([sprime(r, 0.5) for r in Rc]) - np.mean([sprime(r, 1.0) for r in Rc])
    cc = [c for c in c3 if c["arm"] == "K2+c1s" and c["site"] == "lacuna_B" and abs(c["eps"] - eps0) < 1e-12][0]
    ind_ok = abs(cc["by_kappa"]["0.5"]["improvement"] - impI) < 1e-12
report("(3) R3", e3 is None and len(c3) == 12 and cov_same and tol_ok and ind_ok and rd_ok, e3 or {a: round(Q3["reading"][a]["mean_improvement_eps010_kappa05"], 4) for a in R3.ARMS})

# ---------- (4) W3
x = np.array([0, 2, 1, 4]); n = np.array([20, 25, 30, 22]); p = x.sum()/n.sum(); V = 4/3*((x - p*n)**2).sum()/n.sum()**2
from scipy.stats import t as tdist   # noqa: E402
hand = p + tdist.ppf(0.9, 3)*np.sqrt(V)
Wid = np.repeat(np.eye(4)*4, 1, axis=0)
cb = W3.cb_ub(x, n, Wid)
Wr = np.random.default_rng(5).multinomial(4, np.full(4, 0.25), size=500).astype(float); same_sum = bool(np.array_equal((Wr * x).sum(1), Wr @ x) and np.array_equal((Wr * n).sum(1), Wr @ n))
x1, n1 = np.array([1, 2, 1, 3, 2, 1]), np.array([40]*6); bb1 = W3.bb_ub(x1, n1); bb2 = W3.bb_ub(np.tile(x1, 4), np.tile(n1, 4))   # same proportions, 4x the fields
e4 = run(W3, T / "w3" / "w3.json"); Q4 = json.loads((T / "w3" / "w3.json").read_text()) if e4 is None else {}
refs_ok = True
for c in Q4.get("cells", []):
    if c["rule"] in ("S3", "B1u_c0", "B1u_c1s"):
        st = json.load(open(stored[c["site"]]))["agg"][L.STORED_KEY[c["site"]].format(rule=c["rule"], eps=c["eps"])]
        refs_ok &= abs(c["coverage"] - st["cov"]) < 1e-12 and abs(c["ds_own"] - st["ds_own"]) < 1e-12
report("(4) W3", same_sum and abs(W3.gee_ub(x, n) - hand) < 1e-12 and min(x/n) - 1e-12 <= cb <= max(x/n) + 1e-12 and bb1 > 10/240 and bb2 < bb1 and e4 is None
       and len(Q4["cells"]) == 3*2*9 and refs_ok, (round(bb1, 4), round(bb2, 4), e4))

# ---------- (4b) W3 exception counting: BB made to raise on every call -> each BB draw stops at its first grid point (not certified), every draw counted once,
#      type ValueError, flag set; CB and GEE have no exceptions; S3/B1u unaffected
_bb = W3.bb_ub


def _boom(x, n, delta=0.10): raise ValueError("forced")


W3.bb_ub = _boom; e4b = run(W3, T / "w3b" / "w3.json"); W3.bb_ub = _bb
Q4b = json.loads((T / "w3b" / "w3.json").read_text()) if e4b is None else {}
cb4 = Q4b.get("cells", []); bbc = [c for c in cb4 if c["rule"].startswith("BB_")]; oth = [c for c in cb4 if c["rule"].startswith(("CB_", "GEE_"))]
ref_same = all(abs(c["coverage"] - d["coverage"]) < 1e-12 for c in cb4 for d in Q4.get("cells", []) if c["rule"] in ("S3", "B1u_c0", "B1u_c1s") and (c["site"], c["eps"], c["rule"]) == (d["site"], d["eps"], d["rule"]))
report("(4b) W3 exception counting", e4b is None and len(bbc) == 12 and all(c["n_bound_exceptions"] == c["n"] and c["exception_types"] == [["ValueError", c["n"]]]
       and c["flag"] == "exceptions>1%" and c["cannot_certify"] == 1.0 for c in bbc) and all(c["n_bound_exceptions"] == 0 and c["flag"] is None for c in oth) and ref_same,
       (e4b, bbc[0]["n_bound_exceptions"] if bbc else None))

# ---------- (5) refusals
e5a = run(R3, T / "r5a" / "x.json", [a.replace(str(lk), str(T / "iml_link")) for a in ARGS])
keep = L.PIN["sums"]; L.PIN["sums"] = "0"*64; e5b = run(R2, T / "r5b" / "x.json"); L.PIN["sums"] = keep
bl = T / "a2_link_bad"; shutil.copytree(lk, bl); np.savez(bl / "a2_s42_bbbc041_all_logits.npz", logits=np.zeros((5, 2)))
e5c = run(R2, T / "r5c" / "x.json", [a.replace(str(lk), str(bl)) for a in ARGS])
e5d = run(R3, T / "r3" / "r3.json")
agg = json.load(open(stored["mpidb_full"])); k0 = next(iter(agg["agg"])); agg["agg"][k0]["ds_own"] += 1e-6; bad = T / "rev1_mh_bad.json"; json.dump(agg, open(bad, "w"))
keep = L.PIN["rev1_mh"]; L.PIN["rev1_mh"] = h(bad); args_bad = [a.replace(str(stored["mpidb_full"]), str(bad)) for a in ARGS]
e5e = run(R3, T / "r5e" / "x.json", args_bad); e5f = run(W3, T / "r5f" / "x.json", args_bad); L.PIN["rev1_mh"] = keep
report("(5) refusals, nothing written", e5a and "confirmatory" in e5a and e5b and "pinned" in e5b and e5c and "sha" in e5c and e5d and "exists" in e5d
       and e5e and "self-check" in e5e and e5f and "self-check" in e5f and not any((T / d).exists() for d in ("r5a", "r5b", "r5c", "r5e", "r5f")),
       " | ".join(str(v)[:45] for v in (e5a, e5b, e5c, e5d, e5e, e5f)))
print("ALL PASS" if ok else "FAIL")

# Local smoke for m3_k12_dev_eval.py (draft5 b3318a71 + addendum A).  SYNTHETIC score files, manifests and stored "REV-1" aggregates in an empty scratch dir; the
# only real file read is results_src/mh/mh_rev1.py (code; verbatim b1u check).  WRITTEN, NOT RUN, until 实验管理 approves.  Library set to 2 seeds x 6 draws here.
# (1) interface (addendum A2): to_seq(UB) <= eps iff UB <= 0, for UB in {0, -1e-12, +1e-12, +inf, NaN} and eps in {0.05, 0.10}.
# (2) AllCells: sum_F w_f(t) and the population sd of w_f(t) equal a direct per-field computation at several thresholds; g_lab and c_lab equal direct counts.
# (3) identity (addendum A3): on a synthetic field, e_f = z_f - w_f equals -P_f (TPR_f(t) - eps).
# (4) K1: the t form equals a hand formula; on 200 normal residuals the bootstrap-t bound is within 0.3 SE of the t bound; s_e = 0 (constant residual,
#     floating-point sd ~ 1e-16) -> +inf and counted as se0, for K1 and for the control K1t;
#     a weight matrix whose resamples are mostly degenerate (> 5% dropped) -> +inf and counted.
# (5) K2: the chi-square sd bound equals the hand formula; on 200 normal residuals the bootstrap-t sd bound is within 10% of the chi-square one; s_e = 0 -> +inf
#     and counted as se0, for K2 and for the control K2chi.
# (6) main on synthetic sites: writes 6 arms x 3 sites x 2 eps cells; S3 and B1u_c1s rows equal the synthetic stored aggregates; judgement has K1, K2 and choice;
#     every K arm cell reports n_se0, n_drop_over_5pct and exception counts.
# (7) exceptions and refusals: k1_ub forced to raise -> K1 and K1t cells count one exception per draw and candidate (3 per draw) and are flagged; output exists -> refused; a stored aggregate
#     changed by 1e-6 -> self-check failure, nothing written.
# Usage: python -I -B m3_k12_smoke_local.py <empty scratch dir>
import hashlib, importlib.util, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import chi2, t as tdist

HERE = Path(__file__).resolve().parent; T = Path(sys.argv[1]); assert T.is_dir() and not any(T.iterdir())
MH_REV1 = HERE.parents[1] / "mh" / "mh_rev1.py"
sys.path.insert(0, str(HERE))
import m3_replay_lib as L   # noqa: E402
L.SEEDS = [41, 42]; L.NDRAW = 6
spec = importlib.util.spec_from_file_location("ev", HERE / "m3_k12_dev_eval.py"); E = importlib.util.module_from_spec(spec); spec.loader.exec_module(E)
ok = True


def report(tag, cond, extra=""):
    global ok
    ok &= bool(cond); print(tag, "ok" if cond else "FAIL", extra, flush=True)


def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


# (1) interface
vals = [0.0, -1e-12, 1e-12, np.inf, np.nan]
report("(1) to_seq interface", all((E.to_seq(u) <= e) == (bool(np.isfinite(u)) and u <= 0) for u in vals for e in (0.05, 0.10)), [E.to_seq(u) for u in vals])

# synthetic development sites (same construction as m3_k2prep_smoke_local.py)
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

# (2) AllCells vs direct computation
pop = L.population(inp, 41, "lacuna_B"); okc = True
for k in ("c0", "c1s"):
    A = E.AllCells(pop, k); sc, gg = pop["sc"][k], pop["g"]; fl = list(pop["fl"])
    for q in (0.5, 0.9, 0.97, 0.995):
        t = float(np.quantile(sc, q)); eps = 0.05
        gf = np.array([(sc[gg == f] >= t).sum() for f in fl], float); cf = np.array([(gg == f).sum() for f in fl], float); w = gf - eps*cf
        Sw, sw = A.w_stats(t, eps); okc &= abs(Sw - w.sum()) < 1e-9 and abs(sw - w.std()) < 1e-9
        Fl = fl[:5]; okc &= np.array_equal(A.g_lab(Fl, t), gf[:5]) and np.array_equal(A.c_lab(Fl), cf[:5])
report("(2) AllCells exact", okc)

# (3) identity e_f = -P_f (TPR_f - eps)
k = "c0"; f0 = list(pop["fl"])[3]; sc, gg, yy = pop["sc"][k], pop["g"], pop["y"]; t = float(np.quantile(sc, 0.9)); eps = 0.10
neg = sc[(gg == f0) & (yy == 0)]; pos = sc[(gg == f0) & (yy == 1)]; allc = sc[gg == f0]
zf = (neg >= t).sum() - eps*len(neg); wf = (allc >= t).sum() - eps*len(allc); Pf = len(pos); tprf = (pos >= t).mean() if Pf else 0.0
report("(3) identity e = -P (TPR - eps)", abs((zf - wf) - (-Pf*(tprf - eps))) < 1e-9, (Pf, round(zf - wf, 4)))

# (4) K1
r4 = np.random.default_rng(4); m = 22; e = r4.normal(0, 1, m); F = 400; Sw = -50.0
SE = F*e.std(ddof=1)/np.sqrt(m)*np.sqrt(1 - m/F); hand = Sw + F*e.mean() + tdist.ppf(1 - E.ALPHA, m - 1)*SE
W = r4.multinomial(m, np.full(m, 1/m), size=E.B).astype(float); E.DEG.clear()
okt = abs(E.k1_ub(e, Sw, F, m, W, "t", "x") - hand) < 1e-9
eB = r4.normal(0, 1, 200); WB = r4.multinomial(200, np.full(200, 1/200), size=E.B).astype(float); SEB = F*eB.std(ddof=1)/np.sqrt(200)*np.sqrt(1 - 200/F)
close = abs(E.k1_ub(eB, Sw, F, 200, WB, "bt", "x") - E.k1_ub(eB, Sw, F, 200, WB, "t", "x")) < 0.3*SEB
inf0 = E.k1_ub(np.full(m, 0.7), Sw, F, m, W, "bt", "deg") == np.inf
Wdeg = np.zeros_like(W); Wdeg[:, 0] = m; Wdeg[:10] = W[:10]   # >5% of resamples put all weight on one field -> sd* = 0
infd = E.k1_ub(e, Sw, F, m, Wdeg, "bt", "deg2") == np.inf
inf_t = E.k1_ub(np.full(m, 0.7), Sw, F, m, W, "t", "deg3") == np.inf                      # control arm K1t on a constant residual
report("(4) K1 bounds", okt and close and inf0 and infd and inf_t and E.DEG["deg"]["se0"] == 1 and E.DEG["deg2"]["drop>5%"] == 1 and E.DEG["deg3"]["se0"] == 1,
       (okt, close, inf0, infd, inf_t, dict(E.DEG["deg"])))

# (5) K2 sd bound
sd_hand = e.std(ddof=1)*np.sqrt((m - 1)/chi2.ppf(E.A_D, m - 1))
okchi = abs(E.sd_upper(e, m, W, "chi", "y") - sd_hand) < 1e-12
bt200, chi200 = E.sd_upper(eB, 200, WB, "bt", "y"), E.sd_upper(eB, 200, WB, "chi", "y")
okbt = abs(bt200 - chi200) < 0.10*chi200
ok0 = E.k2_ub(np.zeros(m), np.full(m, 0.3), 0.1, F, m, W, "bt", "y2") == np.inf
ok0c = E.k2_ub(np.zeros(m), np.full(m, 0.3), 0.1, F, m, W, "chi", "y3") == np.inf      # control arm K2chi on a constant residual
report("(5) K2 sd bound", okchi and okbt and ok0 and ok0c and E.DEG["y2"]["se0"] == 1 and E.DEG["y3"]["se0"] == 1, (round(bt200, 4), round(chi200, 4)))


def run(out, args=None):
    sys.argv = ["ev"] + (args or ARGS) + ["--out", str(out)]
    try: E.main(); return None
    except (AssertionError, SystemExit) as ex: return repr(ex)[:180]


# (6) main
e6 = run(T / "o6" / "ev.json"); Q = json.loads((T / "o6" / "ev.json").read_text()) if e6 is None else {}
cells = Q.get("cells", []); ref = True
for c in cells:
    if c["arm"] in ("S3", "B1u_c1s"):
        st = json.load(open(stored[c["site"]]))["agg"][L.STORED_KEY[c["site"]].format(rule=c["arm"], eps=c["eps"])]
        ref &= abs(c["coverage"] - st["cov"]) < 1e-12 and abs(c["ds_own"] - st["ds_own"]) < 1e-12 and abs(c["cannot_certify"] - st["cannot_certify"]) < 1e-12
kfields = all(all(k_ in c for k_ in ("n_se0", "n_drop_over_5pct", "n_bound_exceptions")) for c in cells if c["arm"] in E.K_ARMS)
jd = Q.get("judgement", {})
report("(6) main on synthetic", e6 is None and len(cells) == 36 and ref and kfields and "K1" in jd and "K2" in jd and "choice" in jd,
       e6 or {a: (jd[a]["passes"], round(jd[a]["mean_gain_eps010"], 4)) for a in ("K1", "K2")})

# (7) exceptions and refusals
_k1 = E.k1_ub


def _boom(*a, **k): raise ValueError("forced")


E.k1_ub = _boom; e7 = run(T / "o7" / "ev.json"); E.k1_ub = _k1
Q7 = json.loads((T / "o7" / "ev.json").read_text()) if e7 is None else {}
k1c = [c for c in Q7.get("cells", []) if c["arm"] in ("K1", "K1t")]; k2c = [c for c in Q7.get("cells", []) if c["arm"] in ("K2", "K2chi")]
exc_ok = e7 is None and len(k1c) == 12 and all(c["n_bound_exceptions"] == 3*c["n"] and c["flag"] == "exceptions>1%" and c["cannot_certify"] == 1.0 for c in k1c) \
    and all(c["n_bound_exceptions"] == 0 for c in k2c) and Q7["judgement"]["K1"]["passes"] is False
e7b = run(T / "o6" / "ev.json")
agg = json.load(open(stored["mpidb_full"])); k0 = next(iter(agg["agg"])); agg["agg"][k0]["ds_own"] += 1e-6; bad = T / "rev1_mh_bad.json"; json.dump(agg, open(bad, "w"))
keep = L.PIN["rev1_mh"]; L.PIN["rev1_mh"] = h(bad); e7c = run(T / "o7c" / "ev.json", [a.replace(str(stored["mpidb_full"]), str(bad)) for a in ARGS]); L.PIN["rev1_mh"] = keep
report("(7) exceptions and refusals", exc_ok and e7b and "exists" in e7b and e7c and "self-check" in e7c and not (T / "o7c").exists(),
       " | ".join(str(v)[:50] for v in (e7, e7b, e7c)))
print("ALL PASS" if ok else "FAIL")

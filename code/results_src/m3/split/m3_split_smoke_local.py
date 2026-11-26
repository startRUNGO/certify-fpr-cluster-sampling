# Local smoke for m3_split_dev_eval.py.  SYNTHETIC sites (3 rebuilt-style + old), synthetic NIH-validation files, synthetic stored "REV-1" aggregates from the same
# replay; only real file read: results_src/mh/mh_rev1.py.  5 seeds x 4 draws (this process only).  WRITTEN, NOT RUN, until 实验管理 approves.
# (1) split_to_seq: 0.0 only for (True, True); (True, False), (False, True), (False, False) -> 1.0.
# (2) anova_icc equals an independent one-way ANOVA estimate; select_bound widths equal hand formulas; ties resolve b1u > K1 > K1-v.
# (3) QCells.sum_v and v_fields equal brute-force sums over cells at several thresholds (ties included).
# (4) every certified SPLIT row has p_s <= eps and C1 true at the certified threshold (Proposition 2' conditions recorded per row).
# (5) the copied loop: K1 and K1vB rows equal m3_k1vb_dev_eval.run_site_seed rows; S3 rows equal the library replay (the split stream does not disturb REV-1).
# (6) main writes 7 arms x 4 sites x 2 eps = 56 cells; SPLIT cells report bound_chosen (counts sum to the draws) and c2_bottleneck_points; judgement for SPLIT6 and
#     report-only judgements for SPLIT4 / SPLIT8.
# (7) refusals, nothing written: output exists; a stored aggregate changed by 1e-6 -> self-check failure.
# Usage: python -I -B m3_split_smoke_local.py <empty scratch dir>
import argparse, hashlib, importlib.util, json, sys
from pathlib import Path
import numpy as np
from scipy.special import expit
from scipy.stats import norm

HERE = Path(__file__).resolve().parent; T = Path(sys.argv[1]); assert T.is_dir() and not any(T.iterdir())
MH_REV1 = HERE.parents[1] / "mh" / "mh_rev1.py"
for sub in ("k2prep", "stress", "k1vb"): sys.path.insert(0, str(HERE.parents[0] / sub))
import m3_replay_lib as L   # noqa: E402
L.SEEDS = [41, 42, 43, 44, 45]; L.NDRAW = 4
spec = importlib.util.spec_from_file_location("sp", HERE / "m3_split_dev_eval.py"); P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
ST, KV, VB, E = P.ST, P.KV, P.VB, P.E
ok = True


def report(tag, cond, extra=""):
    global ok
    ok &= bool(cond); print(tag, "ok" if cond else "FAIL", extra, flush=True)


def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


# (1)
report("(1) split_to_seq truth table", [P.split_to_seq(a, b) for a, b in ((True, True), (True, False), (False, True), (False, False))] == [0.0, 1.0, 1.0, 1.0])

# (2)
r = np.random.default_rng(5); x = r.integers(0, 8, 6).astype(float); n = r.integers(40, 80, 6).astype(float)
N = n.sum(); ph = x.sum()/N; pi = x/n; n0 = (N - (n**2).sum()/N)/5; msb = (n*(pi - ph)**2).sum()/5; msw = (n*pi*(1 - pi)).sum()/(N - 6); th = (msb/msw - 1)/n0
icc_ind = min(max(th/(1 + th), 0.0), 0.999) if th > 0 else 0.0
e = r.normal(0, 2, 6); v = -e + r.normal(0, 0.3, 6); b, W = P.select_bound(e, v, x, n, 16, 400)
zs = norm.ppf(1 - 0.10/3); den = np.sqrt(16)*n.mean(); fpc = np.sqrt(1 - 16/400)
hand = {"K1": zs*np.std(e, ddof=1)/den*fpc, "K1-v": zs*np.std(e + v, ddof=1)/den*fpc, "b1u": zs*np.sqrt(ph*(1 - ph)*(1 + (n.mean() - 1)*icc_ind)/(16*n.mean()))}
ties = P.select_bound(np.zeros(6), np.zeros(6), np.zeros(6), np.full(6, 50.0), 16, 400)[0]
report("(2) anova_icc, widths, ties", abs(P.anova_icc(x, n) - icc_ind) < 1e-12 and all(abs(W[k] - hand[k]) < 1e-12 for k in hand) and b == min(hand, key=hand.get) and ties == "b1u",
       (b, ties))

# synthetic sites
rng = np.random.default_rng(29); lk = T / "a2_link"; lk.mkdir(); mhlh = T / "frozen"; mhlh.mkdir(); old = T / "old"; old.mkdir(); sums = []; ids = {}; man = {}
for site in L.DEV + ("old",):
    nf, nc = 40, 45; g = np.repeat(np.arange(nf), nc); y = (rng.random(nf*nc) < (0.3 if site == "old" else 0.05)).astype(int); mu = rng.normal(0, 1.0, nf)[g]
    for sd in L.SEEDS:
        sc = mu + rng.normal(0, 1, nf*nc) + 2.5*y; z = np.c_[np.zeros(len(sc)), sc]
        if site == "bbbc041_all":
            f = lk / f"a2_s{sd}_bbbc041_all_logits.npz"; np.savez(f, logits=z); man.setdefault(sd, {})["file"] = h(f)
        elif site == "old":
            if sd == 41: f = T / "iter1_a2_s41_scores.npz"; np.savez(f, mpidb_full_logits=z, nih_locked_test_logits=np.ones((3, 2))); man.setdefault(41, {})["src"] = h(f)
            else: f = old / f"a2_s{sd}_mpidb_full_logits.npz"; np.savez(f, logits=z)
        else:
            f = mhlh / f"a2_s{sd}_{L.SHORT[site]}_logits.npz"; np.savez(f, logits=z); sums.append(f"{h(f)}  ./{f.name}")
    if site in L.SHORT: fi = mhlh / f"ids_{L.SHORT[site]}.npz"; np.savez(fi, y=y, field=g); sums.append(f"{h(fi)}  ./{fi.name}")
    ids[site] = (y, g)
(old / "SHA256SUMS").write_text("".join(f"{h(old / f'a2_s{sd}_mpidb_full_logits.npz')}  ./a2_s{sd}_mpidb_full_logits.npz\n" for sd in (42, 43, 44, 45)))
(lk / "MANIFEST.tsv").write_text("source\tsource_sha256\tmode\tfile_sha256\n" + "".join(
    f"/x/iter1/a2_s{sd}_scores.npz[bbbc041_all_logits]\t{man[sd].get('src', 'abc')}\twritten\t{man[sd]['file']}\n" for sd in L.SEEDS))
(mhlh / "SHA256SUMS").write_text("\n".join(sums) + "\n")
np.savez(T / "ids_iter1.npz", bbbc041_all_y=ids["bbbc041_all"][0], bbbc041_all_group=ids["bbbc041_all"][1], mpidb_full_y=ids["old"][0], mpidb_full_group=ids["old"][1],
         nih_locked_test_y=np.ones(3, int))
nv = T / "a2ref"; nv.mkdir(); yv = (rng.random(800) < 0.5).astype(int)
np.savez(nv / "ids_nih_validation.npz", name=np.arange(800).astype(str), y=yv, field=np.zeros(800, int))
for sd in L.SEEDS: np.savez(nv / f"a2ref_s{sd}_nih_validation_logits.npz", logits=np.c_[np.zeros(800), rng.normal(0, 1, 800) + 3*yv])
KV.NIHVAL_PIN.update({sd: h(nv / f"a2ref_s{sd}_nih_validation_logits.npz") for sd in L.SEEDS}); KV.NIHVAL_PIN["ids"] = h(nv / "ids_nih_validation.npz")
stored = {k: T / f"rev1_{k}.json" for k in ("bbbc", "mh", "lh")}
for p in stored.values(): json.dump({"agg": {}}, open(p, "w"))


def pins():
    L.PIN.update(manifest=h(lk / "MANIFEST.tsv"), sums=h(mhlh / "SHA256SUMS"), ids_iter1=h(T / "ids_iter1.npz"), mh_rev1_py=h(MH_REV1),
                 rev1_bbbc=h(stored["bbbc"]), rev1_mh=h(stored["mh"]), rev1_lh=h(stored["lh"]))


def ns(out="x", rev1_bbbc=None):
    return argparse.Namespace(a2link=str(lk), ids_iter1=str(T / "ids_iter1.npz"), mhlh=str(mhlh), rev1_bbbc=rev1_bbbc or str(stored["bbbc"]), rev1_mh=str(stored["mh"]),
                              rev1_lh=str(stored["lh"]), mh_rev1_py=str(MH_REV1), old_s41=str(T / "iter1_a2_s41_scores.npz"), old_dir=str(old),
                              nihval_dir=str(nv), out=out)


pins(); inp = ST.StressInputs(ns()); dev0 = tuple(L.DEV); L.SITECFG[P.OLD] = (8, 1); L.DEV = dev0 + (P.OLD,); L.M = 22
agg = {k: {} for k in stored}; replay_rows = {}
for site, fk in (("bbbc041_all", "bbbc"), ("mpidb_full", "mh"), ("lacuna_B", "lh"), (P.OLD, "bbbc")):
    rows = [x_ for sd in L.SEEDS for x_ in L.replay(inp, sd, site)[0]]; replay_rows[site] = rows
    for eps, _ in L.SETTINGS:
        for rule in ("S3", "B1u_c0", "B1u_c1s"):
            R = [x_ for x_ in rows if abs(x_["eps"] - eps) < 1e-12 and x_["rule"] == rule]; fpr = np.array([x_["fpr"] for x_ in R])
            agg[fk][ST.PATTERN[site].format(m=22, rule=rule, eps=eps)] = dict(cov=float((fpr <= eps).mean()), fpr_mean=float(fpr.mean()), cannot_certify=float(1-np.mean([x_["cert"] for x_ in R])),
                                                                             ds_eps=float(np.mean([x_["dse"] for x_ in R])), ds_own=float(np.mean([x_["dso"] for x_ in R])), n=len(R))
for k, p in stored.items(): json.dump({"agg": agg[k]}, open(p, "w"))
pins(); inp = ST.StressInputs(ns())

# (3) QCells
pop = L.population(inp, 41, "mpidb_full"); q = expit(1.1*pop["s"] - 0.4); okq = True
for k in ("c0", "c1s"):
    QC = P.QCells(pop["sc"][k], q, pop["g"], pop["fl"]); sk = pop["sc"][k]
    for t in (float(np.quantile(sk, 0.5)), float(np.quantile(sk, 0.93)), float(sk[17])):
        brute_all = float((q*((sk >= t).astype(float) - 0.1)).sum()); fs = list(pop["fl"])[:6]
        brute_f = np.array([(q[pop["g"] == f]*((sk[pop["g"] == f] >= t).astype(float) - 0.1)).sum() for f in fs])
        okq &= abs(QC.sum_v(t, 0.1) - brute_all) < 1e-9 and np.allclose(QC.v_fields(fs, t, 0.1), brute_f, atol=1e-9)
report("(3) QCells exact", okq)

# (4) + (5)
nvd = KV.load_nihval(nv); platt = {sd: KV.platt(*nvd[sd]) for sd in L.SEEDS}
sp_rows = [x_ for s in ("lacuna_B", P.OLD) for sd in L.SEEDS for x_ in P.run_site_seed(inp, sd, s, platt[sd])]
cert_split = [x_ for x_ in sp_rows if x_["rule"].startswith("SPLIT") and x_["cert"]]
ok4 = len(cert_split) > 0 and all(x_["split_ps_at"] is not None and x_["split_ps_at"] <= x_["eps"] + 1e-15 and x_["split_c1_at"] is True for x_ in cert_split)
report("(4) Proposition 2' conditions hold at every certified SPLIT threshold", ok4, len(cert_split))
key = lambda R: [(x_["site"], x_["seed"], x_["draw"], x_["eps"], x_["cand"], x_["cert"], round(x_["fpr"], 12)) for x_ in R]
vb_rows = [x_ for s in ("lacuna_B", P.OLD) for sd in L.SEEDS for x_ in VB.run_site_seed(inp, sd, s)]
same_k = all(key([x_ for x_ in sp_rows if x_["rule"] == a]) == key([x_ for x_ in vb_rows if x_["rule"] == a]) for a in ("K1", "K1vB"))
same_s3 = key([x_ for x_ in sp_rows if x_["rule"] == "S3"]) == key([x_ for x_ in replay_rows["lacuna_B"] + replay_rows[P.OLD] if x_["rule"] == "S3"])
report("(5) copied loop: K1/K1vB == VB loop, S3 == library replay", same_k and same_s3, (same_k, same_s3))
L.DEV = dev0; del L.SITECFG[P.OLD]


def run(a):
    sys.argv = ["sp"] + sum([[f"--{k}", v] for k, v in vars(a).items()], [])
    keep = (tuple(L.DEV), dict(L.SITECFG), L.M)
    try: P.main(); return None
    except (AssertionError, SystemExit) as ex: return repr(ex)[:180]
    finally: L.DEV, L.M = keep[0], keep[2]; L.SITECFG.clear(); L.SITECFG.update(keep[1])


# (6)
e6 = run(ns(out=str(T / "o6" / "sp.json"))); Q = json.loads((T / "o6" / "sp.json").read_text()) if e6 is None else {}
sc = [c for c in Q.get("cells", []) if c["arm"].startswith("SPLIT")]
cnt_ok = all(sum(c["bound_chosen"][k].values()) == c["n"] for c in sc for k in L.CEN) and all("c2_bottleneck_points" in c for c in sc)
report("(6) main", e6 is None and len(Q["cells"]) == 56 and cnt_ok and Q["judgement"]["arm"] == "SPLIT6" and [j["arm"] for j in Q["sensitivity_judgement_report_only"]] == ["SPLIT4", "SPLIT8"],
       e6 or (Q["judgement"]["passes"], round(Q["judgement"]["mean_gain_rebuilt_eps010"], 4)))
# (7)
e7a = run(ns(out=str(T / "o6" / "sp.json")))
ag = json.load(open(stored["mh"])); kk = next(iter(ag["agg"])); ag["agg"][kk]["ds_own"] += 1e-6; bad = T / "rev1_mh_bad.json"; json.dump(ag, open(bad, "w"))
keep = L.PIN["rev1_mh"]; L.PIN["rev1_mh"] = h(bad); a7 = ns(out=str(T / "o7" / "sp.json")); a7.rev1_mh = str(bad); e7b = run(a7); L.PIN["rev1_mh"] = keep
report("(7) refusals, nothing written", e7a and "exists" in e7a and e7b and "self-check" in e7b and not (T / "o7").exists(), (str(e7a)[:40], str(e7b)[:50]))
print("ALL PASS" if ok else "FAIL")

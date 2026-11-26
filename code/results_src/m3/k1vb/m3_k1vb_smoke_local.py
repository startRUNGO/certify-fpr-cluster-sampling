# Local smoke for m3_k1vb_dev_eval.py.  SYNTHETIC 3 rebuilt-style sites + synthetic old MP-IDB files and stored "REV-1" aggregates (built from the same replay);
# the only real file read is results_src/mh/mh_rev1.py (b1u check).  5 seeds x 6 draws (module attributes, this process only).  WRITTEN, NOT RUN, until approved.
# (1) k1_ub_alpha(alpha = delta/3) equals m3_k12_dev_eval.k1_ub(..., "bt") exactly on random residuals; constant residual -> +inf for both.
# (2) the OR rule, checked against independent references:
#     (a) K1 part disabled (k1_ub_alpha -> +inf): K1vB rows equal S3 rows of the library replay run with delta/3 replaced by delta/6 (cand, cert, fpr per draw);
#     (b) b1u part disabled at delta/6 (b1u at level delta/6 -> 1.0): K1vB rows equal the K1 arm run with alpha delta/6 (cand, cert, fpr per draw).
# (3) the copied loop: K1 rows of this script equal K1 rows of m3_k12_dev_eval.run_site_seed on the rebuilt sites (cand, cert, fpr per draw).
# (4) main: writes 4 arms x 4 sites x 2 eps = 32 cells; all 4 self-check statuses passed; judgement has gate (8 cells), main and passes.
# (5) refusals, nothing written: old dir without seed 45 -> refused; output exists; stored aggregate changed by 1e-6 -> self-check failure.
# Usage: python -I -B m3_k1vb_smoke_local.py <empty scratch dir>
import argparse, hashlib, importlib.util, json, shutil, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; T = Path(sys.argv[1]); assert T.is_dir() and not any(T.iterdir())
MH_REV1 = HERE.parents[1] / "mh" / "mh_rev1.py"
sys.path.insert(0, str(HERE.parents[0] / "k2prep")); sys.path.insert(0, str(HERE.parents[0] / "stress"))
import m3_replay_lib as L   # noqa: E402
L.SEEDS = [41, 42, 43, 44, 45]; L.NDRAW = 6
spec = importlib.util.spec_from_file_location("vb", HERE / "m3_k1vb_dev_eval.py"); V = importlib.util.module_from_spec(spec); spec.loader.exec_module(V)
E, ST = V.E, V.ST
ok = True


def report(tag, cond, extra=""):
    global ok
    ok &= bool(cond); print(tag, "ok" if cond else "FAIL", extra, flush=True)


def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


# (1)
r = np.random.default_rng(1); m = 22; Wb = r.multinomial(m, np.full(m, 1/m), size=E.B).astype(float); same = True
for _ in range(20):
    e = r.normal(0, 3, m); Sw = float(r.normal(-40, 10)); F = 300
    same &= V.k1_ub_alpha(e, Sw, F, m, Wb, E.ALPHA, "a") == E.k1_ub(e, Sw, F, m, Wb, "bt", "b")
const = V.k1_ub_alpha(np.full(m, 0.3), 1.0, 300, m, Wb, V.A_VB, "c") == np.inf and E.k1_ub(np.full(m, 0.3), 1.0, 300, m, Wb, "bt", "d") == np.inf
report("(1) k1_ub_alpha == k1_ub at delta/3", same and const)

# synthetic sites (as m3_stress_smoke_local.py)
rng = np.random.default_rng(17); lk = T / "a2_link"; lk.mkdir(); mhlh = T / "frozen"; mhlh.mkdir(); old = T / "old"; old.mkdir(); sums = []; ids = {}; man = {}
for site in L.DEV + ("old",):
    nf, nc = 40, 50; g = np.repeat(np.arange(nf), nc); pos_rate = 0.30 if site == "old" else 0.05
    y = (rng.random(nf*nc) < pos_rate).astype(int); mu = rng.normal(0, 1.0, nf)[g]
    for sd in L.SEEDS:
        s = mu + rng.normal(0, 1, nf*nc) + 2.5*y; z = np.c_[np.zeros(len(s)), s]
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
stored = {k: T / f"rev1_{k}.json" for k in ("bbbc", "mh", "lh")}
for p in stored.values(): json.dump({"agg": {}}, open(p, "w"))


def pins():
    L.PIN.update(manifest=h(lk / "MANIFEST.tsv"), sums=h(mhlh / "SHA256SUMS"), ids_iter1=h(T / "ids_iter1.npz"), mh_rev1_py=h(MH_REV1),
                 rev1_bbbc=h(stored["bbbc"]), rev1_mh=h(stored["mh"]), rev1_lh=h(stored["lh"]))


def ns(old_dir=None, out="x", rev1_bbbc=None):
    return argparse.Namespace(a2link=str(lk), ids_iter1=str(T / "ids_iter1.npz"), mhlh=str(mhlh), rev1_bbbc=rev1_bbbc or str(stored["bbbc"]), rev1_mh=str(stored["mh"]),
                              rev1_lh=str(stored["lh"]), mh_rev1_py=str(MH_REV1), old_s41=str(T / "iter1_a2_s41_scores.npz"), old_dir=old_dir or str(old), out=out)


pins(); inp = ST.StressInputs(ns()); dev0 = tuple(L.DEV); L.SITECFG[V.OLD] = (8, 1); L.DEV = dev0 + (V.OLD,); L.M = 22
agg = {k: {} for k in stored}
for site, fk in (("bbbc041_all", "bbbc"), ("mpidb_full", "mh"), ("lacuna_B", "lh"), (V.OLD, "bbbc")):
    rows = [x for sd in L.SEEDS for x in L.replay(inp, sd, site)[0]]
    for eps, _ in L.SETTINGS:
        for rule in ("S3", "B1u_c0", "B1u_c1s"):
            R = [x for x in rows if abs(x["eps"] - eps) < 1e-12 and x["rule"] == rule]; fpr = np.array([x["fpr"] for x in R])
            agg[fk][ST.PATTERN[site].format(m=22, rule=rule, eps=eps)] = dict(cov=float((fpr <= eps).mean()), fpr_mean=float(fpr.mean()), cannot_certify=float(1-np.mean([x["cert"] for x in R])),
                                                                             ds_eps=float(np.mean([x["dse"] for x in R])), ds_own=float(np.mean([x["dso"] for x in R])), n=len(R))
for k, p in stored.items(): json.dump({"agg": agg[k]}, open(p, "w"))
pins(); inp = ST.StressInputs(ns())
key = lambda R: [(x["site"], x["seed"], x["draw"], x["eps"], x["cand"], x["cert"], round(x["fpr"], 12)) for x in R]
sites = ("bbbc041_all", V.OLD)

# (2a) K1 disabled -> K1vB == S3 at delta/6
_k1 = V.k1_ub_alpha; V.k1_ub_alpha = lambda *a, **k: np.inf
vb_a = [x for s in sites for sd in L.SEEDS for x in V.run_site_seed(inp, sd, s) if x["rule"] == "K1vB"]; V.k1_ub_alpha = _k1
_b = L.b1u; L.b1u = lambda x, n, d: _b(x, n, 0.10/6 if abs(d - 0.10/3) < 1e-15 else d)
ref_a = [x for s in sites for sd in L.SEEDS for x in L.replay(inp, sd, s)[0] if x["rule"] == "S3"]; L.b1u = _b
ok_a = key(vb_a) == key(ref_a)
# (2b) b1u disabled at delta/6 -> K1vB == K1 at delta/6
L.b1u = lambda x, n, d: 1.0 if abs(d - 0.10/6) < 1e-15 else _b(x, n, d); keepA = V.A_K1; V.A_K1 = V.A_VB
rows_b = [x for s in sites for sd in L.SEEDS for x in V.run_site_seed(inp, sd, s)]; L.b1u = _b; V.A_K1 = keepA
ok_b = key([x for x in rows_b if x["rule"] == "K1vB"]) == key([x for x in rows_b if x["rule"] == "K1"])
report("(2) OR rule vs independent references", ok_a and ok_b and len(vb_a) == 2*5*6*2, (ok_a, ok_b))

# (3) copied loop vs m3_k12_dev_eval.run_site_seed (K1 arm)
E.K_ARMS = ("K1",)
mine = [x for s in ("bbbc041_all", "mpidb_full") for sd in L.SEEDS for x in V.run_site_seed(inp, sd, s) if x["rule"] == "K1"]
theirs = [x for s in ("bbbc041_all", "mpidb_full") for sd in L.SEEDS for x in E.run_site_seed(inp, sd, s) if x["rule"] == "K1"]
report("(3) K1 rows equal the approved evaluation loop", key(mine) == key(theirs) and len(mine) == 2*5*6*2)
L.DEV = dev0; del L.SITECFG[V.OLD]


def run(a):
    sys.argv = ["vb"] + sum([[f"--{k}", v] for k, v in vars(a).items()], [])
    keep = (tuple(L.DEV), dict(L.SITECFG), L.M)
    try: V.main(); return None
    except (AssertionError, SystemExit, FileNotFoundError) as ex: return repr(ex)[:180]
    finally: L.DEV, L.M = keep[0], keep[2]; L.SITECFG.clear(); L.SITECFG.update(keep[1])


# (4)
e4 = run(ns(out=str(T / "o4" / "vb.json"))); Q = json.loads((T / "o4" / "vb.json").read_text()) if e4 is None else {}
jd = Q.get("judgement", {})
report("(4) main", e4 is None and len(Q["cells"]) == 32 and all(v.startswith("self-check passed") for v in Q["self_check"].values()) and len(jd["gate"]) == 8
       and "passes" in jd and "old_not_worse" in jd, e4 or (jd["passes"], round(jd["mean_gain_rebuilt_eps010"], 4), jd["old_not_worse"]))
# (5)
o5 = T / "old_missing"; o5.mkdir(); [shutil.copy(q, o5 / q.name) for q in old.iterdir() if "s45" not in q.name]
(o5 / "SHA256SUMS").write_text("".join(l + "\n" for l in (old / "SHA256SUMS").read_text().splitlines() if "s45" not in l))
e5a = run(ns(old_dir=str(o5), out=str(T / "o5a" / "vb.json")))
e5b = run(ns(out=str(T / "o4" / "vb.json")))
ag = json.load(open(stored["bbbc"])); kk = ST.PATTERN[V.OLD].format(m=22, rule="S3", eps=0.1); ag["agg"][kk]["cov"] += 1e-6; bad = T / "rev1_bbbc_bad.json"; json.dump(ag, open(bad, "w"))
keep = L.PIN["rev1_bbbc"]; L.PIN["rev1_bbbc"] = h(bad); e5c = run(ns(rev1_bbbc=str(bad), out=str(T / "o5c" / "vb.json"))); L.PIN["rev1_bbbc"] = keep
report("(5) refusals, nothing written", e5a is not None and e5b and "exists" in e5b and e5c and "self-check" in e5c
       and not (T / "o5a").exists() and not (T / "o5c").exists(), " | ".join(str(v)[:50] for v in (e5a, e5b, e5c)))
print("ALL PASS" if ok else "FAIL")

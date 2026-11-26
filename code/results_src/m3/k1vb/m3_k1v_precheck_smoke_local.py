# Local smoke for m3_k1v_precheck.py.  SYNTHETIC sites (3 rebuilt-style + old), synthetic NIH-validation files, synthetic stored "REV-1" aggregates (from the same
# replay); only real file read: results_src/mh/mh_rev1.py (b1u check).  5 seeds x 4 draws (this process only).  WRITTEN, NOT RUN, until approved.
# (1) platt() recovers known (a, b) = (1.5, -2.0) within 0.1 on 20,000 simulated cells.
# (2) field_stats(): identity e_f = -P_f (TPR_f - eps) holds (max abs < 1e-9); with the oracle q = label, e + v = 0 exactly (ratio 0.0).
# (3) load_nihval(): only key "logits" (score files) and "y" (ids) are accessed; a file whose sha differs from its pin -> refused; a path containing "iml" -> refused.
# (4) main writes 4 sites x 2 candidates x 2 eps = 16 cells, Platt parameters for 5 seeds, old_mpidb_condition_met as a bool, identity_max_abs < 1e-9 everywhere.
# (5) output exists -> refused, nothing overwritten.
# Usage: python -I -B m3_k1v_precheck_smoke_local.py <empty scratch dir>
import argparse, hashlib, importlib.util, json, sys
from pathlib import Path
import numpy as np
from scipy.special import expit

HERE = Path(__file__).resolve().parent; T = Path(sys.argv[1]); assert T.is_dir() and not any(T.iterdir())
MH_REV1 = HERE.parents[1] / "mh" / "mh_rev1.py"
sys.path.insert(0, str(HERE.parents[0] / "k2prep")); sys.path.insert(0, str(HERE.parents[0] / "stress"))
import m3_replay_lib as L   # noqa: E402
L.SEEDS = [41, 42, 43, 44, 45]; L.NDRAW = 4
spec = importlib.util.spec_from_file_location("kv", HERE / "m3_k1v_precheck.py"); K = importlib.util.module_from_spec(spec); spec.loader.exec_module(K)
ST = K.ST
ok = True


def report(tag, cond, extra=""):
    global ok
    ok &= bool(cond); print(tag, "ok" if cond else "FAIL", extra, flush=True)


def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


# (1)
r = np.random.default_rng(2); s = r.normal(0, 2, 20000); y = (r.random(20000) < expit(1.5*s - 2.0)).astype(int); a_, b_ = K.platt(s, y)
report("(1) platt recovers (1.5, -2.0)", abs(a_ - 1.5) < 0.1 and abs(b_ + 2.0) < 0.1, (round(a_, 3), round(b_, 3)))

# synthetic sites
rng = np.random.default_rng(19); lk = T / "a2_link"; lk.mkdir(); mhlh = T / "frozen"; mhlh.mkdir(); old = T / "old"; old.mkdir(); sums = []; ids = {}; man = {}
for site in L.DEV + ("old",):
    nf, nc = 30, 40; g = np.repeat(np.arange(nf), nc); y = (rng.random(nf*nc) < (0.3 if site == "old" else 0.05)).astype(int); mu = rng.normal(0, 1.0, nf)[g]
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
for sd in L.SEEDS:
    sv = rng.normal(0, 1, 800) + 3*yv; np.savez(nv / f"a2ref_s{sd}_nih_validation_logits.npz", logits=np.c_[np.zeros(800), sv])
K.NIHVAL_PIN.update({sd: h(nv / f"a2ref_s{sd}_nih_validation_logits.npz") for sd in L.SEEDS}); K.NIHVAL_PIN["ids"] = h(nv / "ids_nih_validation.npz")
stored = {k: T / f"rev1_{k}.json" for k in ("bbbc", "mh", "lh")}
for p in stored.values(): json.dump({"agg": {}}, open(p, "w"))


def pins():
    L.PIN.update(manifest=h(lk / "MANIFEST.tsv"), sums=h(mhlh / "SHA256SUMS"), ids_iter1=h(T / "ids_iter1.npz"), mh_rev1_py=h(MH_REV1),
                 rev1_bbbc=h(stored["bbbc"]), rev1_mh=h(stored["mh"]), rev1_lh=h(stored["lh"]))


def ns(out="x", nihval=None):
    return argparse.Namespace(a2link=str(lk), ids_iter1=str(T / "ids_iter1.npz"), mhlh=str(mhlh), rev1_bbbc=str(stored["bbbc"]), rev1_mh=str(stored["mh"]),
                              rev1_lh=str(stored["lh"]), mh_rev1_py=str(MH_REV1), old_s41=str(T / "iter1_a2_s41_scores.npz"), old_dir=str(old),
                              nihval_dir=nihval or str(nv), out=out)


pins(); inp = ST.StressInputs(ns()); dev0 = tuple(L.DEV); L.SITECFG[K.OLD] = (8, 1); L.DEV = dev0 + (K.OLD,); L.M = 22
agg = {k: {} for k in stored}
for site, fk in (("bbbc041_all", "bbbc"), ("mpidb_full", "mh"), ("lacuna_B", "lh"), (K.OLD, "bbbc")):
    rows = [x for sd in L.SEEDS for x in L.replay(inp, sd, site)[0]]
    for eps, _ in L.SETTINGS:
        for rule in ("S3", "B1u_c0", "B1u_c1s"):
            fpr = np.array([x["fpr"] for x in rows if abs(x["eps"] - eps) < 1e-12 and x["rule"] == rule])
            agg[fk][ST.PATTERN[site].format(m=22, rule=rule, eps=eps)] = dict(fpr_median=float(np.median(fpr)), cov=float((fpr <= eps).mean()))
for k, p in stored.items(): json.dump({"agg": agg[k]}, open(p, "w"))
L.DEV = dev0; del L.SITECFG[K.OLD]; pins(); inp = ST.StressInputs(ns())

# (2)
pop = L.population(inp, 41, "mpidb_full"); t = float(np.quantile(pop["sneg"]["c0"], 0.92))
fs = K.field_stats(pop, "c0", t, 0.10, expit(pop["s"])); fo = K.field_stats(pop, "c0", t, 0.10, pop["y"].astype(float))
report("(2) identity and oracle q", fs["identity_max_abs"] < 1e-9 and fo["var_ratio_e_plus_v"] is not None and fo["var_ratio_e_plus_v"] < 1e-20,
       (fs["identity_max_abs"], fo["var_ratio_e_plus_v"]))

# (3)
ACC = []; _load = np.load


class Proxy:
    def __init__(self, z): self.z = z
    def __getitem__(self, k): ACC.append(k); return self.z[k]
    def __getattr__(self, k): ACC.append("ATTR:" + k); return getattr(self.z, k)


K.np.load = lambda p, *a, **k: Proxy(_load(p, *a, **k)) if "a2ref" in str(p) else _load(p, *a, **k)
try: K.load_nihval(nv)
finally: K.np.load = _load
keys_ok = sorted(set(ACC)) == ["logits", "y"] and ACC.count("y") == 1 and ACC.count("logits") == 5
keep = K.NIHVAL_PIN[43]; K.NIHVAL_PIN[43] = "0"*64
try: K.load_nihval(nv); e3a = None
except AssertionError as ex: e3a = repr(ex)
K.NIHVAL_PIN[43] = keep
try: K.load_nihval(T / "iml_a2ref"); e3b = None
except AssertionError as ex: e3b = repr(ex)
report("(3) NIH-validation key access and refusals", keys_ok and e3a and "sha" in e3a and e3b and "confirmatory" in e3b, (sorted(set(ACC)), str(e3a)[:40]))


def run(a):
    sys.argv = ["kv"] + sum([[f"--{k}", v] for k, v in vars(a).items()], [])
    keep_ = (tuple(L.DEV), dict(L.SITECFG))
    try: K.main(); return None
    except (AssertionError, SystemExit) as ex: return repr(ex)[:160]
    finally: L.DEV = keep_[0]; L.SITECFG.clear(); L.SITECFG.update(keep_[1])


# (4)
e4 = run(ns(out=str(T / "o4" / "kv.json"))); Q = json.loads((T / "o4" / "kv.json").read_text()) if e4 is None else {}
report("(4) main", e4 is None and len(Q["cells"]) == 16 and len(Q["platt"]) == 5 and isinstance(Q["old_mpidb_condition_met"], bool)
       and all(c["identity_max_abs"] < 1e-9 for c in Q["cells"]), e4 or Q["old_mpidb_ratios"])
# (5)
e5 = run(ns(out=str(T / "o4" / "kv.json")))
report("(5) output exists refused", e5 and "exists" in e5, e5)
print("ALL PASS" if ok else "FAIL")

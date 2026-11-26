# Local smoke for m3_stress.py.  SYNTHETIC development sites, synthetic old MP-IDB files (iter1-style seed-41 file holding mpidb_full_logits plus NIH-named keys,
# and an extraction dir with seeds 42-45 + SHA256SUMS), synthetic stored "REV-1" aggregates; the only real file read is results_src/mh/mh_rev1.py (b1u check).
# WRITTEN, NOT RUN, until 实验管理 approves.  5 seeds x 6 draws (module attribute in this process only) so that the 5-seed self-check is exercised.
# (1) inputs: seed-41 old file -> only key "mpidb_full_logits" accessed; a pulled file that differs from its SHA256SUMS line -> refused; a seed-41 file whose sha
#     differs from the MANIFEST source_sha256 -> refused; an old path containing "iml" -> refused.
# (2) main with all 5 old seeds: 60 cells (3 sites x 3 m x 3 arms x 2 eps + old 3 arms x 2 eps); stream status: bbbc m 10 and 22, mh and lh m 22, old m 22 "self-check
#     passed"; bbbc m 15 and mh/lh m 10, 15 "new random stream"; S3 and B1u_c1s cells at checked (site, m) equal the stored aggregates.
# (3) risk_cells are exactly the cells with coverage < 0.90.
# (4) without --old_dir (seed 41 only): old cells labelled "fewer than 5 seeds".
# (5) refusals, nothing written: a stored aggregate (bbbc, m 10) changed by 1e-6 -> self-check failure; output exists.
# Usage: python -I -B m3_stress_smoke_local.py <empty scratch dir>
import argparse, hashlib, importlib.util, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; T = Path(sys.argv[1]); assert T.is_dir() and not any(T.iterdir())
MH_REV1 = HERE.parents[1] / "mh" / "mh_rev1.py"
sys.path.insert(0, str(HERE.parents[0] / "k2prep"))
import m3_replay_lib as L   # noqa: E402
L.SEEDS = [41, 42, 43, 44, 45]; L.NDRAW = 6
spec = importlib.util.spec_from_file_location("st", HERE / "m3_stress.py"); S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
ok = True


def report(tag, cond, extra=""):
    global ok
    ok &= bool(cond); print(tag, "ok" if cond else "FAIL", extra, flush=True)


def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


rng = np.random.default_rng(13); lk = T / "a2_link"; lk.mkdir(); mhlh = T / "frozen"; mhlh.mkdir(); old = T / "old"; old.mkdir(); sums = []; ids = {}; man_rows = {}
for site in L.DEV + ("old",):
    nf, nc = 40, 50; g = np.repeat(np.arange(nf), nc); y = (rng.random(nf*nc) < 0.06).astype(int); mu = rng.normal(0, 1.0, nf)[g]
    for sd in L.SEEDS:
        s = mu + rng.normal(0, 1, nf*nc) + 2.5*y; z = np.c_[np.zeros(len(s)), s]
        if site == "bbbc041_all":
            f = lk / f"a2_s{sd}_bbbc041_all_logits.npz"; np.savez(f, logits=z); man_rows.setdefault(sd, {})["file"] = h(f)
        elif site == "old":
            if sd == 41:
                f = T / "iter1_a2_s41_scores.npz"; np.savez(f, mpidb_full_logits=z, nih_locked_test_logits=np.ones((3, 2))); man_rows.setdefault(41, {})["src"] = h(f)
            else:
                f = old / f"a2_s{sd}_mpidb_full_logits.npz"; np.savez(f, logits=z)
        else:
            f = mhlh / f"a2_s{sd}_{L.SHORT[site]}_logits.npz"; np.savez(f, logits=z); sums.append(f"{h(f)}  ./{f.name}")
    if site in L.SHORT: fi = mhlh / f"ids_{L.SHORT[site]}.npz"; np.savez(fi, y=y, field=g); sums.append(f"{h(fi)}  ./{fi.name}")
    ids[site] = (y, g)
(old / "SHA256SUMS").write_text("".join(f"{h(old / f'a2_s{sd}_mpidb_full_logits.npz')}  ./a2_s{sd}_mpidb_full_logits.npz\n" for sd in (42, 43, 44, 45)))
(lk / "MANIFEST.tsv").write_text("source\tsource_sha256\tmode\tfile_sha256\n" + "".join(
    f"/x/iter1/a2_s{sd}_scores.npz[bbbc041_all_logits]\t{man_rows[sd].get('src', 'abc')}\twritten\t{man_rows[sd]['file']}\n" for sd in L.SEEDS))
(mhlh / "SHA256SUMS").write_text("\n".join(sums) + "\n")
np.savez(T / "ids_iter1.npz", bbbc041_all_y=ids["bbbc041_all"][0], bbbc041_all_group=ids["bbbc041_all"][1], mpidb_full_y=ids["old"][0], mpidb_full_group=ids["old"][1],
         nih_locked_test_y=np.ones(3, int))
stored = {k: T / f"rev1_{k}.json" for k in ("bbbc", "mh", "lh")}
for p in stored.values(): json.dump({"agg": {}}, open(p, "w"))


def pins():
    L.PIN.update(manifest=h(lk / "MANIFEST.tsv"), sums=h(mhlh / "SHA256SUMS"), ids_iter1=h(T / "ids_iter1.npz"), mh_rev1_py=h(MH_REV1),
                 rev1_bbbc=h(stored["bbbc"]), rev1_mh=h(stored["mh"]), rev1_lh=h(stored["lh"]))


def ns(old_dir=True, out="x"):
    return argparse.Namespace(a2link=str(lk), ids_iter1=str(T / "ids_iter1.npz"), mhlh=str(mhlh), rev1_bbbc=str(stored["bbbc"]), rev1_mh=str(stored["mh"]),
                              rev1_lh=str(stored["lh"]), mh_rev1_py=str(MH_REV1), old_s41=str(T / "iter1_a2_s41_scores.npz"), old_dir=(str(old) if old_dir else None), out=out)


# build synthetic stored aggregates from the same replay (bbbc m 10 and 22; mh, lh m 22; old m 22 inside the bbbc file)
pins(); inp = S.StressInputs(ns()); L.SITECFG[S.OLD] = (8, 1); dev0 = tuple(L.DEV); L.DEV = dev0 + (S.OLD,)
want = [("bbbc041_all", 10, "bbbc"), ("bbbc041_all", 22, "bbbc"), ("mpidb_full", 22, "mh"), ("lacuna_B", 22, "lh"), (S.OLD, 22, "bbbc")]
agg = {k: {} for k in stored}
for site, m, fk in want:
    L.M = m; rows = [r for sd in L.SEEDS for r in L.replay(inp, sd, site)[0]]
    for eps, _ in L.SETTINGS:
        for rule in ("S3", "B1u_c0", "B1u_c1s"):
            R = [r for r in rows if abs(r["eps"] - eps) < 1e-12 and r["rule"] == rule]; fpr = np.array([r["fpr"] for r in R])
            agg[fk][S.PATTERN[site].format(m=m, rule=rule, eps=eps)] = dict(cov=float((fpr <= eps).mean()), fpr_mean=float(fpr.mean()), cannot_certify=float(1-np.mean([r["cert"] for r in R])),
                                                                            ds_eps=float(np.mean([r["dse"] for r in R])), ds_own=float(np.mean([r["dso"] for r in R])), n=len(R))
L.DEV = dev0; del L.SITECFG[S.OLD]; L.M = 22
for k, p in stored.items(): json.dump({"agg": agg[k]}, open(p, "w"))
pins()

# (1) inputs
ACC = []; _load = np.load


class Proxy:
    def __init__(self, z): self.z = z
    def __getitem__(self, k): ACC.append(k); return self.z[k]
    def __getattr__(self, k): ACC.append("ATTR:" + k); return getattr(self.z, k)


inp = S.StressInputs(ns()); S.np.load = lambda p, *a, **k: (Proxy(_load(p, *a, **k)) if Path(p).name == "iter1_a2_s41_scores.npz" else _load(p, *a, **k))
try: inp.load_site(41, S.OLD)
finally: S.np.load = _load
errs = []
bad_old = T / "old_bad"; bad_old.mkdir(); [Path(bad_old / q.name).write_bytes(q.read_bytes()) for q in old.iterdir()]; np.savez(bad_old / "a2_s43_mpidb_full_logits.npz", logits=np.zeros((4, 2)))
for kw in (dict(old_dir=str(bad_old)), dict(old_s41=str(T / "ids_iter1.npz")), dict(old_dir=str(T / "iml_old"))):
    a = ns(); [setattr(a, k, v) for k, v in kw.items()]
    try: S.StressInputs(a); errs.append(None)
    except (AssertionError, FileNotFoundError) as e: errs.append(repr(e)[:80])
report("(1) inputs", ACC == ["mpidb_full_logits"] and errs[0] and "sha" in errs[0] and errs[1] and "sha" in errs[1] and errs[2] and "confirmatory" in errs[2], (ACC, errs))


def run(out, old_dir=True, rev1_bbbc=None):
    a = ns(old_dir=old_dir, out=str(out))
    if rev1_bbbc: a.rev1_bbbc = rev1_bbbc
    sys.argv = ["st"] + sum([[f"--{k}", v] for k, v in vars(a).items() if v is not None], [])
    keep = (tuple(L.DEV), dict(L.SITECFG), list(L.SEEDS), L.M)
    try: S.main(); return None
    except (AssertionError, SystemExit) as e: return repr(e)[:180]
    finally: L.DEV, L.SEEDS, L.M = keep[0], keep[2], keep[3]; L.SITECFG.clear(); L.SITECFG.update(keep[1])


# (2) + (3)
e2 = run(T / "o2" / "st.json"); Q = json.loads((T / "o2" / "st.json").read_text()) if e2 is None else {}
cells = Q.get("cells", []); st = {(c["site"], c["m"]): c["stream"] for c in cells}
passed = {("bbbc041_all", 10), ("bbbc041_all", 22), ("mpidb_full", 22), ("lacuna_B", 22), (S.OLD, 22)}
stream_ok = all(st[k].startswith("self-check passed") for k in passed) and all(v.startswith("new random stream") for k, v in st.items() if k not in passed)
ref_ok = True
for c in cells:
    if (c["site"], c["m"]) in passed and c["arm"] in ("S3", "B1u_c1s"):
        fk = {"bbbc041_all": "bbbc", "mpidb_full": "mh", "lacuna_B": "lh", S.OLD: "bbbc"}[c["site"]]
        sv = agg[fk][S.PATTERN[c["site"]].format(m=c["m"], rule=c["arm"], eps=c["eps"])]; ref_ok &= abs(c["coverage"] - sv["cov"]) < 1e-12 and abs(c["ds_own"] - sv["ds_own"]) < 1e-12
risk_ok = e2 is None and sorted((r["site"], r["m"], r["arm"], r["eps"]) for r in Q["risk_cells"]) == sorted((c["site"], c["m"], c["arm"], c["eps"]) for c in cells if c["coverage"] < 0.90)
report("(2) main, streams, references", e2 is None and len(cells) == 60 and stream_ok and ref_ok, e2 or sorted(set(st.values())))
report("(3) risk cells", risk_ok, len(Q.get("risk_cells", [])))
# (4) seed 41 only
e4 = run(T / "o4" / "st.json", old_dir=False); Q4 = json.loads((T / "o4" / "st.json").read_text()) if e4 is None else {}
oc = [c for c in Q4.get("cells", []) if c["site"] == S.OLD]
report("(4) seed 41 only", e4 is None and len(oc) == 6 and all(c["stream"].startswith("fewer than 5 seeds (1)") for c in oc) and Q4["old_mpidb_seeds"] == [41], e4)
# (5) refusals
ag = json.load(open(stored["bbbc"])); kk = S.PATTERN["bbbc041_all"].format(m=10, rule="S3", eps=0.1); ag["agg"][kk]["ds_own"] += 1e-6; bad = T / "rev1_bbbc_bad.json"; json.dump(ag, open(bad, "w"))
keep = L.PIN["rev1_bbbc"]; L.PIN["rev1_bbbc"] = h(bad); e5 = run(T / "o5" / "st.json", rev1_bbbc=str(bad)); L.PIN["rev1_bbbc"] = keep
e6 = run(T / "o2" / "st.json")
report("(5) refusals, nothing written", e5 and "self-check" in e5 and not (T / "o5").exists() and e6 and "exists" in e6, (str(e5)[:60], str(e6)[:40]))
print("ALL PASS" if ok else "FAIL")

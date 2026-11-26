# P5a D-d (实验管理 2026-11-25, decision 1): per-seed coverage and Delta Sens_eps of every rule of the development-site forest plot F-d, and their two-sided 90 %
#   seed-cluster intervals (ebatch_lib.seed_cluster_ci, 359b6bce).  WRITTEN, NOT RUN, until 实验管理 reviews it.  Development sites only (BBBC-ALL, MP-IDB, Lacuna-B), m = 22.
# Nature: post hoc descriptive intervals for a figure; no criterion, no selection.  Nothing is trained.  Four groups of rules:
#   npz    S3, R1_2, B1, B1u_c0, B1u_c1s, R1_all: the stored REV-1 per-draw files (runs/c5b_eval_20261027T044702Z/eval_a2/<s>_rev1/c5b_<s>_rev1_perdraw_a2_m22.npz),
#          1000 rows = 5 seeds x 200 draws in seed order; no rule is run.  Delta Sens_eps per draw = sens - orac_unc[eps] of the seed (m3_replay_lib.population).
#   w3     CB_c1s, GEE_c1s, BB_c1s: deterministic replay of the REV-1 draws (m3_replay_lib.replay) with the W3 bounds (m3_w3_baselines.make_rules, 249be397) as extra rules;
#          the bootstrap weights of CB come from W3's own stream (default_rng([20261031, site, seed, draw, eps, candidate])), independent of the other rules.
#   split  SPLIT6 (Split-Select): m3_split_dev_eval.run_site_seed (cbf1e7f4) with the frozen Platt of the seed (m3_k1v_precheck).
#   ssrcps SSRCPS_WSR, SSRCPS_CLT: m3_ssrcps_dev_eval.run_new_arms (d252ad3a) with the seed's Youden source threshold and Platt.
# Anchors (exit 3, nothing written; tolerance 1e-12):
#   every module and input pinned; npz: the 5-seed aggregate cov, fpr_mean, ds_eps, cannot_certify equal the stored REV-1 aggregates (inp.rev1, keys of m3_replay_lib.STORED_KEY),
#   and the per-seed cov / dsens of S3, R1_2, B1u_c0, B1u_c1s equal g3_intervals.json per_seed (cd9543cf); w3: coverage, fpr_mean, ds_eps, cannot_certify equal w3.json (2f015508);
#   split / ssrcps: the same four fields equal ssrcps.json (3ab40370); in every group n == 1000 per (rule, site, eps).
# Resources (实验管理): at most 4 processes (--workers <= 4), nice 19 and one BLAS thread in run.sh; MemAvailable > 20 GB is checked here before anything starts.
# Modes: --timing <group>:<site>:<seed> runs one task, prints its wall time and per-seed values, and writes nothing (used first for CB, as 实验管理 asked);
#        otherwise every task of --groups (default all four) is run and one NEW directory --out is written ('x' mode): d_d_dev_intervals.json, SHA256SUMS.
# Run (VM, pyproject env): python -I -B d_d_dev_intervals.py --m3 <zdeploy results_src/m3> --ebatch_lib <ebatch_lib.py> --g3 <g3_intervals.json> --w3 <w3.json>
#   --ssrcps <ssrcps.json> --perdraw_dir <runs/c5b_eval_20261027T044702Z/eval_a2> --nihval_dir <esparse_inputs_.../a2ref> + the 7 Inputs arguments of m3_replay_lib
import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[_k] = "1"
import argparse, hashlib, importlib.util, json, re, sys, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np

PINS = {"m3_replay_lib.py": "ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de", "m3_w3_baselines.py": "249be3977f9c7e06dbcaea0830baef92a0fb65b7430d7c802a457cb0c809c20e",
        "m3_split_dev_eval.py": "cbf1e7f4f33dd1c74cbfd0f3ac4d34a397e8a7ef44b82c9365ca117733bf6457", "m3_ssrcps_dev_eval.py": "d252ad3ad4db19eec827d5ed223772f6edb64e9e456885670695702129d402f8",
        "m3_k12_dev_eval.py": "0091a3337074a227e404f0585966b87a9f5e9e55cfbb8dd766df7e4618893b67", "m3_k1v_precheck.py": "0e6aa9f8f423dbaed286db969afa8a14c9ce477e37bba962732d3b411002147c",
        "m3_stress.py": "75eb9afd3a130b00170757c4027f8b158bf95deb0ad6124229e9292c60e2d1ad", "m3_k1vb_dev_eval.py": "7b7e6ecda22b789c573451c757b6e197d747f77e10bcfc201ccdb85921e25abd",
        "ebatch_lib.py": "359b6bce2dfdbf50a2acb8b909d631ab40d445447275915e2318e70981ef007f", "g3": "cd9543cfc0cfe91b583b28431191f14858cda2b712f8f39b49b4539e33b58b75",
        "w3": "2f015508424452d276f6b570cf89c29d2eb0b4d36637140e736efa87a29f5776", "ssrcps": "3ab403706d8417a0e438831c7c9e89f257c6b8b055c37d0a98ee18360327c1b8",
        "bbbc041_all": "1910e656bc62bc9ee75e2816fd31dd8390c0e71a7ba4ca4d4e2dd887b91144d0", "mpidb_full": "69c51a383a963e75ce9c1ba950bc5eb1ed931b8604e1525d8e7ae40af6f941d0",
        "lacuna_B": "00cec64b41a69dfa7ce5ace9dcbcb2d4b282950f664157bf18af3f4585464654"}
NPZ_NAME = {"bbbc041_all": "bbbc_rev1/c5b_bbbc_rev1_perdraw_a2_m22.npz", "mpidb_full": "mh_rev1/c5b_mh_rev1_perdraw_a2_m22.npz", "lacuna_B": "lh_rev1/c5b_lh_rev1_perdraw_a2_m22.npz"}
GROUPS = {"npz": ("S3", "R1_2", "B1", "B1u_c0", "B1u_c1s", "R1_all"), "w3": ("CB_c1s", "GEE_c1s", "BB_c1s"), "split": ("SPLIT6",), "ssrcps": ("SSRCPS_WSR", "SSRCPS_CLT")}
G3_RULES = ("S3", "R1_2", "B1u_c0", "B1u_c1s")
SITES = ("bbbc041_all", "mpidb_full", "lacuna_B"); SEEDS = (41, 42, 43, 44, 45); EPSS = (0.1, 0.05); NDRAW = 200; TOL = 1e-12
FORBID = re.compile(r"plasmo|iml|r8|m5", re.I)
A = None; M = {}


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def mem_available_gb():
    for l in open("/proc/meminfo"):
        if l.startswith("MemAvailable:"): return int(l.split()[1]) / 1024 ** 2
    return 0.0


def load_modules(m3, ebl):
    """Import the frozen modules from the zdeploy layout (their own sys.path rules) after checking every file's sha."""
    files = {"m3_replay_lib.py": m3 / "k2prep" / "m3_replay_lib.py", "m3_w3_baselines.py": m3 / "k2prep" / "m3_w3_baselines.py", "m3_k12_dev_eval.py": m3 / "k2prep" / "m3_k12_dev_eval.py",
             "m3_split_dev_eval.py": m3 / "split" / "m3_split_dev_eval.py", "m3_ssrcps_dev_eval.py": m3 / "ssrcps" / "m3_ssrcps_dev_eval.py",
             "m3_k1v_precheck.py": m3 / "k1vb" / "m3_k1v_precheck.py", "m3_k1vb_dev_eval.py": m3 / "k1vb" / "m3_k1vb_dev_eval.py", "m3_stress.py": m3 / "stress" / "m3_stress.py", "ebatch_lib.py": ebl}
    for n, p in files.items():
        if sha(p) != PINS[n]: raise ValueError(f"module sha differs: {n}")
    for sub in ("k2prep", "stress", "k1vb", "split", "ssrcps"): sys.path.insert(0, str(m3 / sub))
    sys.path.insert(0, str(ebl.parent))
    import m3_replay_lib as L, m3_w3_baselines as W3, m3_split_dev_eval as SP, m3_ssrcps_dev_eval as SS, m3_k1v_precheck as KV, ebatch_lib as EB
    for mod, n in ((L, "m3_replay_lib.py"), (W3, "m3_w3_baselines.py"), (SP, "m3_split_dev_eval.py"), (SS, "m3_ssrcps_dev_eval.py"), (KV, "m3_k1v_precheck.py"), (EB, "ebatch_lib.py")):
        if Path(mod.__file__).resolve() != files[n].resolve(): raise ValueError(f"imported {n} from {mod.__file__}, not from the pinned path")
    return dict(L=L, W3=W3, SP=SP, SS=SS, KV=KV, EB=EB)


def inputs():
    return M["L"].Inputs(A.a2link, A.ids_iter1, A.mhlh, A.rev1_bbbc, A.rev1_mh, A.rev1_lh, A.mh_rev1_py)


def per_seed(rows_of_seed):
    """rows of one seed for one rule -> {eps: (cov, ds_eps, fpr_mean, cannot_certify, n)} over its draws."""
    out = {}
    for e in EPSS:
        R = [r for r in rows_of_seed if abs(r["eps"] - e) < 1e-12]
        fpr = np.array([r["fpr"] for r in R]); out[e] = (float((fpr <= e).mean()), float(np.mean([r["dse"] for r in R])), float(fpr.mean()),
                                                         float(1 - np.mean([r["cert"] for r in R])), len(R))
    return out


def task(t):
    """(group, site, seed) -> (key, {rule: {eps: (cov, ds, fpr_mean, cc, n)}}, seconds)."""
    g, site, sd = t; L = M["L"]; t0 = time.time(); inp = inputs(); res = {}
    if g == "npz":
        z = np.load(Path(A.perdraw_dir) / NPZ_NAME[site]); i = SEEDS.index(sd); sl = slice(i * NDRAW, (i + 1) * NDRAW); pop = L.population(inp, sd, site)
        for r in GROUPS["npz"]:
            rows = [dict(eps=e, fpr=f, cert=c, dse=s - pop["orac_unc"][e]) for e in EPSS
                    for f, c, s in zip(z[f"{r}|{e}|fpr"][sl], z[f"{r}|{e}|cert"][sl], z[f"{r}|{e}|sens"][sl])]
            res[r] = per_seed(rows)
    elif g == "w3":
        rules = {k: v for k, v in M["W3"].make_rules().items() if k in GROUPS["w3"]}
        rows, _ = L.replay(inp, sd, site, extra=rules)
        for r in GROUPS["w3"]: res[r] = per_seed([x for x in rows if x["rule"] == r])
    else:
        KV, SS = M["KV"], M["SS"]; nv = KV.load_nihval(A.nihval_dir); platt = KV.platt(*nv[sd])
        if g == "split":
            rows = M["SP"].run_site_seed(inp, sd, site, platt)
        else:
            thr = SS.youden(nv[sd][1], nv[sd][0]); rows, _ = SS.run_new_arms(inp, sd, site, thr, platt)
        for r in GROUPS[g]: res[r] = per_seed([x for x in rows if x["rule"] == r])
    return (g, site, sd), res, time.time() - t0


def init(a):
    global A, M
    A = a; M = load_modules(Path(a.m3), Path(a.ebatch_lib))


def main():
    ap = argparse.ArgumentParser()
    for k in ("m3", "ebatch_lib", "g3", "w3", "ssrcps", "perdraw_dir", "nihval_dir", "a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py"): ap.add_argument("--" + k, required=True)
    ap.add_argument("--out"); ap.add_argument("--timing"); ap.add_argument("--groups", default="npz,w3,split,ssrcps"); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if any(FORBID.search(str(v)) for v in vars(a).values()): print("GATE (exit 4): forbidden argument", file=sys.stderr); return 4
    if not 1 <= a.workers <= 4: return fail("--workers must be 1 ... 4")
    if (a.out is None) == (a.timing is None): return fail("give exactly one of --out and --timing")
    if a.out and Path(a.out).exists(): return fail("output directory exists", a.out)
    mem = mem_available_gb()
    if mem <= 20: return fail(f"MemAvailable {mem:.1f} GB <= 20 GB")
    try:
        init(a)
        for k in ("g3", "w3", "ssrcps"):
            if sha(getattr(a, k)) != PINS[k]: raise ValueError(f"input sha differs: {k}")
        for s in SITES:
            if sha(Path(a.perdraw_dir) / NPZ_NAME[s]) != PINS[s]: raise ValueError(f"per-draw file sha differs: {s}")
    except (ValueError, AssertionError, OSError) as e: return fail(e)
    if a.timing:
        g, site, sd = a.timing.split(":"); t = (g, site, int(sd))
        if g not in GROUPS or site not in SITES or t[2] not in SEEDS: return fail("bad --timing", a.timing)
        _, res, sec = task(t); print(json.dumps(dict(task=a.timing, seconds=sec, mem_available_gb_at_start=mem, per_seed=res), indent=1)); return 0
    groups = [g for g in a.groups.split(",") if g]
    if any(g not in GROUPS for g in groups): return fail("unknown group", groups)
    tasks = [(g, s, sd) for g in groups for s in SITES for sd in SEEDS]
    with get_context("fork").Pool(a.workers, initializer=init, initargs=(a,)) as pool: out = pool.map(task, tasks, chunksize=1)
    per = {}; secs = {}
    for (g, s, sd), res, sec in out:
        secs[f"{g}|{s}|{sd}"] = sec
        for r, byeps in res.items():
            for e, v in byeps.items(): per.setdefault((r, s, e), {})[sd] = v
    # anchors
    inp = inputs(); L = M["L"]; bad = []
    g3 = json.load(open(a.g3))["cells"]; w3 = {(c["rule"], c["site"], round(c["eps"], 6)): c for c in json.load(open(a.w3))["cells"]}
    ss = {(c["arm"], c["site"], round(c["eps"], 6)): c for c in json.load(open(a.ssrcps))["cells"]}
    cells = []
    for (r, s, e), bysd in sorted(per.items()):
        if sorted(bysd) != list(SEEDS): bad.append(("seeds", r, s, e)); continue
        v = [bysd[sd] for sd in SEEDS]; n = sum(x[4] for x in v)
        agg = dict(coverage=float(np.mean([x[0] for x in v])), ds_eps=float(np.mean([x[1] for x in v])), fpr_mean=float(np.mean([x[2] for x in v])),
                   cannot_certify=float(np.mean([x[3] for x in v])), n=n)
        if n != len(SEEDS) * NDRAW: bad.append(("n", r, s, e, n))
        g = next(k for k, rs in GROUPS.items() if r in rs)
        if g == "npz":
            st = json.load(open(inp.rev1[s]))["agg"][L.STORED_KEY[s].format(rule=r, eps=e)]
            ref = dict(coverage=st["cov"], ds_eps=st["ds_eps"], fpr_mean=st["fpr_mean"], cannot_certify=st["cannot_certify"])
            if r in G3_RULES:
                for kind, idx in (("cov", 0), ("dsens", 1)):
                    ps = g3[f"{s}|{r}|{e}|{kind}"]["per_seed"]
                    if any(abs(x[idx] - y) > TOL for x, y in zip(v, ps)): bad.append(("g3 per_seed", r, s, e, kind))
        elif g == "w3": c = w3[(r, s, round(e, 6))]; ref = {k: c[k] for k in ("coverage", "ds_eps", "fpr_mean", "cannot_certify")}
        else: c = ss[(r, s, round(e, 6))]; ref = {k: c[k] for k in ("coverage", "ds_eps", "fpr_mean", "cannot_certify")}
        bad += [(r, s, e, k, agg[k], ref[k]) for k in ref if abs(agg[k] - ref[k]) > TOL]
        ci_c = M["EB"].seed_cluster_ci([x[0] for x in v]); ci_d = M["EB"].seed_cluster_ci([x[1] for x in v])
        cells.append(dict(rule=r, site=s, eps=e, group=g, **agg, per_seed_coverage=[x[0] for x in v], per_seed_ds_eps=[x[1] for x in v],
                          coverage_ci90_seedcluster=ci_c, ds_eps_ci90_seedcluster=ci_d))
    if bad: return fail("anchors", bad[:10])
    rep = dict(nature="P5a D-d: post hoc descriptive two-sided 90% seed-cluster intervals (ebatch_lib.seed_cluster_ci) of coverage and Delta Sens_eps on the development sites; "
               "no criterion", script_sha256=sha(__file__), pins=PINS, score_files_sha256=inp.used, groups=groups, workers=a.workers, seconds=secs, cells=cells)
    payload = json.dumps(rep, indent=1, allow_nan=False); out = Path(a.out); out.mkdir(parents=True)
    with open(out / "d_d_dev_intervals.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh: fh.write(f"{sha(out / 'd_d_dev_intervals.json')}  d_d_dev_intervals.json\n")
    print("written", out, "| cells", len(cells), "| total task seconds", round(sum(secs.values()))); return 0


if __name__ == "__main__": sys.exit(main())

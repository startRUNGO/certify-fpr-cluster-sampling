# P5a SIMW3: the three cluster-aware baselines of Section 6 (cluster bootstrap CB, GEE cluster-robust GEE, beta-binomial BB) added to the simulation study of Section 5 (实验管理, 2026-11-24, after the BJ benchmark).
# WRITTEN, NOT RUN, until 实验管理 reviews it.  Synthetic only: no project data, no image, no classifier.
# NATURE: comparison baselines in the simulation; no rule, no constant and no baseline definition is changed; nothing is selected; every cell is reported.
# Baseline definitions: the functions gee_ub, cb_ub and bb_ub of results_src/m3/k2prep/m3_w3_baselines.py (249be397, the W3 run behind the Cluster bootstrap / GEE cluster-robust / Beta-binomial rows of Table 3)
#   are imported and called unchanged: single candidate c1s (as Table 3), level delta = 0.10, the REV-1 fixed-sequence test on the same descending grid (m3_replay_lib.fixed_sequence; an exception in a bound
#   gives "not certified", is counted, and a cell with exceptions above 1% of its draws is flagged), CB with B = 2000 field resamples drawn once per (realisation, draw, eps, candidate).
#   The one difference from W3: the CB resampling stream.  W3 used default_rng([20261031, site index, seed, draw, eps code, candidate code]); here the site index is 0 for every synthetic cell, so three
#   elements [13, design code (10 or 12), cell index] are appended to keep the streams of different cells apart.  This stream does not touch the draw stream.
# Data generation, streams, sizes (unchanged): E10 single-factor design (10 cells, e10_synth.cells(False), make_site with default_rng([20261103, 10, cell, seed, 0])) and E12 stress design (24 cells, e12_stress.cells(),
#   make_site2 with e12_stress.stream(cell, seed, 0); the four C2 cells use e12_stress.two_stage_sampler); 5 realisations (seeds 41..45), 200 draws of m = 22 fields per realisation, eps in {0.10, 0.05}.
#   Uniform-sampler cells run m3_replay_lib.replay itself with the baselines as its `extra` rules (exactly as W3); the C2 cells run replay_hook, a copy of replay whose only change is the field-sampler hook.
# Anchors (exit 3, nothing written): pins of every imported module and both anchor files; (b) make_site2 with neutral parameters equals make_site bit for bit; (c) replay_hook with the uniform sampler equals
#   m3_replay_lib.replay row for row on E10 cell 0, seed 41, 5 draws, every rule including the baselines; (d) every S3, B1u_c0 and B1u_c1s row of every E12 cell that was run equals the stored E12 row
#   (e12_rows.json; cert, cand, fpr, sens, dse, dso exactly); (e) for every E10 cell that was run with 5 realisations x 200 draws: S3 and B1u_c1s per-realisation coverage, cannot-certify and ds_own equal
#   e10.json within 1e-12.  (d) and (e) check that the draws and the existing rules are the published ones before the baselines are reported.
# Output per design, cell, rule (S3, B1u_c0, B1u_c1s, CB_c1s, GEE_c1s, BB_c1s), eps: n, coverage, cannot_certify, sens, ds_eps, ds_own, fpr_mean, each with its 5 per-realisation values and its Monte Carlo SE
#   (SD of the per-realisation values / sqrt(5)); coverage also with the seed-cluster 90% interval of E10/E12 (ebatch_lib.seed_cluster_ci) and the binomial SE sqrt(c(1 - c) / n); bound exceptions and flag.
# Cost (one run, defaults): workers processes x 1 BLAS thread (default 16, at most 32); 170 tasks (34 cells x 5 realisations), each one synthetic site; memory about 0.3 GB per worker plus about 2 GB in the parent
#   (stored E12 rows for anchor d and the new rows), so about 7 GB at 16 workers; disk about 80 MB (rows) + 1 MB (summary).  Time: dominated by the BB profile likelihood (W3 note: 1-2 h for 12,000
#   draw x eps x candidate units); here 68,000 units, so about 6-11 CPU hours, about 0.4-0.8 h wall at 16 workers.  A timed pilot (one task) is proposed before the full run.
# Output: one NEW directory named simw3_<UTC yyyymmddThhmmssZ>: simw3_summary.json, simw3_rows.json (per draw, without ni and fields), SHA256SUMS; 'x' mode, allow_nan=False.
# Usage: python3 -I -B simw3_baselines.py --e10_json <ebatch_e10.json> --e12_rows <out_e12/e12_rows.json> --out <parent>/simw3_<UTC> [--workers 16] [--cells e10:0,3;e12:19] [--reps 5] [--ndraw 200]
import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[_k] = "1"
import argparse, hashlib, json, re, sys, time
from collections import Counter, defaultdict
from multiprocessing import get_context
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; R = HERE.parent
sys.path[:0] = [str(R / "e10"), str(R / "e12"), str(R / "e_batch"), str(R / "m3" / "k2prep")]
PINS = {"m3_replay_lib": "ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de", "m3_w3_baselines": "249be3977f9c7e06dbcaea0830baef92a0fb65b7430d7c802a457cb0c809c20e",
        "e10_synth": "00be97b1f15a0fd5ac7c03b06d41f6adc0ac44a2c11b165e54748fbe68f7e08e", "e12_stress": "aef3570877fdf7427fb469d4e99d3341547342bb36fb86366e47c69b971ccd59",
        "ebatch_lib": "359b6bce2dfdbf50a2acb8b909d631ab40d445447275915e2318e70981ef007f"}
FILE_PINS = {"e10_json": "a25bc85c40a9b85c742444deb16b704107deb399f49b1cb15d9c94d6a79e9bc8", "e12_rows": "ece7c768005f0a8b79e0e4a883a1a50dc9f6f70e5b70995e4ece4d9949d54a20"}
EPSS = (0.10, 0.05); M = 22; SITE = "bbbc041_all"; CAND = "c1s"; DESIGN_CODE = {"e10": 10, "e12": 12}; N_CELLS = {"e10": 10, "e12": 24}
ANCHOR_RULES = ("S3", "B1u_c0", "B1u_c1s"); CORE = ("cert", "cand", "fpr", "sens", "dse", "dso"); METRICS = ("coverage", "cannot_certify", "sens", "ds_eps", "ds_own", "fpr_mean")
OUT_RE = re.compile(r"^simw3_\d{8}T\d{6}Z$"); MAX_WORKERS = 32
_MOD = {}


class AnchorError(Exception): pass


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def mods():
    """Imports once per process.  e10_synth sets the thread variables to 8 when imported; numpy is already loaded with 1 thread, and the variables are set back to 1 for any child."""
    if not _MOD:
        import e10_synth as E10, e12_stress as E12, m3_replay_lib as L, m3_w3_baselines as W3, ebatch_lib as EB
        for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[_k] = "1"
        if W3.L is not L: raise AnchorError("m3_w3_baselines did not import the same m3_replay_lib")
        _MOD.update(E10=E10, E12=E12, L=L, W3=W3, EB=EB)
    return _MOD


def check_pins():
    m = mods(); files = {"m3_replay_lib": m["L"], "m3_w3_baselines": m["W3"], "e10_synth": m["E10"], "e12_stress": m["E12"], "ebatch_lib": m["EB"]}
    bad = [n for n, mod in files.items() if sha(mod.__file__) != PINS[n]]
    if bad: raise AnchorError(("module differs from its pin", bad))


class Synth:
    used = {}

    def __init__(self, sites): self.sites = sites

    def load_site(self, sd, st): return self.sites[sd]


def cell_spec(design, ci):
    m = mods()
    return m["E10"].cells(False)[ci] if design == "e10" else m["E12"].cells()[ci]


def site_for(design, ci, sd):
    """The synthetic site of one realisation, exactly as E10 / E12 generate it, and the field sampler (None = uniform, the library draw)."""
    m = mods(); E10, E12 = m["E10"], m["E12"]; cell = cell_spec(design, ci)
    if design == "e10":
        s, y, g = E10.make_site(np.random.default_rng([20261103, 10, ci, sd, 0]), **cell); return (s, y, g), None
    s, y, g, slide_of, _ = E12.make_site2(E10, np.random.default_rng(E12.stream(ci, sd, 0)), cell["tau"], cell["r"], cell["pi"], cell["tail"], cell["F"], cell["size_sigma"], cell["rare"], cell["slides"])
    smp = None if cell["two_stage"] is None else E12.two_stage_sampler(slide_of, cell["slides"][0], cell["two_stage"])
    return (s, y, g), smp


def make_extra(design, ci, exc):
    """The W3 baselines for candidate c1s with the W3 functions; CB weights drawn once per (seed, draw, eps, candidate) and shared along the grid, as in W3."""
    m = mods(); L, W3 = m["L"], m["W3"]; wc = {}

    def weights(sd, d_, eps, k):
        key = (sd, d_, eps, k)
        if key not in wc:
            rng = np.random.default_rng([W3.BOOT_SEED, L.DEV.index(SITE), sd, d_, int(round(eps * 100)), W3.CAND_CODE[k], 13, DESIGN_CODE[design], ci])
            wc.clear(); wc[key] = rng.multinomial(L.M, np.full(L.M, 1 / L.M), size=W3.B).astype(float)
        return wc[key]

    def counted(name, fn):
        def g(x, n, ctx):
            try: return fn(x, n, ctx)
            except Exception as e:
                exc[(ctx[3], name)][type(e).__name__] += 1; raise
        return g
    k = CAND
    return {"CB_" + k: (k, counted("CB_" + k, lambda x, n, ctx: W3.cb_ub(x, n, weights(ctx[1], ctx[2], ctx[3], k)))),
            "GEE_" + k: (k, counted("GEE_" + k, lambda x, n, ctx: W3.gee_ub(x, n))),
            "BB_" + k: (k, counted("BB_" + k, lambda x, n, ctx: W3.bb_ub(x, n)))}


def replay_hook(L, inp, sd, site, extra, sampler):
    """Copy of m3_replay_lib.replay (ed9d1f17).  The only change: the m fields are drawn by sampler(rng, fl, m) instead of rng.choice(len(fl), m, replace=False); with e12_stress.uniform_sampler (the same call)
    it equals the library row for row (anchor c)."""
    pop = L.population(inp, sd, site); fl, fl_neg, fl_posv = pop["fl"], pop["fl_neg"], pop["fl_posv"]; sens_at, fpr_at = pop["sens_at"], pop["fpr_at"]
    key8, sidx = L.SITECFG[site]; ai = L.AI; rows = []; m = L.M
    rng = np.random.default_rng([20261022, key8, ai, sd, m, sidx])
    for d_ in range(L.NDRAW):
        sel_ = sampler(rng, fl, m); F = [fl[i] for i in sel_]; clicks = sum(len(fl_posv["c0"][f]) for f in F); ni = np.array([len(fl_neg["c0"][f]) for f in F])
        for eps, dl in L.SETTINGS:
            r_ = np.random.default_rng(int(rng.integers(1 << 30))); [r_.choice(n_, 2, replace=False) if n_ >= 2 else np.arange(n_) for n_ in ni]
        D = {}
        for k in L.CEN:
            negs = [fl_neg[k][f] for f in F]; pooled = np.concatenate(negs); grid = np.unique(np.quantile(pooled, 1 - np.arange(L.G) / L.G))[::-1]
            X = np.stack([ni[i] - np.searchsorted(negs[i], grid, side="left") for i in range(m)]).astype(float)
            posl = np.concatenate([fl_posv[k][f] for f in F]) if clicks > 0 else np.array([]); D[k] = (negs, pooled, grid, X, posl)
        cache = {}

        def ubf(rule, k, j):
            if (rule, k, j) not in cache: cache[(rule, k, j)] = L.b1u(D[k][3][:, j], ni, 0.10 if rule == "U1" else 0.10 / 3)
            return cache[(rule, k, j)]

        def rec(name, eps, k, t, cert):
            sn = sens_at(k, t)
            rows.append(dict(site=site, seed=sd, draw=d_, eps=eps, rule=name, cand=k, cert=bool(cert), fpr=float(fpr_at(k, t)), sens=float(sn),
                             dse=float(sn - pop["orac_unc"][eps]), dso=float(sn - pop["orac_own"][(k, eps)]), ni=ni.tolist(), fields=[str(f) for f in F]))
        for eps, dl in L.SETTINGS:
            T3 = {k: L.fixed_sequence(D[k][2], lambda j, k=k: ubf("U3", k, j), eps) for k in L.CEN}
            cert_c = [k for k in L.CEN if T3[k][1]]
            if not cert_c: rec("S3", eps, "c0", T3["c0"][0], False)
            else:
                if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
                else:
                    sl = [(D[k][4] >= T3[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
                rec("S3", eps, ch, T3[ch][0], True)
            for k in ("c0", "c1s"):
                t_, c_, _ = L.fixed_sequence(D[k][2], lambda j, k=k: ubf("U1", k, j), eps); rec("B1u_" + k, eps, k, t_, c_)
            for name, (k, fn) in (extra or {}).items():
                t_, c_, _ = L.fixed_sequence(D[k][2], lambda j, k=k, fn=fn: fn(D[k][3][:, j], ni, (site, sd, d_, eps, j)), eps); rec(name, eps, k, t_, c_)
    return rows


def run_task(args):
    """One realisation of one cell: all rules on its draws.  Returns compact rows (without ni and fields), the bound exceptions and the run time."""
    design, ci, sd, ndraw = args; m = mods(); L = m["L"]; L.M = M; L.NDRAW = ndraw; t0 = time.time()
    site, smp = site_for(design, ci, sd); exc = defaultdict(Counter); extra = make_extra(design, ci, exc); inp = Synth({sd: site})
    rows = L.replay(inp, sd, SITE, extra=extra)[0] if smp is None else replay_hook(L, inp, sd, SITE, extra, smp)
    rows = [dict({k: v for k, v in r.items() if k not in ("ni", "fields", "site")}, design=design, cell=ci) for r in rows]
    return dict(design=design, cell=ci, seed=sd, rows=rows, exc={f"{e}|{n}": dict(c) for (e, n), c in exc.items()}, seconds=time.time() - t0)


def anchor_b():
    m = mods(); E10, E12 = m["E10"], m["E12"]
    for kw in (dict(tau=0.5, r=100, pi=0.05, tail="normal", F=300), dict(tau=0.0, r=20, pi=0.2, tail="t5", F=60)):
        x1 = E10.make_site(np.random.default_rng([1, 2, 3]), kw["tau"], kw["r"], kw["pi"], kw["tail"], F=kw["F"]); x2 = E12.make_site2(E10, np.random.default_rng([1, 2, 3]), **kw)[:3]
        if not all(np.array_equal(u, v) for u, v in zip(x1, x2)): raise AnchorError(("b: make_site2 with neutral parameters differs from make_site", kw))


def anchor_c():
    m = mods(); L, E12 = m["L"], m["E12"]; L.M = M; keep = L.NDRAW; L.NDRAW = 5
    try:
        site, _ = site_for("e10", 0, 41)
        lib = L.replay(Synth({41: site}), 41, SITE, extra=make_extra("e10", 0, defaultdict(Counter)))[0]
        mine = replay_hook(L, Synth({41: site}), 41, SITE, make_extra("e10", 0, defaultdict(Counter)), E12.uniform_sampler)
    finally: L.NDRAW = keep
    kk = lambda r: (r["draw"], r["eps"], r["rule"]); A = {kk(r): r for r in lib}; B = {kk(r): r for r in mine}
    if set(A) != set(B) or len(A) != len(lib) or not all(all(A[k][f] == B[k][f] for f in CORE) for k in A): raise AnchorError(("c: replay_hook differs from m3_replay_lib.replay", len(A), len(B)))
    return len(A)


def anchor_d(results, e12_rows_path):
    with open(e12_rows_path) as fh: stored = json.load(fh)
    ref = {(r["cell"], r["seed"], r["draw"], round(r["eps"], 6), r["rule"]): r for r in stored if r["rule"] in ANCHOR_RULES}
    n = 0
    for t in results:
        if t["design"] != "e12": continue
        for r in t["rows"]:
            if r["rule"] not in ANCHOR_RULES: continue
            s = ref.get((r["cell"], r["seed"], r["draw"], round(r["eps"], 6), r["rule"]))
            if s is None or not all(s[f] == r[f] for f in CORE): raise AnchorError(("d: E12 row differs from the stored row", r["cell"], r["seed"], r["draw"], r["eps"], r["rule"]))
            n += 1
    return n


def anchor_e(summary, e10_json_path, reps, ndraw):
    if reps != 5 or ndraw != 200: return "not_applicable (non-default size)"
    with open(e10_json_path) as fh: e10 = json.load(fh)
    if e10["reps"] != 5 or e10["ndraw"] != 200: raise AnchorError("e: E10 size")
    n = 0
    for key, v in summary.items():
        design, ci, rule, eps = key.split("|")
        if design != "e10" or rule not in ("S3", "B1u_c1s"): continue
        ref = e10["cells"][ci]["library_rules"][f"{rule}|{float(eps)}"]
        mine = v["coverage"]["per_realisation"]; refc = ref["coverage_per_realisation"]
        ok = len(mine) == len(refc) and all(abs(a - b) <= 1e-12 for a, b in zip(mine, refc))
        ok = ok and abs(v["cannot_certify"]["mean"] - ref["cannot_certify"]) <= 1e-12 and abs(v["ds_own"]["mean"] - ref["ds_own"]) <= 1e-12
        if not ok: raise AnchorError(("e: E10 values differ from e10.json", ci, rule, eps))
        n += 1
    return f"passed ({n} cell x rule x eps rows of S3 and B1u_c1s equal e10.json within 1e-12)"


def summarise(results, EB):
    by = defaultdict(lambda: defaultdict(list)); exc = defaultdict(Counter)
    for t in results:
        for r in t["rows"]: by[(t["design"], t["cell"], r["rule"], r["eps"])][t["seed"]].append(r)
        for k, c in t["exc"].items():
            e, name = k.split("|"); exc[(t["design"], t["cell"], name, float(e))].update(c)
    out = {}
    for (design, ci, rule, eps), per in sorted(by.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], -kv[0][3])):
        sds = sorted(per); vals = {mname: [] for mname in METRICS}
        for sd in sds:
            R_ = per[sd]; fpr = np.array([r["fpr"] for r in R_])
            vals["coverage"].append(float((fpr <= eps).mean())); vals["cannot_certify"].append(float(1 - np.mean([r["cert"] for r in R_])))
            vals["sens"].append(float(np.mean([r["sens"] for r in R_]))); vals["ds_eps"].append(float(np.mean([r["dse"] for r in R_])))
            vals["ds_own"].append(float(np.mean([r["dso"] for r in R_]))); vals["fpr_mean"].append(float(fpr.mean()))
        n = int(sum(len(per[sd]) for sd in sds)); cell = {"n": n, "n_realisations": len(sds)}
        for mname, v in vals.items():
            cell[mname] = dict(mean=float(np.mean(v)), per_realisation=v, mcse=(float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None))
        c = cell["coverage"]["mean"]; cell["coverage"]["mcse_binomial"] = float(np.sqrt(c * (1 - c) / n))
        cell["coverage"]["ci90_seedcluster"] = EB.seed_cluster_ci(vals["coverage"]) if len(sds) == 5 else None
        ex = exc.get((design, ci, rule, eps), Counter()); nex = int(sum(ex.values()))
        cell.update(n_bound_exceptions=nex, exception_types=[[k, v] for k, v in ex.most_common(3)], flag=("exceptions>1%" if nex > 0.01 * n else None))
        out[f"{design}|{ci}|{rule}|{eps}"] = cell
    return out


def parse_cells(s):
    if not s: return [(d, ci) for d in ("e10", "e12") for ci in range(N_CELLS[d])]
    out = []
    for part in s.split(";"):
        d, idx = part.split(":"); assert d in N_CELLS, ("unknown design", d)
        for x in idx.split(","):
            ci = int(x); assert 0 <= ci < N_CELLS[d], ("cell out of range", d, ci); out.append((d, ci))
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser(); ap.add_argument("--e10_json", required=True); ap.add_argument("--e12_rows", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=16); ap.add_argument("--cells", default=None); ap.add_argument("--reps", type=int, default=5); ap.add_argument("--ndraw", type=int, default=200)
    a = ap.parse_args(argv); out = Path(a.out)
    m = mods(); m["E12"].gate(argv)
    if out.exists() or not OUT_RE.match(out.name) or not out.parent.is_dir(): raise SystemExit("--out must be a new directory simw3_<UTC yyyymmddThhmmssZ> in an existing parent")
    if not (1 <= a.workers <= MAX_WORKERS) or not (1 <= a.reps <= 5) or not (1 <= a.ndraw <= 200): raise SystemExit("workers 1..32, reps 1..5, ndraw 1..200")
    # all anchors before anything is written
    check_pins()
    for k, p in (("e10_json", a.e10_json), ("e12_rows", a.e12_rows)):
        if sha(p) != FILE_PINS[k]: raise AnchorError(("input file differs from its pin", k))
    anchor_b(); n_c = anchor_c()
    tasks = [(d, ci, 41 + rid, a.ndraw) for d, ci in parse_cells(a.cells) for rid in range(a.reps)]; t0 = time.time()
    with get_context("spawn").Pool(a.workers) as pool:
        results = []
        for i, res in enumerate(pool.imap_unordered(run_task, tasks), 1):
            results.append(res); print(f"task {i}/{len(tasks)} {res['design']} cell {res['cell']} seed {res['seed']} {res['seconds']:.0f} s", flush=True)
    n_d = anchor_d(results, a.e12_rows); summary = summarise(results, m["EB"]); e_txt = anchor_e(summary, a.e10_json, a.reps, a.ndraw)
    cells = {f"{d}|{ci}": cell_spec(d, ci) for d, ci in parse_cells(a.cells)}
    rep = dict(version="simw3_r1", request="实验管理 2026-11-24: Section 6 cluster-aware baselines in the Section 5 simulation", script_sha256=sha(__file__), pins=PINS, file_pins=FILE_PINS,
               settings=dict(reps=a.reps, ndraw=a.ndraw, m=M, eps=list(EPSS), delta=m["W3"].DELTA, candidate=CAND, cb_B=m["W3"].B,
                             cb_stream="default_rng([20261031, 0, seed, draw, round(100 eps), candidate code 1, 13, design code, cell])", workers=a.workers, blas_threads=1),
               anchors=dict(b="passed (make_site2 neutral == make_site bit for bit)", c=f"passed (replay_hook == m3_replay_lib.replay, {n_c} rows, every rule incl. baselines)",
                            d=f"passed ({n_d} S3 / B1u_c0 / B1u_c1s rows equal the stored E12 rows)", e=e_txt),
               mcse_definition="SD (ddof 1) of the per-realisation values / sqrt(number of realisations); coverage also: binomial SE sqrt(c(1 - c) / n) over all draws (ignores clustering by realisation) and the seed-cluster 90% interval (" + m["EB"].CI_DEFINITION + ")",
               cells_spec=cells, results=summary, task_seconds={f"{t['design']}|{t['cell']}|{t['seed']}": t["seconds"] for t in results}, seconds=time.time() - t0,
               note="synthetic; no project data; comparison baselines with the W3 definitions; no rule, constant or baseline definition changed; every cell is reported")
    rows = [r for t in sorted(results, key=lambda t: (t["design"], t["cell"], t["seed"])) for r in t["rows"]]
    files = {"simw3_summary.json": json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False), "simw3_rows.json": json.dumps(rows, allow_nan=False)}
    out.mkdir(parents=False, exist_ok=False)
    for n, txt in files.items():
        with open(out / n, "x") as fh: fh.write(txt)
    with open(out / "SHA256SUMS", "x") as fh: fh.write("".join(f"{sha(out / n)}  {n}\n" for n in files))
    print(out); return 0


if __name__ == "__main__":
    try: sys.exit(main())
    except AnchorError as e: print("ANCHOR MISMATCH (exit 3, nothing written):", e, file=sys.stderr); sys.exit(3)

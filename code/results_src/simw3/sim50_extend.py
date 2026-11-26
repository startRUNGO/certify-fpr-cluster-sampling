# P5a SIM50: the Section 5 simulation extended from 5 to 50 synthetic sites per cell (user-approved, 实验管理 2026-11-24) for S3, SPLIT6 and the three Section 6 baselines CB_c1s, GEE_c1s, BB_c1s.
# WRITTEN, NOT RUN, until 实验管理 reviews it.  Synthetic only: no project data, no image, no classifier.  No rule, constant or baseline definition is changed; every cell is reported.
# Realisations: the existing 5 (seeds 41..45) are NOT rerun for new data; 45 new realisations, seeds 46..90, with exactly the E10 / E12 site generation (E10: make_site, default_rng([20261103, 10, cell, seed, 0]);
#   E12: make_site2 with e12_stress.stream(cell, seed, 0), two_stage_sampler for C2), the same draw stream, m = 22, eps in {0.10, 0.05}, 200 draws per realisation.
# Rules and implementations (reused, never copied here):
#   pass 1  S3 and SPLIT6: e12_stress.run_site_seed_e12 (aef35708; the E12 copy of m3_split_dev_eval.run_site_seed with the field-sampler hook, anchored in E12 to the library), with the Platt fit of the
#           realisation's synthetic source set exactly as E10 / E12 (e10_synth.source_set with stream purpose 1, e10_synth.platt_synth with m3_k1v_precheck).  The function also computes S3c and B1u_c0 / B1u_c1s
#           on the same draws; they cost well under 1% of a task and are discarded (removing them would need a copy of the function).
#   pass 2  CB_c1s, GEE_c1s, BB_c1s: simw3_baselines.make_extra through m3_replay_lib.replay (uniform sampler) or simw3_baselines.replay_hook (C2), exactly as the simw3 run (45f1ebcd).  Its S3 rows are kept
#           only for anchor (c).
# The original 5 realisations in the combined result: E12 S3 / SPLIT6 from the stored rows (e12_rows.json); the three baselines of every cell from the simw3 full run rows (simw3_rows.json); E10 S3 / SPLIT6
#   have no stored rows, so pass 1 is REPLAYED for seeds 41..45 of the E10 cells (cheap, no new data) and must equal e10.json (anchor b) before it is used.
# Anchors (exit 3, nothing written): (p) pins of every module and of e10.json, e12_rows.json, simw3_rows.json (the simw3 pins are empty until its full run is done: the script refuses to run);
#   (a) seed 41 of ANCHOR_CELLS, run first: pass-1 S3 / SPLIT6 rows of E12 cells equal e12_rows.json and pass-2 CB / GEE / BB rows of every anchor cell equal simw3_rows.json (cert, cand, fpr, sens, dse, dso exactly);
#   (b) E10 replay: per-realisation coverage of S3 and SPLIT6 equal e10.json and cannot_certify / ds_own equal e10.json within 1e-12 (all 5 seeds, 200 draws);
#   (c) every task: pass-2 S3 rows equal pass-1 S3 rows.
# Output per design, cell, rule, eps and block (orig5 = seeds 41..45, new45 = seeds 46..90, all50): coverage, cannot_certify, sens, ds_eps, ds_own, fpr_mean, each with per-realisation values and MCSE
#   (SD ddof 1 / sqrt(number of realisations)); coverage also with the binomial SE sqrt(c(1 - c) / n) and the seed-cluster 90% interval (ebatch_lib.seed_cluster_ci); SPLIT6 / baseline bound exceptions.
# Cost: workers x 1 BLAS thread (default 24, at most 32); tasks are (cell, seed); the BB profile likelihood dominates (simw3 full run: about 700 s per 200-draw task); pass 1 adds about 20 s per task.
#   All 34 cells: 1530 new tasks + 50 E10 replay tasks (pass 1 only) + 5 anchor tasks, about 310 CPU hours, about 13 h wall at 24 workers.  --cells restricts the cells (option A) or splits them into batches (option B);
#   each run is self-contained (its cells get their full 50-realisation result).  Memory: about 0.3 GB per worker plus about 3 GB in the parent; disk about 150 MB per full run (compressed rows).
# Output: one NEW directory fifty_<UTC yyyymmddThhmmssZ>: sim50_summary.json, sim50_rows.npz (new-seed and E10-replay rows, columnar), SHA256SUMS; 'x' mode, allow_nan=False.
# v1.2 (实验管理 2026-11-24): output directory fifty_<UTC> (was sim50_<UTC>): e12_stress.gate forbids the substring 'm5' (any case) in every argument and 'sim50' contains it;
#   the gate is unchanged; package and output parents are named p5fifty_pkg_<UTC> and <DATA>/p5_fifty/.  FORBID_MOD (^(m5_|e8_|r8_|e7_) on loaded module names) is not hit by sim50_extend or simw3_baselines.
# Usage: python3 -I -B sim50_extend.py --e10_json <ebatch_e10.json> --e12_rows <e12_rows.json> --simw3_rows <simw3_rows.json> --out <parent>/fifty_<UTC> [--workers 24] [--cells "e10:0,3;e12:11-20"]
#        [--new_seeds 46-90] [--ndraw 200]
import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[_k] = "1"
import argparse, hashlib, io, json, re, sys, time
from collections import Counter, defaultdict
from multiprocessing import get_context
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; R = HERE.parent
sys.path[:0] = [str(HERE), str(R / "e10"), str(R / "e12"), str(R / "e_batch"), str(R / "m3" / "split"), str(R / "m3" / "k2prep"), str(R / "m3" / "k1vb"), str(R / "m3" / "stress"), str(R / "c5b")]
PINS = {"simw3_baselines": "45f1ebcdc1cae62943b4ae1a3d23ee31d8ecc1b2c25c302ec0447fe8291001a3", "m3_replay_lib": "ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de",
        "m3_w3_baselines": "249be3977f9c7e06dbcaea0830baef92a0fb65b7430d7c802a457cb0c809c20e", "e10_synth": "00be97b1f15a0fd5ac7c03b06d41f6adc0ac44a2c11b165e54748fbe68f7e08e",
        "e12_stress": "aef3570877fdf7427fb469d4e99d3341547342bb36fb86366e47c69b971ccd59", "ebatch_lib": "359b6bce2dfdbf50a2acb8b909d631ab40d445447275915e2318e70981ef007f",
        "m3_split_dev_eval": "cbf1e7f4f33dd1c74cbfd0f3ac4d34a397e8a7ef44b82c9365ca117733bf6457", "m3_k12_dev_eval": "0091a3337074a227e404f0585966b87a9f5e9e55cfbb8dd766df7e4618893b67",
        "m3_stress": "75eb9afd3a130b00170757c4027f8b158bf95deb0ad6124229e9292c60e2d1ad", "m3_k1vb_dev_eval": "7b7e6ecda22b789c573451c757b6e197d747f77e10bcfc201ccdb85921e25abd",
        "m3_k1v_precheck": "0e6aa9f8f423dbaed286db969afa8a14c9ce477e37bba962732d3b411002147c", "s8_bounds": "54a2af8d3704c16b0994dcdc584d4f407cd9564cc07d49ad74cdf7c4cad5ae7f"}
FILE_PINS = {"e10_json": "a25bc85c40a9b85c742444deb16b704107deb399f49b1cb15d9c94d6a79e9bc8", "e12_rows": "ece7c768005f0a8b79e0e4a883a1a50dc9f6f70e5b70995e4ece4d9949d54a20",
             "simw3_rows": "bed0438a1f0717d98482502a0baa93d2fddfefb88e3179048e44ea8bff2d3794"}
SIMW3_RUN = "simw3_20261124T194046Z"
RULES = ("S3", "SPLIT6", "CB_c1s", "GEE_c1s", "BB_c1s"); P1 = ("S3", "SPLIT6"); P2 = ("CB_c1s", "GEE_c1s", "BB_c1s"); CANDS = ("c0", "c1s", "c2s"); DESIGNS = ("e10", "e12")
CORE = ("cert", "cand", "fpr", "sens", "dse", "dso"); METRICS = ("coverage", "cannot_certify", "sens", "ds_eps", "ds_own", "fpr_mean")
ORIG = (41, 42, 43, 44, 45); N_CELLS = {"e10": 10, "e12": 24}; EPSS = (0.1, 0.05); M = 22; SITE = "bbbc041_all"
ANCHOR_CELLS = (("e10", 0), ("e10", 3), ("e12", 5), ("e12", 14), ("e12", 19))
OUT_RE = re.compile(r"^fifty_\d{8}T\d{6}Z$"); MAX_WORKERS = 32
_MOD = {}


class AnchorError(Exception): pass


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def mods():
    if not _MOD:
        import simw3_baselines as W3S
        m = W3S.mods()
        import m3_split_dev_eval as SP, m3_k1v_precheck as KV, m3_k12_dev_eval as K12, m3_stress as ST, m3_k1vb_dev_eval as VB, s8_bounds as B8
        for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[_k] = "1"
        if SP.L is not m["L"]: raise AnchorError("m3_split_dev_eval did not import the same m3_replay_lib")
        _MOD.update(m, W3S=W3S, SP=SP, KV=KV, K12=K12, ST=ST, VB=VB, B8=B8)
    return _MOD


def check_pins():
    m = mods(); names = {"simw3_baselines": "W3S", "m3_replay_lib": "L", "m3_w3_baselines": "W3", "e10_synth": "E10", "e12_stress": "E12", "ebatch_lib": "EB", "m3_split_dev_eval": "SP",
                         "m3_k12_dev_eval": "K12", "m3_stress": "ST", "m3_k1vb_dev_eval": "VB", "m3_k1v_precheck": "KV", "s8_bounds": "B8"}
    bad = [n for n, k in names.items() if sha(m[k].__file__) != PINS[n]]
    if bad: raise AnchorError(("module differs from its pin", bad))


def site_platt(design, ci, sd):
    """Site, sampler (None = uniform) and Platt of one realisation, exactly as E10 / E12."""
    m = mods(); E10, E12, KV = m["E10"], m["E12"], m["KV"]; cell = m["W3S"].cell_spec(design, ci)
    site, smp = m["W3S"].site_for(design, ci, sd)
    src = [20261103, 10, ci, sd, 1] if design == "e10" else E12.stream(ci, sd, 1)
    ss, ys, _ = E10.source_set(np.random.default_rng(src), cell["r"], cell["pi"], cell["tail"]); platt, pinfo = E10.platt_synth(KV, ss, ys)
    return site, smp, platt, bool(pinfo["fallback"])


def to_cols(rows, design, ci):
    keep = [r for r in rows if r["rule"] in RULES]
    return dict(design=np.array([DESIGNS.index(design)] * len(keep), np.int8), cell=np.array([ci] * len(keep), np.int16), seed=np.array([r["seed"] for r in keep], np.int16),
                draw=np.array([r["draw"] for r in keep], np.int16), eps=np.array([r["eps"] for r in keep], float), rule=np.array([RULES.index(r["rule"]) for r in keep], np.int8),
                cand=np.array([CANDS.index(r["cand"]) for r in keep], np.int8), cert=np.array([r["cert"] for r in keep], bool), fpr=np.array([r["fpr"] for r in keep], float),
                sens=np.array([r["sens"] for r in keep], float), dse=np.array([r["dse"] for r in keep], float), dso=np.array([r["dso"] for r in keep], float))


def run_task(args):
    """kind: 'new' (passes 1 and 2), 'e10replay' (pass 1 only), 'anchor' (passes 1 and 2, seed 41).  Returns columnar rows, pass-1 / pass-2 S3 mismatches, exceptions, Platt fallback, time."""
    kind, design, ci, sd, ndraw = args; m = mods(); L, SP, E12, W3S = m["L"], m["SP"], m["E12"], m["W3S"]; L.M = M; L.NDRAW = ndraw; t0 = time.time()
    site, smp, platt, fb = site_platt(design, ci, sd); inp = W3S.Synth({sd: site})
    SP.E.EXC.clear(); SP.E.DEG.clear(); SP.E.DROP.clear(); SP.CHOSEN.clear(); SP.C2STOP.clear()
    r1 = [r for r in E12.run_site_seed_e12(SP, inp, sd, SITE, platt, smp or E12.uniform_sampler) if r["rule"] in P1]
    exc = {f"SPLIT6|{k}": dict(v) for k, v in SP.E.EXC.items()}
    rows, s3_mismatch = list(r1), 0
    if kind != "e10replay":
        ex2 = defaultdict(Counter); extra = W3S.make_extra(design, ci, ex2)
        r2 = L.replay(inp, sd, SITE, extra=extra)[0] if smp is None else W3S.replay_hook(L, inp, sd, SITE, extra, smp)
        a = {(r["draw"], r["eps"]): r for r in r1 if r["rule"] == "S3"}; b = {(r["draw"], r["eps"]): r for r in r2 if r["rule"] == "S3"}
        s3_mismatch = int(set(a) != set(b)) + sum(1 for k in a if k in b and not all(a[k][f] == b[k][f] for f in CORE))
        rows += [r for r in r2 if r["rule"] in P2]; exc.update({f"{e}|{n}": dict(c) for (e, n), c in ex2.items()})
    for r in rows: r["seed"] = sd
    return dict(kind=kind, design=design, cell=ci, seed=sd, cols=to_cols(rows, design, ci), rows_core=([{k: r[k] for k in ("draw", "eps", "rule") + CORE} for r in rows] if kind == "anchor" else None),
                s3_mismatch=s3_mismatch, exc=exc, platt_fallback=fb, seconds=time.time() - t0)


def parse_cells(s):
    if not s: return [(d, ci) for d in DESIGNS for ci in range(N_CELLS[d])]
    out = []
    for part in s.split(";"):
        d, idx = part.split(":"); assert d in N_CELLS, ("unknown design", d)
        for x in idx.split(","):
            lo, hi = (x.split("-") + [x])[:2] if "-" in x else (x, x)
            for ci in range(int(lo), int(hi) + 1):
                assert 0 <= ci < N_CELLS[d], ("cell out of range", d, ci)
                if (d, ci) not in out: out.append((d, ci))
    return out


def parse_seeds(s):
    lo, hi = (s.split("-") + [s])[:2] if "-" in s else (s, s); seeds = list(range(int(lo), int(hi) + 1))
    if not seeds or min(seeds) < 46 or max(seeds) > 90: raise SystemExit("--new_seeds must lie in 46..90")
    return seeds


def anchor_a(res, e12_rows, w3_rows):
    n = 0
    for t in res:
        for r in t["rows_core"]:
            if r["rule"] in P1 and t["design"] == "e12": ref = e12_rows.get((t["cell"], 41, r["draw"], round(r["eps"], 6), r["rule"]))
            elif r["rule"] in P2: ref = w3_rows.get((t["design"], t["cell"], 41, r["draw"], round(r["eps"], 6), r["rule"]))
            else: continue
            if ref is None or not all(ref[f] == r[f] for f in CORE): raise AnchorError(("a: seed-41 row differs from the stored row", t["design"], t["cell"], r["draw"], r["eps"], r["rule"]))
            n += 1
    return n


def anchor_b(e10rep, e10_json):
    e10 = json.load(open(e10_json)); n = 0
    by = defaultdict(lambda: defaultdict(list))
    for t in e10rep:
        c = t["cols"]
        for i in range(len(c["seed"])): by[(t["cell"], RULES[c["rule"][i]], float(c["eps"][i]))][t["seed"]].append(i)
        t["_c"] = c
    cols = {t["cell"]: {} for t in e10rep}
    for t in e10rep: cols[t["cell"]][t["seed"]] = t["cols"]
    for (ci, rule, eps), per in by.items():
        ref = e10["cells"][str(ci)]["library_rules"][f"{rule}|{eps}"]
        cov = [float(np.mean(cols[ci][sd]["fpr"][per[sd]] <= eps)) for sd in ORIG]; cc = [float(1 - np.mean(cols[ci][sd]["cert"][per[sd]])) for sd in ORIG]
        dso = [float(np.mean(cols[ci][sd]["dso"][per[sd]])) for sd in ORIG]
        if any(abs(u - v) > 1e-12 for u, v in zip(cov, ref["coverage_per_realisation"])) or abs(float(np.mean(cc)) - ref["cannot_certify"]) > 1e-12 or abs(float(np.mean(dso)) - ref["ds_own"]) > 1e-12:
            raise AnchorError(("b: E10 replay differs from e10.json", ci, rule, eps))
        n += 1
    return n


def per_real(fpr, cert, sens, dse, dso, eps):
    return dict(coverage=float(np.mean(np.asarray(fpr) <= eps)), cannot_certify=float(1 - np.mean(cert)), sens=float(np.mean(sens)), ds_eps=float(np.mean(dse)), ds_own=float(np.mean(dso)), fpr_mean=float(np.mean(fpr)))


def block(per, EB):
    """per: {seed: per-realisation dict}; returns metrics with per-realisation values, MCSE, binomial SE and seed-cluster interval for coverage."""
    sds = sorted(per); out = {"n_realisations": len(sds), "seeds": [int(s) for s in sds]}
    for mn in METRICS:
        v = [per[s][mn] for s in sds]; out[mn] = dict(mean=float(np.mean(v)), per_realisation=v, mcse=(float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None))
    n = int(sum(per[s]["_n"] for s in sds)); c = out["coverage"]["mean"]; out["n"] = n
    out["coverage"]["mcse_binomial"] = float(np.sqrt(c * (1 - c) / n)); out["coverage"]["ci90_seedcluster"] = EB.seed_cluster_ci(out["coverage"]["per_realisation"])
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser()
    for k in ("--e10_json", "--e12_rows", "--simw3_rows", "--out"): ap.add_argument(k, required=True)
    ap.add_argument("--workers", type=int, default=24); ap.add_argument("--cells", default=None); ap.add_argument("--new_seeds", default="46-90"); ap.add_argument("--ndraw", type=int, default=200)
    a = ap.parse_args(argv); out = Path(a.out); m = mods(); m["E12"].gate(argv); EB = m["EB"]
    if not FILE_PINS["simw3_rows"]: raise SystemExit("simw3_rows pin is empty: fill it through a reviewed diff after the simw3 full run")
    if out.exists() or not OUT_RE.match(out.name) or not out.parent.is_dir(): raise SystemExit("--out must be a new directory fifty_<UTC yyyymmddThhmmssZ> in an existing parent")
    if not (1 <= a.workers <= MAX_WORKERS) or not (1 <= a.ndraw <= 200): raise SystemExit("workers 1..32, ndraw 1..200")
    cells = parse_cells(a.cells); seeds = parse_seeds(a.new_seeds); t0 = time.time()
    check_pins()
    for k, p in (("e10_json", a.e10_json), ("e12_rows", a.e12_rows), ("simw3_rows", a.simw3_rows)):
        if sha(p) != FILE_PINS[k]: raise AnchorError(("input file differs from its pin", k))
    e12_rows = {(r["cell"], r["seed"], r["draw"], round(r["eps"], 6), r["rule"]): r for r in json.load(open(a.e12_rows)) if r["rule"] in P1}
    w3_rows = {(r["design"], r["cell"], r["seed"], r["draw"], round(r["eps"], 6), r["rule"]): r for r in json.load(open(a.simw3_rows)) if r["rule"] in P2}
    with get_context("spawn").Pool(a.workers) as pool:
        # (a) seed-41 anchors first, full 200 draws
        res_a = pool.map(run_task, [("anchor", d, ci, 41, 200) for d, ci in ANCHOR_CELLS]); n_a = anchor_a(res_a, e12_rows, w3_rows)
        tasks = [("e10replay", d, ci, sd, 200) for d, ci in cells if d == "e10" for sd in ORIG] + [("new", d, ci, sd, a.ndraw) for d, ci in cells for sd in seeds]
        res = []
        for i, r in enumerate(pool.imap_unordered(run_task, tasks), 1):
            res.append(r); print(f"task {i}/{len(tasks)} {r['kind']} {r['design']} cell {r['cell']} seed {r['seed']} {r['seconds']:.0f} s", flush=True)
    bad_c = [(t["design"], t["cell"], t["seed"]) for t in res_a + res if t["s3_mismatch"]]
    if bad_c: raise AnchorError(("c: pass-2 S3 rows differ from pass-1 S3 rows", bad_c[:5]))
    e10rep = [t for t in res if t["kind"] == "e10replay"]; n_b = anchor_b(e10rep, a.e10_json) if e10rep else 0
    # per-realisation values: new seeds from this run; orig5 from the E10 replay, the stored E12 rows and the simw3 rows
    per = defaultdict(dict)
    for t in res:
        c = t["cols"]
        for ri, rule in enumerate(RULES):
            for eps in EPSS:
                sel = (c["rule"] == ri) & np.isclose(c["eps"], eps)
                if sel.any(): per[(t["design"], t["cell"], rule, eps)][t["seed"]] = dict(per_real(c["fpr"][sel], c["cert"][sel], c["sens"][sel], c["dse"][sel], c["dso"][sel], eps), _n=int(sel.sum()))
    for d, ci in cells:
        for eps in EPSS:
            for sd in ORIG:
                for rule in RULES:
                    if rule in P1 and d == "e10": continue
                    src = ([e12_rows[(ci, sd, dr, round(eps, 6), rule)] for dr in range(200)] if rule in P1 else [w3_rows[(d, ci, sd, dr, round(eps, 6), rule)] for dr in range(200)])
                    per[(d, ci, rule, eps)][sd] = dict(per_real([r["fpr"] for r in src], [r["cert"] for r in src], [r["sens"] for r in src], [r["dse"] for r in src], [r["dso"] for r in src], eps), _n=len(src))
    summary = {}
    for (d, ci, rule, eps), p in sorted(per.items(), key=lambda kv: (kv[0][0], kv[0][1], RULES.index(kv[0][2]), -kv[0][3])):
        o5 = {s: v for s, v in p.items() if s in ORIG}; n45 = {s: v for s, v in p.items() if s not in ORIG}
        summary[f"{d}|{ci}|{rule}|{eps}"] = dict(orig5=block(o5, EB), new45=block(n45, EB), all50=block(p, EB))
    exc = defaultdict(Counter)
    for t in res:
        for k, v in t["exc"].items(): exc[f"{t['design']}|{t['cell']}|{k}"].update(v)
    cols = {k: np.concatenate([t["cols"][k] for t in res]) for k in res[0]["cols"]}
    rep = dict(version="sim50_r1", request="实验管理 2026-11-24 (user-approved): 50 synthetic sites per cell for S3, SPLIT6, CB_c1s, GEE_c1s, BB_c1s", script_sha256=sha(__file__), pins=PINS, file_pins=FILE_PINS,
               simw3_run=SIMW3_RUN, settings=dict(new_seeds=seeds, orig_seeds=list(ORIG), ndraw=a.ndraw, m=M, eps=list(EPSS), workers=a.workers, blas_threads=1, cells=[f"{d}:{ci}" for d, ci in cells]),
               anchors=dict(a=f"passed ({n_a} seed-41 rows of {len(ANCHOR_CELLS)} anchor cells equal e12_rows.json / simw3_rows.json)", b=(f"passed ({n_b} E10 cell x rule x eps rows equal e10.json)" if e10rep else "not_applicable (no E10 cell)"),
                            c="passed (pass-2 S3 rows equal pass-1 S3 rows in every task)"),
               mcse_definition="SD (ddof 1) of the per-realisation values / sqrt(number of realisations); coverage also: binomial SE sqrt(c(1 - c) / n) over all draws and the seed-cluster 90% interval (" + EB.CI_DEFINITION + ")",
               sources_orig5=dict(e12_S3_SPLIT6="e12_rows.json", e10_S3_SPLIT6="replay in this run (anchor b)", baselines=f"{SIMW3_RUN}/simw3_rows.json"),
               results=summary, exceptions={k: dict(v) for k, v in exc.items()}, platt_fallback=sorted({(t["design"], t["cell"], t["seed"]) for t in res if t["platt_fallback"]}),
               task_seconds={f"{t['kind']}|{t['design']}|{t['cell']}|{t['seed']}": t["seconds"] for t in res}, seconds=time.time() - t0,
               note="synthetic; no project data; rules and baselines reused unchanged; every cell is reported")
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False); buf = io.BytesIO(); np.savez_compressed(buf, rule_names=np.array(RULES), cand_names=np.array(CANDS), design_names=np.array(DESIGNS), **cols)
    out.mkdir(parents=False, exist_ok=False)
    with open(out / "sim50_summary.json", "x") as fh: fh.write(txt)
    with open(out / "sim50_rows.npz", "xb") as fh: fh.write(buf.getvalue())
    with open(out / "SHA256SUMS", "x") as fh: fh.write("".join(f"{sha(out / n)}  {n}\n" for n in ("sim50_summary.json", "sim50_rows.npz")))
    print(out); return 0


if __name__ == "__main__":
    try: sys.exit(main())
    except AnchorError as e: print("ANCHOR MISMATCH (exit 3, nothing written):", e, file=sys.stderr); sys.exit(3)

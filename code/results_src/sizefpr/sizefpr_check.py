# P5a size-FPR check (实验管理 2026-11-25, priority before submission).  WRITTEN, NOT RUN, until 实验管理 reviews it.  Synthetic only: no project data, no image, no classifier.
# WHY: in P5b G0/G1 (98d2e898) the B1u design-effect bound lost its level when field size correlates with field FPR.  P5a states that DE-Grid reached 1 - delta in 48 of 48 stress settings,
#   and its stress grid has no size-FPR correlation cell.  This run adds four such cells and evaluates P5a's FROZEN rules on them through the SAME code path as sim50.
# GENERATOR: sizecorr_gen.make_site_sizecorr (new module, pinned below; never edited into e12_stress): E12 cell-B sizes (log-normal dispersion size_sigma at mean r) and a field shift
#   correlated with the size latent, mu_f = tau (rho_c z_f + sqrt(1 - rho_c^2) w_f); other factors at the E12 base cell (tau 0.5, r 100, pi 0.05, normal, F 300, SEP 2.5).
#   Cells (design 'sc'): sc0 = (size_sigma 1.0, rho_c 0.4), sc1 = (1.0, 0.8), sc2 = (0.5, 0.8), sc3 = (1.0, -0.8).  Positive rho_c: large fields have higher field FPR.
# RULES (reused, never copied): sim50_extend.run_task (pinned) with the 'new' kind, i.e. pass 1 S3 (DE-Grid) and SPLIT6 (Split-Select) through e12_stress.run_site_seed_e12 with the
#   E10 / E12 Platt fit of the realisation's synthetic source set, pass 2 CB_c1s, GEE_c1s, BB_c1s through m3_replay_lib.replay with simw3_baselines.make_extra; m = 22, eps in {0.10, 0.05},
#   delta 0.10, the frozen grid and K.  The new design enters the reused code only through three registrations made in every process before any task: simw3_baselines.site_for /
#   cell_spec dispatch design 'sc' to the new generator (every other design to the original functions), simw3_baselines.DESIGN_CODE['sc'] = 13 (the CB bootstrap stream), and
#   sim50_extend.DESIGNS gains 'sc' (column coding).  Site stream: e12_stress.stream(100 + j, seed, 0) for cell scj (purpose 1 for the Platt source set, as E12), so no E12 stream is reused.
# REALISATIONS: seeds 41..90 (50 per cell, as sim50 all50), 200 draws each.
# ANCHORS (exit 3, nothing written): (p) pins of sizecorr_gen, sim50_extend and every module sim50_extend pins (its check_pins) and of e12_rows.json / simw3_rows.json;
#   (g) generator: with rho_c = 0, make_site_sizecorr equals e12_stress.make_site2 bit for bit for E12 cell 8 (size_sigma 1.0) and cell 7 (0.5) on the E12 streams of seeds 41..45;
#   (r) rules: E12 cell 8, seed 41, all 200 draws, run through run_task with the 'anchor' kind while site_for dispatches ('e12', 8) to make_site_sizecorr(rho_c = 0): every S3 / SPLIT6
#       row equals e12_rows.json and every CB / GEE / BB row equals simw3_rows.json (sim50_extend.anchor_a, fields cert, cand, fpr, sens, dse, dso);
#   (c) every task: pass-2 S3 rows equal pass-1 S3 rows (as sim50).
# REPORTED (pre-written): per cell, rule (S3, SPLIT6, CB_c1s, GEE_c1s, BB_c1s) and eps: coverage with MCSE (SD of the 50 per-realisation coverages / sqrt 50), cannot_certify, ds_eps, ds_own,
#   fpr_mean, each with per-realisation values; the realised site-level correlation between field size and field FPR at the oracle threshold, per realisation (descriptive).
# HOW IT ENTERS THE PAPER (pre-written, 实验管理):
#   if DE-Grid (S3) stays >= 1 - delta in all 4 cells x 2 targets: the cells join the stress grid of SI Table S3.1 and the 48-of-48 counts are updated to the new total as new keys;
#   if S3 falls below 1 - delta in any of them: it is reported as a measured limit in main Section 6 (limits list) and Section 4.6, and the abstract and the main claim are qualified to the
#   settings where coverage holds.  No change to the rules (a method fix is P5b / D1 material).  Split-Select and the baselines are reported alongside in either case.
# COST: per (cell, seed) task about 720 s at 200 draws (BB profile likelihood dominates, as sim50); 4 x 50 = 200 tasks + 1 anchor task = about 40 CPU hours, about 1.7 h wall at
#   24 workers on the VM (BLAS 1 thread, nice 19, MemAvailable > 20 GB, about 0.3 GB per worker).  Local is over 1 h, so the full run is a VM job (needs 实验管理 approval).
#   --smoke: local, seed 41 only, 4 draws, the generator anchor (g) only (the stored rows are on the VM), output marked smoke; never used for the paper.
# OUTPUT: one NEW directory sizefpr_<UTC>: sizefpr_summary.json, sizefpr_rows.npz, SHA256SUMS ('x' mode, allow_nan=False).
# Usage: python3 -I -B sizefpr_check.py --e12_rows <e12_rows.json> --simw3_rows <simw3_rows.json> --out <parent>/sizefpr_<UTC> [--workers 24]
#        python3 -I -B sizefpr_check.py --smoke --out <parent>/sizefpr_<UTC> [--workers 2]
import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[_k] = "1"
import argparse, hashlib, importlib.util, io, json, re, sys, time
from collections import defaultdict
from multiprocessing import get_context
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; R5 = HERE.parent
sys.path[:0] = [str(R5 / "simw3")]
PINS = {"sizecorr_gen": (HERE / "sizecorr_gen.py", "3c3a7b1edcc6e9ff003321a890420a73d5e58660e49098ff3b8d07152bdf8220"),
        "sim50_extend": (R5 / "simw3" / "sim50_extend.py", "ea3ea087620e8bcda8d08a95a80db6897fd26622ede8c5842262699eb8b602a6")}
SC_CELLS = ((1.0, 0.4), (1.0, 0.8), (0.5, 0.8), (1.0, -0.8))
SC_BASE = 100; DESIGN_CODE_SC = 13; SEEDS = tuple(range(41, 91)); ANCHOR_SEEDS = (41, 42, 43, 44, 45)
OUT_RE = re.compile(r"^sizefpr_\d{8}T\d{6}Z$")
_G = {}


class AnchorError(Exception): pass


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def load_gen():
    p, want = PINS["sizecorr_gen"]
    if sha(p) != want: raise AnchorError("sizecorr_gen differs from its pin")
    spec = importlib.util.spec_from_file_location("sizecorr_gen", p); m = importlib.util.module_from_spec(spec); sys.modules["sizecorr_gen"] = m; spec.loader.exec_module(m); return m


def register(anchor_rho0=False):
    """Per process: the three registrations of the header.  anchor_rho0: ('e12', 8) is generated by make_site_sizecorr with rho_c = 0 (anchor r only)."""
    if _G.get("done") == anchor_rho0: return _G
    import sim50_extend as S5
    m = S5.mods(); W3S = m["W3S"]; E10, E12 = m["E10"], m["E12"]; GEN = load_gen()
    orig_site, orig_spec = _G.get("orig_site", W3S.site_for), _G.get("orig_spec", W3S.cell_spec)

    def cell_spec(design, ci):
        if design != "sc": return orig_spec(design, ci)
        b = dict(E12.cells()[0]); sg, rc = SC_CELLS[ci - SC_BASE]; b.update(name=f"SC_sig{sg}_rho{rc}", size_sigma=sg, rho_c=rc); return b

    def site_for(design, ci, sd):
        if design == "sc":
            c = cell_spec(design, ci)
            return GEN.make_site_sizecorr(E10, np.random.default_rng(E12.stream(ci, sd, 0)), c["tau"], c["r"], c["pi"], c["tail"], c["F"], c["size_sigma"], c["rho_c"]), None
        if anchor_rho0 and design == "e12" and ci == 8:
            c = orig_spec(design, ci)
            return GEN.make_site_sizecorr(E10, np.random.default_rng(E12.stream(ci, sd, 0)), c["tau"], c["r"], c["pi"], c["tail"], c["F"], c["size_sigma"], 0.0), None
        return orig_site(design, ci, sd)
    W3S.site_for, W3S.cell_spec = site_for, cell_spec; W3S.DESIGN_CODE["sc"] = DESIGN_CODE_SC
    if "sc" not in S5.DESIGNS: S5.DESIGNS = tuple(S5.DESIGNS) + ("sc",)
    _G.update(done=anchor_rho0, S5=S5, m=m, GEN=GEN, orig_site=orig_site, orig_spec=orig_spec); return _G


def task(args):
    anchor_rho0, targs = args; G = register(anchor_rho0); return G["S5"].run_task(targs)


def size_fpr_corr(design, ci, sd, eps):
    """Descriptive: Pearson correlation over fields between the number of negatives and the field FPR at the uncentered threshold with site FPR closest to eps."""
    G = register(); (s, y, g), _ = G["m"]["W3S"].site_for(design, ci, sd); neg = y == 0
    t = np.quantile(s[neg], 1 - eps); n = np.bincount(g[neg]); x = np.bincount(g[neg & (s >= t)], minlength=len(n))
    k = n > 0; return float(np.corrcoef(n[k], x[k] / n[k])[0, 1])


def anchor_g(G):
    E10, E12 = G["m"]["E10"], G["m"]["E12"]
    for ci in (8, 7):
        c = G["orig_spec"]("e12", ci)
        for sd in ANCHOR_SEEDS:
            a = E12.make_site2(E10, np.random.default_rng(E12.stream(ci, sd, 0)), c["tau"], c["r"], c["pi"], c["tail"], c["F"], c["size_sigma"], c["rare"], c["slides"])[:3]
            b = G["GEN"].make_site_sizecorr(E10, np.random.default_rng(E12.stream(ci, sd, 0)), c["tau"], c["r"], c["pi"], c["tail"], c["F"], c["size_sigma"], 0.0)
            if not all(np.array_equal(u, v) for u, v in zip(a, b)): raise AnchorError(("g: rho_c = 0 differs from make_site2", ci, sd))
    return f"passed (E12 cells 8 and 7, seeds {ANCHOR_SEEDS[0]}..{ANCHOR_SEEDS[-1]}, bit for bit)"


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser(); ap.add_argument("--e12_rows"); ap.add_argument("--simw3_rows"); ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=24); ap.add_argument("--smoke", action="store_true"); a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() or not OUT_RE.match(out.name) or not out.parent.is_dir(): raise SystemExit("--out must be a new directory sizefpr_<UTC> in an existing parent")
    if not (1 <= a.workers <= 32): raise SystemExit("workers 1..32")
    t0 = time.time(); G = register(); S5 = G["S5"]; S5.mods()["E12"].gate(argv)
    if sha(PINS["sim50_extend"][0]) != PINS["sim50_extend"][1] or sha(S5.__file__) != PINS["sim50_extend"][1]: raise AnchorError("sim50_extend differs from its pin")
    S5.check_pins(); anc = dict(g=anchor_g(G)); ndraw = 4 if a.smoke else 200; seeds = (41,) if a.smoke else SEEDS
    if not a.smoke:
        for k, p in (("e12_rows", a.e12_rows), ("simw3_rows", a.simw3_rows)):
            if not p or sha(p) != S5.FILE_PINS[k]: raise AnchorError(("input file differs from its pin", k))
        e12_rows = {(r["cell"], r["seed"], r["draw"], round(r["eps"], 6), r["rule"]): r for r in json.load(open(a.e12_rows)) if r["rule"] in S5.P1}
        w3_rows = {(r["design"], r["cell"], r["seed"], r["draw"], round(r["eps"], 6), r["rule"]): r for r in json.load(open(a.simw3_rows)) if r["rule"] in S5.P2}
    with get_context("spawn").Pool(a.workers) as pool:
        if not a.smoke:
            res_a = pool.map(task, [(True, ("anchor", "e12", 8, 41, 200))])
            anc["r"] = f"passed ({S5.anchor_a(res_a, e12_rows, w3_rows)} rows of E12 cell 8 seed 41 equal the stored rows with the new generator at rho_c = 0)"
        tasks = [(False, ("new", "sc", SC_BASE + j, sd, ndraw)) for j in range(len(SC_CELLS)) for sd in seeds]
        res = []
        for i, r in enumerate(pool.imap_unordered(task, tasks), 1):
            res.append(r); print(f"task {i}/{len(tasks)} cell {r['cell']} seed {r['seed']} {r['seconds']:.0f} s", flush=True)
        corr = pool.map(_corr, [(SC_BASE + j, sd) for j in range(len(SC_CELLS)) for sd in seeds])
    bad = [(t["cell"], t["seed"]) for t in res if t["s3_mismatch"]]
    if bad: raise AnchorError(("c: pass-2 S3 rows differ from pass-1 S3 rows", bad[:5]))
    anc["c"] = "passed (pass-2 S3 rows equal pass-1 S3 rows in every task)"
    per = defaultdict(dict)
    for t in res:
        c = t["cols"]
        for ri, rule in enumerate(S5.RULES):
            for eps in S5.EPSS:
                sel = (c["rule"] == ri) & np.isclose(c["eps"], eps)
                if sel.any(): per[(t["cell"], rule, eps)][t["seed"]] = dict(S5.per_real(c["fpr"][sel], c["cert"][sel], c["sens"][sel], c["dse"][sel], c["dso"][sel], eps), _n=int(sel.sum()))
    EB = G["m"]["EB"]; summary = {}
    for (ci, rule, eps), p in sorted(per.items(), key=lambda kv: (kv[0][0], S5.RULES.index(kv[0][1]), -kv[0][2])):
        sg, rc = SC_CELLS[ci - SC_BASE]; summary[f"sc|{ci - SC_BASE}|{rule}|{eps}"] = dict(size_sigma=sg, rho_c=rc, all50=S5.block(p, EB))
    cr = defaultdict(list)
    for (ci, sd), v in corr: cr[ci - SC_BASE].append(v)
    hold = {f"sc|{j}|S3|{eps}": summary[f"sc|{j}|S3|{eps}"]["all50"]["coverage"]["mean"] >= 0.90 for j in range(len(SC_CELLS)) for eps in S5.EPSS}
    cols = {k: np.concatenate([t["cols"][k] for t in res]) for k in res[0]["cols"]}
    rep = dict(version="sizefpr_r1" + ("_SMOKE_NOT_FOR_THE_PAPER" if a.smoke else ""), script_sha256=sha(__file__), pins={k: v[1] for k, v in PINS.items()}, sim50_pins=S5.PINS, file_pins=S5.FILE_PINS,
               settings=dict(cells=[dict(index=j, size_sigma=sg, rho_c=rc, base="E12 cell 0") for j, (sg, rc) in enumerate(SC_CELLS)], seeds=list(seeds), ndraw=ndraw, m=S5.M, eps=list(S5.EPSS), workers=a.workers),
               anchors=anc, results=summary, s3_hold_1_minus_delta=hold, s3_holds_everywhere=all(hold.values()),
               size_fpr_correlation={str(j): dict(mean=float(np.mean(v)), per_realisation=v) for j, v in cr.items()},
               platt_fallback=sorted({(t["cell"], t["seed"]) for t in res if t["platt_fallback"]}), seconds=time.time() - t0,
               note="synthetic; rules and baselines reused unchanged through sim50_extend.run_task; every cell is reported")
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False); buf = io.BytesIO(); np.savez_compressed(buf, rule_names=np.array(S5.RULES), **cols)
    out.mkdir(parents=False, exist_ok=False)
    with open(out / "sizefpr_summary.json", "x") as fh: fh.write(txt)
    with open(out / "sizefpr_rows.npz", "xb") as fh: fh.write(buf.getvalue())
    with open(out / "SHA256SUMS", "x") as fh: fh.write("".join(f"{sha(out / n)}  {n}\n" for n in ("sizefpr_summary.json", "sizefpr_rows.npz")))
    print(json.dumps(dict(hold=hold, anchors=anc), indent=1)); return 0


def _corr(arg):
    ci, sd = arg; return (ci, sd), size_fpr_corr("sc", ci, sd, 0.10)


if __name__ == "__main__":
    try: sys.exit(main())
    except AnchorError as e: print("ANCHOR MISMATCH (exit 3, nothing written):", e, file=sys.stderr); sys.exit(3)

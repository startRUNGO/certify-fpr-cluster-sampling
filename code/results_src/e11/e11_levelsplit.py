# P5 E11: level-corrected B1u on the DEVELOPMENT sites (spec P5_规格_E11E12_水平修正与压力测试_开发站与合成_20261108_draft2.md, frozen in PRECHECK_LEDGER).  LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# NATURE: a supplementary analysis added after the internal blind review; it changes no rule and no reported number.  New comparison arms only:
#   S3c        = S3 with every candidate's B1u at b1u(x, n, delta/(2K)) = b1u(., ., 0.10/6): F limit and Clopper-Pearson limit each at delta/6, union over the K = 3 candidates, same grid, fixed sequence, candidate choice as S3;
#   B1u_c1s_c  = B1u_c1s with b1u(., ., 0.10/2) (each step delta/2, single candidate c1s).
#   Controls, same draws (common random numbers): S3 (b1u at delta/3), B1u_c0 and B1u_c1s (b1u at delta), recomputed by the copy of the frozen REV-1 replay below.
# Data boundary: development sites only (bbbc041_all, mpidb_full, lacuna_B), arm a2, seeds 41-45, 200 draws per seed, m = 22, eps in {0.10, 0.05}, delta = 0.10.  No PlasmoCount, IML, M5, NIH test key.  Gate (exit 4, before any file is read):
#   any argument containing plasmo / iml / r8 / m5 (any case), or any already-imported module named m5_* / e8_* / r8_* / e7_*.
# Anchors (exit 3, no output directory):
#   lib   : m3_replay_lib sha pin (ed9d1f17...); L.Inputs pins every input file (as the frozen drivers do).
#   A2    : the S3, B1u_c0 and B1u_c1s rows of THIS copy equal the stored REV-1 aggregates within 1e-12 (L.self_check, 3 sites x 2 eps).
#   A3    : on every 10th draw (d % 10 == 0), for every candidate and every grid point, b1u at 0.10/6 >= b1u at 0.10/3 (and 0.10/2 >= 0.10) - the monotonicity in the level (a smaller level gives a larger upper bound); any violation exits 3.
#           Draws where either bound raises are counted (reported), not compared.
#   A4    : per draw and eps, S3c certified => S3 certified; B1u_c1s_c certified => B1u_c1s certified (a more conservative bound cannot certify more).
# Reporting (per site, eps, rule): n, coverage = share of draws with realised FPR <= eps, coverage per seed and its 5-seed 90% interval (np.percentile 5/95) with min and max, band, near_boundary, fpr_mean, cannot_certify, sens, ds_eps, ds_own.
#   band: coverage >= 0.90 'hold_1-delta'; 0.80 <= coverage < 0.90 'below_1-delta_within_1-2delta'; < 0.80 'fail'.  near_boundary: within 0.01 of 0.90 or of 0.80.
#   Paired (S3c - S3, B1u_c1s_c - B1u_c1s) per site and eps: delta_coverage (pooled), delta_ds_own (pooled mean of per-draw differences), per-seed means and the number of seeds with a negative delta_ds_own.
#   reading_flags: data only, no sentence: material_cost (eps = 0.10: at least 2 of 3 sites with delta_ds_own <= -0.03), material_gain per site and eps (S3 in a lower band than S3c, or delta_coverage >= 0.03).
# Output: one NEW directory with e11_summary.json, e11_rows.json (per draw, without ni and fields), SHA256SUMS (files opened with 'x', allow_nan=False).  Prints the directory only.
# Usage: python -I -B e11_levelsplit.py --a2link .. --ids_iter1 .. --mhlh .. --rev1_bbbc .. --rev1_mh .. --rev1_lh .. --mh_rev1_py .. --out <new dir>
import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ.setdefault(_k, "8")
import argparse, hashlib, json, re, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; R = HERE.parent
FORBID_ARG = re.compile(r"plasmo|iml|r8|m5", re.I); FORBID_MOD = re.compile(r"^(m5_|e8_|r8_|e7_)")
PIN_LIB = "ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de"
LEVEL = {"U1": 0.10, "U3": 0.10 / 3, "U1c": 0.10 / 2, "U3c": 0.10 / 6}   # b1u level per rule name


class AnchorError(Exception): pass


def gate(argv):
    """Exit 4 on any confirmatory-site string in the arguments or any confirmatory-site module already imported."""
    bad = [a for a in argv if FORBID_ARG.search(a)]
    if bad: print("GATE (exit 4): confirmatory-site string in argument:", bad[:3], file=sys.stderr); sys.exit(4)
    mods = [m for m in sys.modules if FORBID_MOD.match(m.split(".")[-1])]
    if mods: print("GATE (exit 4): confirmatory-site module loaded:", mods[:3], file=sys.stderr); sys.exit(4)


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def replay_e11(L, inp, sd, site, a3):
    """Copy of m3_replay_lib.replay (draw loop, streams, centring, grid, fixed sequence, labelled choice) with two added rule families (U3c -> S3c, U1c -> B1u_c1s_c); rows drop ni and fields.
    a3: dict with counters 'checked', 'violations', 'exceptions' (updated in place)."""
    pop = L.population(inp, sd, site); fl, fl_neg, fl_posv = pop["fl"], pop["fl_neg"], pop["fl_posv"]; sens_at, fpr_at = pop["sens_at"], pop["fpr_at"]
    key8, sidx = L.SITECFG[site]; ai = L.AI; rows = []; m = L.M
    rng = np.random.default_rng([20261022, key8, ai, sd, m, sidx])
    for d_ in range(L.NDRAW):
        sel_ = rng.choice(len(fl), m, replace=False); F = [fl[i] for i in sel_]; clicks = sum(len(fl_posv["c0"][f]) for f in F); ni = np.array([len(fl_neg["c0"][f]) for f in F])
        for eps, dl in L.SETTINGS:   # kept: consumes the stream exactly as REV-1 (R1-2 subsets)
            r_ = np.random.default_rng(int(rng.integers(1 << 30))); [r_.choice(n_, 2, replace=False) if n_ >= 2 else np.arange(n_) for n_ in ni]
        D = {}
        for k in L.CEN:
            negs = [fl_neg[k][f] for f in F]; pooled = np.concatenate(negs); grid = np.unique(np.quantile(pooled, 1-np.arange(L.G)/L.G))[::-1]
            X = np.stack([ni[i]-np.searchsorted(negs[i], grid, side="left") for i in range(m)]).astype(float)
            posl = np.concatenate([fl_posv[k][f] for f in F]) if clicks > 0 else np.array([]); D[k] = (negs, pooled, grid, X, posl)
        cache = {}

        def ubf(rule, k, j):
            if (rule, k, j) not in cache: cache[(rule, k, j)] = L.b1u(D[k][3][:, j], ni, LEVEL[rule])
            return cache[(rule, k, j)]
        if d_ % 10 == 0:     # A3: monotonicity in the level, on every grid point of every candidate
            for k in L.CEN:
                for j in range(len(D[k][2])):
                    for hi_rule, lo_rule in (("U3c", "U3"), ("U1c", "U1")):
                        try: u_lo, u_hi = ubf(lo_rule, k, j), ubf(hi_rule, k, j)
                        except Exception: a3["exceptions"] += 1; continue
                        if not (np.isfinite(u_lo) and np.isfinite(u_hi)): a3["exceptions"] += 1; continue
                        a3["checked"] += 1
                        if u_hi < u_lo - 1e-12: a3["violations"] += 1; a3.setdefault("first", [site, sd, d_, k, j, hi_rule, float(u_hi), float(u_lo)])

        def rec(name, eps, k, t, cert):
            sn = sens_at(k, t)
            rows.append(dict(site=site, seed=sd, draw=d_, eps=eps, rule=name, cand=k, cert=bool(cert), fpr=float(fpr_at(k, t)), sens=float(sn),
                             dse=float(sn-pop["orac_unc"][eps]), dso=float(sn-pop["orac_own"][(k, eps)])))

        def s3_like(name, rule, eps):
            T = {k: L.fixed_sequence(D[k][2], lambda j, k=k: ubf(rule, k, j), eps) for k in L.CEN}; cert_c = [k for k in L.CEN if T[k][1]]
            if not cert_c: rec(name, eps, "c0", T["c0"][0], False)
            else:
                if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
                else:
                    sl = [(D[k][4] >= T[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
                rec(name, eps, ch, T[ch][0], True)
        for eps, dl in L.SETTINGS:
            s3_like("S3", "U3", eps)
            for k in ("c0", "c1s"):
                t_, c_, _ = L.fixed_sequence(D[k][2], lambda j, k=k: ubf("U1", k, j), eps); rec("B1u_" + k, eps, k, t_, c_)
            s3_like("S3c", "U3c", eps)
            t_, c_, _ = L.fixed_sequence(D["c1s"][2], lambda j: ubf("U1c", "c1s", j), eps); rec("B1u_c1s_c", eps, "c1s", t_, c_)
    return rows


def band(c):
    return "hold_1-delta" if c >= 0.90 else ("below_1-delta_within_1-2delta" if c >= 0.80 else "fail")


def near(c): return bool(abs(c - 0.90) <= 0.01 or abs(c - 0.80) <= 0.01)


def cell(rows_c, eps):
    fpr = np.array([r["fpr"] for r in rows_c]); cov = float((fpr <= eps).mean()); seeds = sorted({r["seed"] for r in rows_c})
    cps = [float(np.mean([r["fpr"] <= eps for r in rows_c if r["seed"] == s])) for s in seeds]
    return dict(n=len(rows_c), coverage=cov, coverage_per_seed=cps, coverage_ci90_5seeds=[float(np.percentile(cps, 5)), float(np.percentile(cps, 95))], coverage_min=min(cps), coverage_max=max(cps),
                band=band(cov), near_boundary=near(cov), fpr_mean=float(fpr.mean()), cannot_certify=float(1 - np.mean([r["cert"] for r in rows_c])), sens=float(np.mean([r["sens"] for r in rows_c])),
                ds_eps=float(np.mean([r["dse"] for r in rows_c])), ds_own=float(np.mean([r["dso"] for r in rows_c])))


def paired(rows, site, eps, a, b):
    ra = {(r["seed"], r["draw"]): r for r in rows if r["site"] == site and r["eps"] == eps and r["rule"] == a}
    rb = {(r["seed"], r["draw"]): r for r in rows if r["site"] == site and r["eps"] == eps and r["rule"] == b}
    assert ra.keys() == rb.keys(), ("unpaired draws", site, eps, a, b); keys = sorted(ra)
    dc = float(np.mean([(ra[k]["fpr"] <= eps) for k in keys]) - np.mean([(rb[k]["fpr"] <= eps) for k in keys])); d = np.array([ra[k]["dso"] - rb[k]["dso"] for k in keys])
    per = [float(np.mean([ra[k]["dso"] - rb[k]["dso"] for k in keys if k[0] == s])) for s in sorted({k[0] for k in keys})]
    return dict(rule=a, vs=b, delta_coverage=dc, delta_ds_own=float(d.mean()), delta_ds_own_per_seed=per, seeds_with_negative_delta_ds_own=int(sum(x < 0 for x in per)))


def check_a4(rows):
    cert = {(r["site"], r["seed"], r["draw"], r["eps"], r["rule"]): r["cert"] for r in rows}; bad = []
    for (site, sd, d_, eps, rule), c in cert.items():
        base = {"S3c": "S3", "B1u_c1s_c": "B1u_c1s"}.get(rule)
        if base and c and not cert[(site, sd, d_, eps, base)]: bad.append((site, sd, d_, eps, rule))
    return bad


def main(argv=None, _hooks=None):
    """_hooks: test seam used only by e11_smoke_local.py ({'inputs': factory(L, args) -> inputs object, 'self_check': fn(inp, rows) -> list}); None in every real run."""
    argv = sys.argv[1:] if argv is None else argv; gate(argv)
    ap = argparse.ArgumentParser()
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(argv); out = Path(a.out)
    if out.exists(): raise SystemExit("output directory exists: use a new directory")
    sys.path[:0] = [str(R / "m3" / "k2prep")]
    import m3_replay_lib as L
    gate(argv)
    if not _hooks and L.sha(L.__file__) != PIN_LIB: raise AnchorError(("m3_replay_lib differs from its pin", L.sha(L.__file__)))
    inp = _hooks["inputs"](L, a) if _hooks else L.Inputs(a.a2link, a.ids_iter1, a.mhlh, a.rev1_bbbc, a.rev1_mh, a.rev1_lh, a.mh_rev1_py)
    t0 = time.time(); a3 = dict(checked=0, violations=0, exceptions=0); rows = []
    for site in L.DEV:
        for sd in L.SEEDS: rows += replay_e11(L, inp, sd, site, a3)
    bad = (_hooks["self_check"] if _hooks else L.self_check)(inp, rows)
    if bad: raise AnchorError(("A2 replay self-check vs stored REV-1 failed", bad[:5]))
    if a3["violations"]: raise AnchorError(("A3 monotonicity in the level violated", a3))
    bad4 = check_a4(rows)
    if bad4: raise AnchorError(("A4 a more conservative rule certified a draw its control did not", bad4[:5]))
    rules = ("S3", "S3c", "B1u_c1s", "B1u_c1s_c"); cells = {}; pairs = []
    for site in L.DEV:
        for eps, _ in L.SETTINGS:
            for rule in rules: cells[f"{site}|{eps}|{rule}"] = dict(site=site, eps=eps, rule=rule, **cell([r for r in rows if r["site"] == site and r["eps"] == eps and r["rule"] == rule], eps))
            for x, y in (("S3c", "S3"), ("B1u_c1s_c", "B1u_c1s")): pairs.append(dict(site=site, eps=eps, **paired(rows, site, eps, x, y)))
    cost_sites = [p["site"] for p in pairs if p["eps"] == 0.10 and p["rule"] == "S3c" and p["delta_ds_own"] <= -0.03]
    order = {"fail": 0, "below_1-delta_within_1-2delta": 1, "hold_1-delta": 2}
    gain = [dict(site=p["site"], eps=p["eps"], band_S3=cells[f"{p['site']}|{p['eps']}|S3"]["band"], band_S3c=cells[f"{p['site']}|{p['eps']}|S3c"]["band"], delta_coverage=p["delta_coverage"],
                 material=bool(order[cells[f"{p['site']}|{p['eps']}|S3"]["band"]] < order[cells[f"{p['site']}|{p['eps']}|S3c"]["band"]] or p["delta_coverage"] >= 0.03))
            for p in pairs if p["rule"] == "S3c"]
    rep = dict(version="e11_r1", spec="P5_规格_E11E12_水平修正与压力测试_开发站与合成_20261108_draft2", script_sha256=sha(__file__), lib_sha256=L.sha(L.__file__), inputs_sha256=dict(L.PIN), score_files_sha256=inp.used,
               levels=LEVEL, anchors=dict(A2="passed (S3, B1u_c0, B1u_c1s vs stored REV-1 within 1e-12)", A3=a3, A4="passed"), cells=cells, paired=pairs,
               reading_flags=dict(material_cost_eps010=dict(sites_with_delta_ds_own_le_minus_0p03=cost_sites, flag=bool(len(cost_sites) >= 2)), material_gain=gain), seconds=time.time() - t0,
               note="supplementary analysis; changes no rule and no reported number; the sealed S3/SPLIT6 rules and confirmations are unchanged")
    out.mkdir(parents=False, exist_ok=False)
    files = {"e11_summary.json": json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False), "e11_rows.json": json.dumps(rows, allow_nan=False)}
    for n, txt in files.items():
        with open(out / n, "x") as fh: fh.write(txt)
    with open(out / "SHA256SUMS", "x") as fh: fh.write("".join(f"{sha(out / n)}  {n}\n" for n in files))
    print(out)


if __name__ == "__main__":
    try: main()
    except AnchorError as e: print("ANCHOR MISMATCH (exit 3, no output):", e, file=sys.stderr); sys.exit(3)
    except AssertionError as e: print("PINNED INPUT / LIBRARY CHECK FAILED (exit 3, no output):", e, file=sys.stderr); sys.exit(3)   # L.Inputs pin and path assertions

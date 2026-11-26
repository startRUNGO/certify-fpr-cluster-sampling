# P5 K1vB development-site evaluation (problem card draft6 f859dd4e, section 3 K1vB and section 4 criteria; draft5 b3318a71 + addendum A 6e16d981 for K1).
# LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 approves.
# NATURE: development-site mechanism selection on 4 sites (bbbc041_all, mpidb_full rebuilt, lacuna_B rebuilt, old MP-IDB); not confirmatory evidence; no M5, IML,
#   PlasmoCount, NIH key.
# ARMS: S3 (REV-1), B1u_c1s (= single-candidate c1s), K1 (draft5 K1+S3, alpha = delta/3 per candidate), and the new main arm K1vB:
#   for each candidate k (c0, c1s, c2s) and grid point t: "certified at t" iff  b1u(x, n, delta/6) <= eps  OR  UB_Z^{K1}(t; alpha = delta/6) <= 0;
#   fixed-sequence test on the S3 grid; union over the 3 candidates; choice as S3.  Validity per t: Pr(false) <= delta/6 + delta/6 = delta/3 (Bonferroni; no
#   independence needed); fixed sequence + union -> delta.  b1u is evaluated first; K1 only if b1u does not certify (same result, fewer evaluations).
# K1 at a given alpha: k1_ub_alpha copies m3_k12_dev_eval.k1_ub (0091a333) with alpha as a parameter (degenerate tolerance, drop > 5% -> +inf, bootstrap-t with
#   B = 999, weights default_rng([20261101, site index, seed, draw, eps code, candidate code])); the smoke asserts k1_ub_alpha(alpha = delta/3) == k1_ub exactly.
#   Counters (exceptions, se0, drop > 5%) are written into the evaluation module's EXC / DEG / DROP so that its summarise() reports them.
# DATA: development sites via m3_replay_lib; old MP-IDB via m3_stress.StressInputs (seed 41 local iter1 file, only key mpidb_full_logits, sha = MANIFEST
#   source_sha256; seeds 42-45 from runs/a2_iter1_mpidb_old_20261030T163124Z checked against SHA256SUMS; ids_iter1 keys mpidb_full_y/group only).
# SELF-CHECK before anything is written (1e-12): S3, B1u_c0, B1u_c1s of this loop vs stored REV-1 for all 4 sites (old MP-IDB: rev1_results.json "a2|mpidb_full|22|...").
# CRITERIA (draft6 section 4, fixed before any run) for K1vB against S3:
#   gate, every one of 8 cells (4 sites x 2 eps): coverage >= 0.90, seed-cluster bootstrap 90% lower limit >= 0.88, cannot_certify <= S3 + 0.02, no flagged cell;
#   main: mean dSens_own improvement over the 3 rebuilt sites at eps 0.10 >= +0.02 and none of those 6 cells worse than S3 by more than 0.01;
#         old MP-IDB: both eps cells not worse than S3 by more than 0.01.
#   Fail -> back to root-cause analysis; criteria and the delta split are not changed.
# Output: one new JSON in a new directory ("x", allow_nan=False).  Prints only the output path.
# Usage: python -I -B m3_k1vb_dev_eval.py --a2link .. --ids_iter1 .. --mhlh .. --rev1_bbbc .. --rev1_mh .. --rev1_lh .. --mh_rev1_py .. --old_s41 .. --old_dir .. --out ..
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "k2prep")); sys.path.insert(0, str(ROOT / "stress"))
import m3_replay_lib as L   # noqa: E402
import m3_k12_dev_eval as E   # noqa: E402
import m3_stress as ST   # noqa: E402

DELTA = 0.10; A_K1 = DELTA/3; A_VB = DELTA/6; OLD = ST.OLD; REBUILT = ("bbbc041_all", "mpidb_full", "lacuna_B"); ARMS = ("S3", "B1u_c1s", "K1", "K1vB")
GATE_COV, GATE_LO, GATE_CC, MAIN_GAIN, MAIN_WORSE, EXC_SHARE = 0.90, 0.88, 0.02, 0.02, 0.01, 0.01


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def k1_ub_alpha(e, Sw, F, m, W, alpha, key):
    """m3_k12_dev_eval.k1_ub with alpha as a parameter (bootstrap-t form only)."""
    eb = e.mean(); se = e.std(ddof=1); fpc = np.sqrt(1 - m/F); Zh = Sw + F*eb; tol = E.degenerate_tol(e)
    if se <= tol: E.DEG[key]["se0"] += 1; return np.inf
    SE = F*se/np.sqrt(m)*fpc
    mu = (W*e).sum(1)/m; ex2 = (W*e*e).sum(1)/m; s2 = (ex2 - mu**2)*m/(m - 1); ok = s2 > tol**2; drop = 1 - ok.mean(); E.DROP[key].append(drop)
    if drop > E.MAX_DROP: E.DEG[key]["drop>5%"] += 1; return np.inf
    T = (mu[ok] - eb)/(np.sqrt(s2[ok])/np.sqrt(m)); return Zh - np.quantile(T, alpha)*SE


def run_site_seed(inp, sd, site):
    pop = L.population(inp, sd, site); fl, fl_neg, fl_posv = pop["fl"], pop["fl_neg"], pop["fl_posv"]; sens_at, fpr_at = pop["sens_at"], pop["fpr_at"]
    AC = {k: E.AllCells(pop, k) for k in L.CEN}
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

        def ubf(level, k, j):
            if (level, k, j) not in cache: cache[(level, k, j)] = L.b1u(D[k][3][:, j], ni, level)
            return cache[(level, k, j)]

        def rec(name, eps, k, t, cert):
            sn = sens_at(k, t)
            rows.append(dict(site=site, seed=sd, draw=d_, eps=eps, rule=name, cand=k, cert=bool(cert), fpr=float(fpr_at(k, t)), sens=float(sn),
                             dse=float(sn-pop["orac_unc"][eps]), dso=float(sn-pop["orac_own"][(k, eps)]), ni=ni.tolist(), fields=[str(f) for f in F]))

        def choose(T, eps, name):
            cert_c = [k for k in L.CEN if T[k][1]]
            if not cert_c: rec(name, eps, "c0", T["c0"][0], False); return
            if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
            else:
                sl = [(D[k][4] >= T[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
            rec(name, eps, ch, T[ch][0], True)
        for eps, dl in L.SETTINGS:
            choose({k: L.fixed_sequence(D[k][2], lambda j, k=k: ubf(0.10/3, k, j), eps) for k in L.CEN}, eps, "S3")
            for k in ("c0", "c1s"):
                t_, c_, _ = L.fixed_sequence(D[k][2], lambda j, k=k: ubf(0.10, k, j), eps); rec("B1u_" + k, eps, k, t_, c_)
            Wk = {k: np.random.default_rng([E.BOOT_SEED, L.DEV.index(site), sd, d_, int(round(eps*100)), E.CAND_CODE[k]]).multinomial(m, np.full(m, 1/m), size=E.B).astype(float)
                  for k in L.CEN}
            cl = {k: AC[k].c_lab(F) for k in L.CEN}

            def k1_U(j, k, alpha, arm, eps=eps):
                t = D[k][2][j]; z = D[k][3][:, j] - eps*ni; w = AC[k].g_lab(F, t) - eps*cl[k]; Sw, _ = AC[k].w_stats(t, eps)
                return k1_ub_alpha(z - w, Sw, AC[k].F, m, Wk[k], alpha, (site, eps, arm))

            def ub_k1(j, k): return E.to_seq(k1_U(j, k, A_K1, "K1"))

            def ub_vb(j, k, eps=eps):
                if ubf(A_VB, k, j) <= eps: return 0.0
                return E.to_seq(k1_U(j, k, A_VB, "K1vB"))
            choose({k: L.fixed_sequence(D[k][2], E.counted("K1", site, eps, lambda j, k=k: ub_k1(j, k)), eps) for k in L.CEN}, eps, "K1")
            choose({k: L.fixed_sequence(D[k][2], E.counted("K1vB", site, eps, lambda j, k=k: ub_vb(j, k)), eps) for k in L.CEN}, eps, "K1vB")
    return rows


def judge(cells):
    get = lambda arm, s, e: [c for c in cells if c["arm"] == arm and c["site"] == s and abs(c["eps"] - e) < 1e-12][0]
    gate, diffs = [], []
    for s in REBUILT + (OLD,):
        for e, _ in L.SETTINGS:
            a, b = get("K1vB", s, e), get("S3", s, e)
            gate.append(dict(site=s, eps=e, cov_ok=a["coverage"] >= GATE_COV, lo_ok=a["coverage_ci90_seedcluster"][0] >= GATE_LO,
                             cc_ok=a["cannot_certify"] <= b["cannot_certify"] + GATE_CC + 1e-12, flagged=a["flag"] is not None))
            diffs.append(dict(site=s, eps=e, diff=a["ds_own"] - b["ds_own"]))
    g_ok = all(x["cov_ok"] and x["lo_ok"] and x["cc_ok"] and not x["flagged"] for x in gate)
    mg = float(np.mean([d["diff"] for d in diffs if d["site"] in REBUILT and abs(d["eps"] - 0.10) < 1e-12]))
    rebuilt_ok = all(d["diff"] >= -MAIN_WORSE for d in diffs if d["site"] in REBUILT); old_ok = all(d["diff"] >= -MAIN_WORSE for d in diffs if d["site"] == OLD)
    m_ok = mg >= MAIN_GAIN and rebuilt_ok and old_ok
    return dict(gate=gate, gate_pass=g_ok, mean_gain_rebuilt_eps010=mg, per_cell_diff=diffs, rebuilt_no_cell_worse=rebuilt_ok, old_not_worse=old_ok,
                main_pass=m_ok, passes=bool(g_ok and m_ok),
                if_fail="back to root-cause analysis; criteria and the delta split are not changed")


def main():
    ap = argparse.ArgumentParser()
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "old_s41", "old_dir", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); out = Path(a.out); assert not out.exists(), "output exists"
    inp = ST.StressInputs(a); assert sorted(inp.old) == [41, 42, 43, 44, 45], "old MP-IDB needs all 5 seeds"
    L.SITECFG[OLD] = (8, 1); L.DEV = tuple(L.DEV) + (OLD,); L.M = 22
    E.EXC.clear(); E.DEG.clear(); E.DROP.clear(); rows = []
    for site in REBUILT + (OLD,):
        for sd in L.SEEDS: rows += run_site_seed(inp, sd, site)
    status, bad = ST.check(inp, rows, REBUILT + (OLD,), 22)
    if bad or not all(v.startswith("self-check passed") for v in status.values()):
        raise SystemExit("self-check failed or not possible, nothing written: " + repr((bad[:5], status)))
    rng = np.random.default_rng(E.CI_SEED)
    cells = [dict(arm=arm, site=s, eps=e, **E.summarise([r for r in rows if r["site"] == s and abs(r["eps"] - e) < 1e-12 and r["rule"] == arm], e, s, arm, rng))
             for arm in ARMS for s in REBUILT + (OLD,) for e, _ in L.SETTINGS]
    rep = dict(script_sha256=sha(__file__), lib_sha256=sha(L.__file__), eval_sha256=sha(E.__file__), stress_sha256=sha(ST.__file__), inputs_sha256=dict(L.PIN),
               score_files_sha256=inp.used, problem_card="draft6 f859dd4e (K1vB); draft5 b3318a71 + addendum A 6e16d981 (K1)",
               nature="development-site mechanism selection on 4 sites; not confirmatory evidence",
               settings=dict(delta=DELTA, alpha_K1=A_K1, alpha_K1vB_each=A_VB, B=E.B, boot_seed=E.BOOT_SEED, ci_B=E.CI_B, ci_seed=E.CI_SEED),
               self_check=status, cells=cells, judgement=judge(cells))
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False)
    out.parent.mkdir(parents=False, exist_ok=False)
    with open(out, "x") as fh: fh.write(txt)
    print("written", out)


if __name__ == "__main__":
    main()

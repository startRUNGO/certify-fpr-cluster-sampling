# P5 K2 preparation R3 (problem card draft4 fcf90632, section 6): design-side UPPER BOUND on what reducing the across-draw dispersion of the certified FPR could buy,
# for the two K2 arms separately: base S3 (arm K2+S3) and base single-candidate c1s = B1u_c1s (arm K2+c1s).  LOCAL, read-only.  WRITTEN, NOT RUN, until approved.
# STOP RULE (fixed by 实验管理 2026-10-30 and filed in draft4 before this script was written; applied per arm, each arm against its own base):
#   if the three-site mean improvement of dSens_own at eps 0.10 with sigma scaled by kappa = 0.5 is < +0.03 -> stop K2 on that arm (as B-round M1);
#   if >= +0.03 -> the K2 main criterion stays +0.02, it is NOT adjusted to the size of the bound.
# MODEL OF "SIGMA SCALED BY kappa, COVERAGE UNCHANGED" (operational definition, stated here before any number is computed):
#   for each certified draw d of the base rule, realized FPR fpr_d is replaced by fpr'_d = eps - kappa (eps - fpr_d) (contraction toward eps): every draw stays on the
#   same side of eps, so coverage is exactly unchanged, and the spread of fpr' is kappa times that of fpr.  Its sensitivity is the population sensitivity of the
#   same candidate at FPR fpr'_d: sens'_d = sens_at(k, quantile(negatives_k, 1 - fpr'_d)).  Uncertified draws keep their fallback values.
#   improvement = mean_d sens'_d(kappa) - mean_d sens'_d(1): both terms go through the same quantile path, so the discretisation error of that path largely cancels
#   (revised before any real run, because the smoke showed the discretisation error; 实验管理 2026-10-30).  The own oracle of the chosen candidate cancels too, so this
#   is the change in dSens_own.
#   This is the gain of a hypothetical rule with the same coverage and kappa times the dispersion; it is not a proof that K2 attains it.
# CHECKS before anything is written (exit, nothing written):
#   (1) replay self-check: S3, B1u_c0, B1u_c1s equal the stored REV-1 aggregates (cov, fpr_mean, cannot_certify, ds_eps, ds_own, n) within 1e-12 (m3_replay_lib);
#   (2) discretisation SANITY check (not a precision guarantee for the improvement, which uses sens'(1) differences): per cell, over certified draws,
#       mean_d |sens'_d(1) - sens_d| <= 1/n_pos (n_pos = positives of the site); otherwise exit, nothing written.  Also reported per cell, without exit:
#       max_d |sens'_d(1) - sens_d|, flagged when > 3/n_pos (several positives can sit in one negative-score gap, e.g. ties).
# Reported per arm, site, eps, kappa in {0.5, 0.7}: improvement, base mean dSens_own, base and contracted fpr sd; and the stop-rule reading per arm.
# Output: one new JSON in a new directory ("x", allow_nan=False).  Prints only the output path.
# Usage: python -I -B m3_r3_bound.py --a2link .. --ids_iter1 .. --mhlh .. --rev1_bbbc .. --rev1_mh .. --rev1_lh .. --mh_rev1_py .. --out <dir/json>
import argparse, json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent)); import m3_replay_lib as L   # noqa: E402

ARMS = {"K2+S3": "S3", "K2+c1s": "B1u_c1s"}; KAPPAS = (0.5, 0.7); STOP_AT = 0.03; MAIN_CRIT = 0.02; MAX_FLAG = 3


def contracted(rows, pops, eps, kappa):
    out = []
    for r in rows:
        if not r["cert"]: out.append(r["sens"]); continue
        pop = pops[(r["site"], r["seed"])]; k = r["cand"]; fp = eps - kappa*(eps - r["fpr"])
        out.append(float(pop["sens_at"](k, np.quantile(pop["sneg"][k], 1 - fp))))
    return np.asarray(out)


def main():
    ap = argparse.ArgumentParser()
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); out = Path(a.out); assert not out.exists(), "output exists"
    inp = L.Inputs(a.a2link, a.ids_iter1, a.mhlh, a.rev1_bbbc, a.rev1_mh, a.rev1_lh, a.mh_rev1_py)
    rows, pops = [], {}
    for site in L.DEV:
        for sd in L.SEEDS:
            r_, pop = L.replay(inp, sd, site); rows += r_; pops[(site, sd)] = dict(sens_at=pop["sens_at"], sneg=pop["sneg"], n_pos=len(pop["spos"]["c0"]))
    bad = L.self_check(inp, rows)
    if bad: raise SystemExit("self-check failed, nothing written: " + repr(bad[:5]))
    cells, disc = [], []
    for arm, base in ARMS.items():
        for site in L.DEV:
            for eps, _ in L.SETTINGS:
                R = [r for r in rows if r["site"] == site and r["eps"] == eps and r["rule"] == base]; sens = np.array([r["sens"] for r in R])
                fpr = np.array([r["fpr"] for r in R]); cert = np.array([r["cert"] for r in R]); s1 = contracted(R, pops, eps, 1.0)
                npos = min(pops[(site, sd)]["n_pos"] for sd in L.SEEDS); dev = np.abs(s1 - sens)[cert]
                d_mean = float(dev.mean()) if cert.any() else 0.0; d_max = float(dev.max()) if cert.any() else 0.0
                disc.append((arm, site, eps, d_mean, 1.0/npos))
                c = dict(arm=arm, base=base, site=site, eps=eps, n=len(R), base_ds_own=float(np.mean([r["dso"] for r in R])), base_fpr_sd_certified=(float(fpr[cert].std()) if cert.any() else None),
                         n_pos=int(npos), discretisation_mean_abs_kappa1=d_mean, discretisation_tol=1.0/npos, discretisation_max_abs_kappa1=d_max,
                         discretisation_max_flag=("max > 3/n_pos" if d_max > MAX_FLAG/npos else None), by_kappa={})
                for kp in KAPPAS:
                    sp = contracted(R, pops, eps, kp); fp2 = np.where(cert, eps - kp*(eps - fpr), fpr)
                    c["by_kappa"][str(kp)] = dict(improvement=float(sp.mean() - s1.mean()), fpr_sd_certified=(float(fp2[cert].std()) if cert.any() else None),
                                                  coverage_base=float((fpr <= eps).mean()), coverage_contracted=float((fp2 <= eps).mean()))
                cells.append(c)
    fails = [(a_, s_, e_, d_, t_) for a_, s_, e_, d_, t_ in disc if d_ > t_]
    if fails: raise SystemExit("discretisation check failed (mean |sens'(1) - sens| > 1/n_pos): " + repr(fails[:5]) + "; nothing written")
    reading = {}
    for arm in ARMS:
        imp = [c["by_kappa"]["0.5"]["improvement"] for c in cells if c["arm"] == arm and abs(c["eps"] - 0.10) < 1e-12]; ub = float(np.mean(imp))
        reading[arm] = dict(mean_improvement_eps010_kappa05=ub, rule=f"< {STOP_AT} stop K2 on this arm; >= {STOP_AT} main criterion stays +{MAIN_CRIT}",
                            verdict=("stop K2 on this arm" if ub < STOP_AT else f"continue; main criterion stays +{MAIN_CRIT} (not adjusted)"))
    rep = dict(script_sha256=L.sha(__file__), lib_sha256=L.sha(L.__file__), inputs_sha256=dict(L.PIN), score_files_sha256=inp.used,
               nature="design-side upper bound (problem card draft4 R3); stop rule fixed before the script was written; not confirmatory evidence",
               model="fpr' = eps - kappa (eps - fpr) on certified draws; coverage unchanged; sens' from the population ROC of the chosen candidate; improvement = mean sens'(kappa) - mean sens'(1)",
               wording="model gain at the same coverage and kappa times the dispersion; not a strict upper bound and not a proof that K2 attains it",
               self_check="passed (S3, B1u_c0, B1u_c1s vs stored REV-1 within 1e-12)", discretisation_check="passed (mean |sens'(1) - sens| <= 1/n_pos in every cell)",
               cells=cells, reading=reading)
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False)
    out.parent.mkdir(parents=False, exist_ok=False)
    with open(out, "x") as fh: fh.write(txt)
    print("written", out)


if __name__ == "__main__":
    main()

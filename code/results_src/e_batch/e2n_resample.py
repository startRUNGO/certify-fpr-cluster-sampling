# P5 E2-n small-calibration-set resampling (spec E-batch draft3 section E2-n; answers Kladny 2026: is a calibration-conditional guarantee true at small n?).  LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# NATURE: descriptive curves of the violation rate against the number n of labeled NEGATIVE CELLS on the three DEVELOPMENT sites (A2 scores, seeds 41..45).  No pass/fail on any rule.
#   M5, IML, PlasmoCount, NIH test, any image, any classifier: never read.  N0 is not here (it does not depend on n; see E1).  S3 and SPLIT6 are not here (they take the number of fields m as design parameter; E2' covers m).
# n in {10, 20, 50, 100, 200}; eps in {0.10, 0.05}; rules Q, MC, CP (definitions in ebatch_lib: Q the ceil(n(1-eps))-th smallest, MC the ceil((n+1)(1-eps))-th smallest or vacuous, CP the exact threshold with at most
#   x* sample false positives, x* the largest count whose Clopper-Pearson upper limit at level 0.90 is <= eps; x* = 0 puts tau just above the sample maximum; none -> vacuous).
# Modes: iid = n negatives drawn WITHOUT replacement from all negatives of the site (exchangeable cells); field = fields in random order, cells of a field in random order, cumulated until >= n, truncated to n (cluster structure kept).
# Stream per (site, seed, n, draw): default_rng([20261103, 82, site index, seed, n index, draw]), created afresh for each mode (the two modes therefore start from the same stream key); NDRAW = 1000.
# Truth: FPR of tau on all negatives of the site; violation = FPR > eps.  coverage = 1 - violation rate.
# (M2) Implementation check, iid mode: for every rule, n, eps that is not vacuous, the violation rate pooled over the 3 sites and 5 seeds must equal the EXACT value P(Beta(n+1-k, k) > eps)
#   (ebatch_lib.viol_analytic) within 3 Monte-Carlo standard errors + 2/N_min (finite-population slack); CP must also be <= delta within the same tolerance.  Any miss -> exit 5, nothing is read, no output directory.
# The analytic law assumes continuous scores: the share of repeated negative scores per site (1 - unique/N) is stored in the output as `duplicate_fraction` and printed with an exit-5 message; look at it FIRST if the check fires.
# Non-vacuous convention: coverage_nonvacuous / fpr_mean_nonvacuous exclude draws with tau = +inf (see ebatch_lib.VAC_CONVENTION); E2-b uses them only.
# Anchors: master pin and the AUROC anchor of E1.  Usage as E1 without --nihval_dir.
import argparse, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ebatch_lib as E

NS = (10, 20, 50, 100, 200); NDRAW = 1000; MODES = ("iid", "field"); TAG = 82


def draw_iid(pop, rng, n): return np.sort(pop["sneg"][rng.choice(pop["N"], n, replace=False)])


def draw_field(pop, rng, n):
    fl = pop["fl"]; got = []; c = 0
    for i in rng.permutation(len(fl)):
        v = pop["negf"][fl[i]]; got.append(v[rng.permutation(len(v))]); c += len(v)
        if c >= n: break
    return np.sort(np.concatenate(got)[:n]) if c >= n else None


def e2n_cell(pop, si, sd, ni, n, mode, ndraw=NDRAW):
    """-> {eps: {rule: [evaluate recs]}} or None when n cannot be drawn."""
    if n > pop["N"]: return None
    orac = {e: E.oracle_own(pop, e) for e in E.EPSS}; recs = {e: {r: [] for r in E.RULES} for e in E.EPSS}
    for d in range(ndraw):
        rng = np.random.default_rng([20261103, TAG, si, sd, ni, d]); neg = draw_iid(pop, rng, n) if mode == "iid" else draw_field(pop, rng, n)
        if neg is None: return None
        for e in E.EPSS:
            for r, fn in E.RULES.items(): recs[e][r].append(E.evaluate(pop, fn(neg, e), e, orac[e]))
    return recs


def main():
    ap = E.common_args(argparse.ArgumentParser()); ap.add_argument("--ndraw", type=int, default=NDRAW); a = ap.parse_args()
    L, inp = E.open_inputs(a); master = E.load_master(a.master); out = Path(a.out); assert not out.parent.exists(), "output directory exists: use a new directory"
    seeds = [int(x) for x in a.seeds.split(",")]; by, aur, chk, nmin, infeasible, dup = {}, {s: [] for s in E.SITES}, {}, None, [], {}
    for si, site in enumerate(E.SITES):
        for sd in seeds:
            s, y, g = inp.load_site(sd, site); aur[site].append(E.auroc_bases(s, y)); pop = E.make_pop(s, y, g); nmin = pop["N"] if nmin is None else min(nmin, pop["N"]); dup[f"{site}__{sd}"] = float(1 - len(np.unique(pop["sneg"])) / pop["N"])
            for ni, n in enumerate(NS):
                for mode in MODES:
                    recs = e2n_cell(pop, si, sd, ni, n, mode, a.ndraw)
                    if recs is None: infeasible.append((site, mode, n)); continue
                    for e in E.EPSS:
                        for r in E.RULES: chk.setdefault((mode, r, n, e), []).extend(x["viol"] for x in recs[e][r]); by.setdefault((site, mode, n, e), {}).setdefault(r, {})[sd] = recs[e][r]
    res_p = {}
    for (site, mode, n, e), rr in by.items():
        res_p.setdefault(site, {}).setdefault(mode, {}).setdefault(str(n), {})[str(e)] = {r: E.pooled_summary(rr[r]) for r in rr}
    rows = [(r, n, e, float(np.mean(v)), len(v)) for (mode, r, n, e), v in chk.items() if mode == "iid"]
    try: E.check_analytic([(r, n, e, o) for r, n, e, o, _ in rows], min(c for *_, c in rows), 2.0 / nmin)
    except E.CheckError as ex: raise E.CheckError((ex.args[0], {"duplicate_fraction_per_site_seed": dup}))
    table = [dict(rule=r, n=n, eps=e, observed_iid=o, analytic=E.viol_analytic(r, n, e)) for r, n, e, o, _ in rows]
    flags = []
    for site in res_p:
        for n in res_p[site].get("field", {}):
            for e in res_p[site]["field"][n]:
                cf = res_p[site]["field"][n][e].get("CP"); ci = res_p[site].get("iid", {}).get(n, {}).get(e, {}).get("CP")
                if cf and ci and cf["coverage_nonvacuous"] is not None and ci["coverage_nonvacuous"] is not None:   # all-vacuous cells carry no information and are not flagged
                    vf, vi = 1 - cf["coverage_nonvacuous"], 1 - ci["coverage_nonvacuous"]
                    flags.append(dict(site=site, n=int(n), eps=float(e), cp_viol_field_nonvacuous=vf, cp_viol_iid_nonvacuous=vi, vacuous_share_field=cf["vacuous_share"], vacuous_share_iid=ci["vacuous_share"],
                                      E2_b=bool(vf > E.DELTA and vi <= E.DELTA)))
    basis = {site: E.anchor_auroc(master, site, aur[site]) for site in E.SITES} if seeds == E.SEEDS else None
    res = dict(version="e2n_r1", spec="E-batch draft3", ndraw=a.ndraw, seeds=seeds, ns=list(NS), pooled=res_p, infeasible=sorted(set(infeasible)), analytic_check_iid=table, reading_E2_b=flags, anchor_auroc_basis=basis, duplicate_fraction=dup, vacuous_convention=E.VAC_CONVENTION, ci_definition=E.CI_DEFINITION, cp_definition="R1 exact (order statistic), no 200-point grid: differs from the manuscript R1 (grid) and must not be compared with Table 5 R1 rows", inputs_used=dict(inp.used))
    E.write_out(out, res)


if __name__ == "__main__":
    try: main()
    except E.AnchorError as e: print("ANCHOR MISMATCH (exit 3, no output):", e, file=sys.stderr); sys.exit(3)
    except E.CheckError as e: print("ANALYTIC CHECK FAILED (exit 5, nothing is read, no output):", e, file=sys.stderr); sys.exit(5)

# P5 E1 certificate baselines (spec P5_规格_E1E2E3E5E6a_开发站_20261103 draft3 section E1).  LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# NATURE: descriptive comparison of three label-light rules on the three DEVELOPMENT sites (A2 scores, seeds 41..45); no pass/fail on any rule, no constant tuned.  M5, IML, PlasmoCount,
#   the NIH locked test fold, any image, any classifier: never read.  The NIH VALIDATION fold is read only for the anchor below.
# Rules, per draw (m = 22 random fields with all their cells annotated, the pooled negatives n vary):
#   N0  source threshold, no label: tau = ln(t/(1-t)), t = the stored Youden value of the frozen metadata (probability scale).  Only this version (the probability-as-score version of the survey is the erratum).
#   Q   empirical (1-eps) quantile of the pooled negatives: the ceil(n(1-eps))-th smallest score (plug-in, never vacuous).
#   MC  marginal split-conformal: the ceil((n+1)(1-eps))-th smallest; > n -> tau = +inf, counted as vacuous.
#   CP, S3, SPLIT6: NOT recomputed.  S3 (and R1_all when present) coverage is read from the stored REV-1 aggregates only to print the E1-a reading; CP exact for other n is in e2n_resample.py.
# Truth: FPR of tau on ALL negatives of the site (cell weighted); violation = FPR > eps; coverage = 1 - violation rate; sens on all positives; ds_own = sens(tau) - sens(own-oracle tau, c0).
# Stream per (site, seed, eps, draw): default_rng([20261103, 81, site index, seed, eps index, draw]); perm = permutation(fields with >= 1 negative); fields = perm[:22].  NDRAW = 1000 per (site, seed, eps).
# Anchors (all must hold, else exit 3 and no output directory):
#   (M1, r4) per seed, on the NIH VALIDATION fold only (a2ref_s{seed}_eval.json holds just source_threshold {from, rule, value} and nih_validation {auroc, n}; only those two keys are read):
#        (a) n equals the number of validation logits; (b) the AUROC of s = logit1 - logit0 equals nih_validation/auroc within 1e-10; (c) the Youden J threshold recomputed on the probability scale
#        p = sigmoid(s) equals source_threshold/value within 1e-7.  Rule of the stored value: results_src/c1a6/c1a6_infer.py lines 116-118 (youden_threshold: roc_curve(y, p), th[(tpr - fpr).argmax()]) applied
#        to p = softmax(logits)[:, 1] in float32 at lines 181-184 (same rule: code/p5_eval_phase1.py lines 229-231); the float32 softmax is why the tolerance is 1e-7 and not tighter.  The locked fold is not loaded.
#   AUROC of each development site (mean over the five seeds) equals results_master.json survey_arm/<site>|a2|auroc_mean within 1e-12 on one of the two historical score scales (r5; diagnosis of 实验管理:
#   bbbc041_all = AUROC on the float32 probability p32 (323 saturated ties at seed 41), lacuna_B = AUROC on s, mpidb_full = both agree); the matching scale is written to the output (anchor_auroc_basis).
# Reading flags printed (descriptive).  E1-a: Q or MC below 0.90 (non-vacuous draws) on >= 2 sites at one eps, while the two FIELD-AWARE rules S3 and SPLIT6 are >= 0.90 on those sites (S3 from the pinned
#   REV-1 file, SPLIT6 from the pinned master key above; a missing value counts as not satisfied).  R1_all (exact rule that treats the cells as exchangeable) is NOT part of the condition, because under
#   field clustering it under-covers by design (stored coverage 0.57-0.68 on all three sites at the main target); it is reported unchanged next to the flag as "R1_all_stored" (correction made before any run:
#   the earlier condition "S3, R1_all and SPLIT6 >= 0.90" could never hold on the real data).  E1-b: Q and MC >= 0.90 at every site and eps.  E1-c: N0 is reported without choosing among versions.
# Caveat: cov_nv falls back to the all-draw coverage (vacuous draws count as covered, i.e. 1.0) when a rule is vacuous in EVERY draw of a cell, so E1-b can be "satisfied" by an all-vacuous cell; at m = 22 neither
#   Q (never vacuous) nor MC (needs (n+1)(1-eps) <= n, n in the hundreds) is expected to be all-vacuous, and vacuous_share is printed for every cell.
# Usage: python -I -B e1_baselines.py --a2link .. --ids_iter1 .. --mhlh .. --rev1_bbbc .. --rev1_mh .. --rev1_lh .. --mh_rev1_py .. --nihval_dir .. --master paper_en/results_master.json --out <new dir>/e1.json
import argparse, json, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ebatch_lib as E

NDRAW = 1000; STREAM_TAG = 81
K_PIN = "0e6aa9f8f423dbaed286db969afa8a14c9ce477e37bba962732d3b411002147c"   # results_src/m3/k1vb/m3_k1v_precheck.py (NIH validation loader with its own sha checks)
SPLIT6_KEY = "ssrcps/cells|SPLIT6|{site}|{eps}|coverage"                      # pinned master; the family the manuscript uses for SPLIT6 on the development sites


def youden_prob(y, p):
    """Youden J threshold on the PROBABILITY scale, the rule of the stored value: results_src/c1a6/c1a6_infer.py lines 116-118 (and code/p5_eval_phase1.py lines 229-231), roc_curve(y, p) and
    th[(tpr - fpr).argmax()].  Here the same argmax is taken over the full list of distinct score values, sorted descending (first maximiser): sklearn's drop_intermediate only removes collinear
    interior points, whose J lies between two retained corner points, so the first maximiser is the same.  Returns (threshold, J)."""
    p = np.asarray(p, float); y = np.asarray(y); o = np.argsort(-p, kind="stable"); ps, ys = p[o], y[o]
    last = np.r_[np.where(np.diff(ps) != 0)[0], len(ps) - 1]; tps = np.cumsum(ys == 1)[last]; fps = np.cumsum(ys == 0)[last]
    J = tps / (ys == 1).sum() - fps / (ys == 0).sum(); j = int(np.argmax(J)); return float(ps[last[j]]), float(J[j])


def nih_val_anchor(eval_json, s, y, tol_thr=1e-7, tol_auc=1e-10):
    """(M1, r4) Three anchors on the NIH VALIDATION fold only; the json has the keys source_threshold {value} and nih_validation {auroc, n}, and only these two are read.
      (a) n equals the number of validation logits;  (b) AUROC of s = logit1 - logit0 equals nih_validation/auroc within 1e-10;
      (c) the Youden J threshold recomputed on the probability scale p = sigmoid(s) equals source_threshold/value within 1e-7 (the stored value comes from float32 softmax probabilities).
    Any miss -> AnchorError.  Returns (tau_logit, t, info)."""
    d = json.load(open(eval_json)); t = float(d["source_threshold"]["value"]); st = d["nih_validation"]; del d
    s = np.asarray(s, float); y = np.asarray(y)
    if int(st["n"]) != len(s) or len(s) != len(y): raise E.AnchorError(("nih_validation n", int(st["n"]), len(s), len(y)))
    au = E.auroc(s, y)
    if abs(au - float(st["auroc"])) > tol_auc: raise E.AnchorError(("nih_validation auroc", au, float(st["auroc"]), tol_auc))
    tre, J = youden_prob(y, 1.0 / (1.0 + np.exp(-s)))
    if abs(tre - t) > tol_thr: raise E.AnchorError(("Youden threshold", tre, t, tol_thr))
    return math.log(t / (1 - t)), t, dict(n=len(s), auroc=au, auroc_stored=float(st["auroc"]), youden_threshold_recomputed=tre, source_threshold_stored=t, youden_J=J)


def e1_cell(pop, si, sd, ei, eps, tau_n0, ndraw=NDRAW, m=E.M_FIELDS):
    fl = pop["fl"]; assert len(fl) >= m, ("fewer fields with negatives than m", len(fl), m)
    orac = E.oracle_own(pop, eps); recs = {"N0": [], "Q": [], "MC": []}; ns = []; r0 = E.evaluate(pop, tau_n0, eps, orac); taus = {"Q": [], "MC": []}
    for d in range(ndraw):
        rng = np.random.default_rng([20261103, STREAM_TAG, si, sd, ei, d]); F = rng.permutation(len(fl))[:m]
        neg = np.sort(np.concatenate([pop["negf"][fl[i]] for i in F])); ns.append(len(neg))
        tq, tm = E.tau_q(neg, eps), E.tau_mc(neg, eps); taus["Q"].append(tq); taus["MC"].append(tm)
        recs["N0"].append(r0); recs["Q"].append(E.evaluate(pop, tq, eps, orac)); recs["MC"].append(E.evaluate(pop, tm, eps, orac))
    return recs, ns, taus


def cov_nv(c): return c["coverage_nonvacuous"] if c["coverage_nonvacuous"] is not None else c["coverage"]   # a rule that is vacuous in every draw has no non-vacuous coverage; its vacuous_share is read separately


def reading(pooled, stored_cov):
    """pooled[site][eps][rule]; stored_cov[site][eps][rule] from the stored REV-1 aggregates and the pinned master (None when absent).  Coverage = non-vacuous draws only."""
    below = lambda rule, e: [s for s in E.SITES if cov_nv(pooled[s][e][rule]) < 0.90]
    a = []
    for rule in ("Q", "MC"):
        for e in E.EPSS:
            bs = below(rule, e)
            if len(bs) >= 2 and all(all((stored_cov[s][e].get(k) or 0) >= 0.90 for k in ("S3", "SPLIT6")) for s in bs): a.append(dict(rule=rule, eps=e, sites=bs))
    b = all(cov_nv(pooled[s][e][r]) >= 0.90 for s in E.SITES for e in E.EPSS for r in ("Q", "MC"))
    return dict(E1_a=bool(a), E1_a_cells=a, R1_all_stored={s: {str(e): stored_cov[s][e].get("R1_all") for e in E.EPSS} for s in E.SITES}, E1_b=bool(b), E1_c="N0 reported for the logit(t) version only; no version chosen",
                note="E1-a requires S3 and SPLIT6 >= 0.90 on the sites where Q or MC is below 0.90: S3 from the pinned REV-1 file, SPLIT6 from the pinned master key " + SPLIT6_KEY + "; a missing value counts as not satisfied. R1_all (cells treated as exchangeable) is reported unchanged and is not a condition")


def main():
    ap = E.common_args(argparse.ArgumentParser(), with_nihval=True); ap.add_argument("--ndraw", type=int, default=NDRAW); a = ap.parse_args()
    E.refuse_path(a.nihval_dir)
    L, inp = E.open_inputs(a); master = E.load_master(a.master); out = Path(a.out); assert not out.parent.exists(), "output directory exists: use a new directory"
    seeds = [int(x) for x in a.seeds.split(",")]
    kp = E.ROOT / "m3" / "k1vb" / "m3_k1v_precheck.py"; assert E.sha(kp) == K_PIN, "m3_k1v_precheck.py differs from its pin"
    sys.path.insert(0, str(kp.parent)); import m3_k1v_precheck as K
    val = K.load_nihval(a.nihval_dir); taus = {}; m1 = {}
    for sd in seeds:
        tau, t, m1[str(sd)] = nih_val_anchor(Path(a.nihval_dir) / f"a2ref_s{sd}_eval.json", *val[sd]); taus[sd] = (tau, t)
    per, pooled, stored_cov, aur, arrays = {}, {}, {}, {}, {}
    for si, site in enumerate(E.SITES):
        by = {e: {r: {} for r in ("N0", "Q", "MC")} for e in E.EPSS}; aur[site] = []; nmean = {}
        for sd in seeds:
            s, y, g = inp.load_site(sd, site); aur[site].append(E.auroc_bases(s, y)); pop = E.make_pop(s, y, g)
            for ei, eps in enumerate(E.EPSS):
                recs, ns, tt = e1_cell(pop, si, sd, ei, eps, taus[sd][0], a.ndraw)
                nmean.setdefault(eps, []).extend(ns); arrays[f"{site}__{sd}__{eps}__n_neg"] = np.array(ns)
                for r in tt: arrays[f"{site}__{sd}__{eps}__tau_{r}"] = np.array(tt[r], float)
                for r in recs: by[eps][r][sd] = recs[r]
        per[site] = {str(sd): {str(e): {r: E.summarise(by[e][r][sd]) for r in by[e]} for e in E.EPSS} for sd in seeds}
        pooled[site] = {e: {r: E.pooled_summary(by[e][r]) for r in by[e]} for e in E.EPSS}
        for e in E.EPSS: pooled[site][e]["n_neg_mean"] = float(np.mean(nmean[e]))
        agg = json.load(open(inp.rev1[site]))["agg"]
        stored_cov[site] = {e: {r: (agg.get(L.STORED_KEY[site].format(rule=r, eps=e)) or {}).get("cov") for r in ("S3", "R1_all")} for e in E.EPSS}
        for e in E.EPSS: stored_cov[site][e]["SPLIT6"] = master.get(SPLIT6_KEY.format(site=site, eps=e))
    basis = {site: E.anchor_auroc(master, site, aur[site]) for site in E.SITES} if seeds == E.SEEDS else None
    res = dict(version="e1_r1", spec="E-batch draft3", ndraw=a.ndraw, seeds=seeds, eps=list(E.EPSS), delta=E.DELTA, m=E.M_FIELDS, t_stored={str(k): v[1] for k, v in taus.items()},
               tau_logit={str(k): v[0] for k, v in taus.items()}, m1_anchor=m1, anchor_auroc_basis=basis, per_seed=per, pooled={s: {str(e): v for e, v in d.items()} for s, d in pooled.items()},
               stored_coverage_for_reading=stored_cov, reading=reading(pooled, stored_cov), inputs_used=dict(inp.used),
               vacuous_convention=E.VAC_CONVENTION, ci_definition=E.CI_DEFINITION, k_pin=K_PIN, per_draw_thresholds="<stem>_taus.npz: tau_Q, tau_MC (+inf = vacuous) and n_neg per draw, keys site__seed__eps__name")
    E.write_out(out, res, arrays)


if __name__ == "__main__":
    try: main()
    except E.AnchorError as e: print("ANCHOR MISMATCH (exit 3, no output):", e, file=sys.stderr); sys.exit(3)

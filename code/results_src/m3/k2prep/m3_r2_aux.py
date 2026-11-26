# P5 K2 preparation R2 (problem card draft4 fcf90632, Revision Roadmap R2): design-side DESCRIPTION of how well the unlabeled all-cell exceedance a_f(t) tracks the
# negative exceedance p_f(t) across fields, on the three development sites.  LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 approves.
# NATURE: description only; no judgement line, no selection, not confirmatory evidence.  It tells whether K2 (and the K1 control arm) can reduce variance at all.
# DATA: population of each development site and seed (all fields), via m3_replay_lib (score files checked against their manifests; no IML, no PlasmoCount).
#   No draws are replayed (population quantities only), so no self-check against REV-1 is needed here.
# THRESHOLD: per candidate k in {c0, c1s} and eps in {0.10, 0.05}: t = quantile(population negatives of k, 1 - eps'), eps' = stored REV-1 fpr_median of the matching
#   single-candidate rule (B1u_c0 for c0, B1u_c1s for c1s), i.e. the typical realized FPR at which the bound operates.
# PER FIELD f (fields with >= 1 negative, as the replay): N_f negatives, n_f all cells; p_f = #negatives >= t / N_f; a_f = #all cells >= t / n_f (scores of candidate k,
#   no label needed); d_f = p_f - a_f; pi_f = positives / n_f.
# REPORTED per site, k, eps (per seed and mean over the 5 seeds): number of fields; sd(p), sd(a), sd(d), sd(pi) (unweighted and N_f-weighted); Pearson rho(p, a)
#   (unweighted and N_f-weighted); Var(d)/Var(p) (unweighted and weighted); skewness and excess kurtosis of d_f (unweighted, scipy, bias=False);
#   and, for reference, (sd(a) + sd(d))^2 / Var(p): the K2 variance bound relative to the true between-field variance (>= 1 by the triangle inequality).
# Output: one new JSON in a new directory ("x", allow_nan=False; undefined statistics -> None).  Prints only the output path.
# Usage: python -I -B m3_r2_aux.py --a2link .. --ids_iter1 .. --mhlh .. --rev1_bbbc .. --rev1_mh .. --rev1_lh .. --mh_rev1_py .. --out <dir/json>
import argparse, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import kurtosis, skew

sys.path.insert(0, str(Path(__file__).resolve().parent)); import m3_replay_lib as L   # noqa: E402

CANDS = ("c0", "c1s"); RULE_OF = {"c0": "B1u_c0", "c1s": "B1u_c1s"}


def fin(v):
    v = float(v); return v if np.isfinite(v) else None


def wstats(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float); mu = (w*x).sum()/w.sum(); return mu, float(np.sqrt((w*(x-mu)**2).sum()/w.sum()))


def wcorr(x, z, w):
    mx, sx = wstats(x, w); mz, sz = wstats(z, w)
    return None if sx == 0 or sz == 0 else float((w*(np.asarray(x)-mx)*(np.asarray(z)-mz)).sum()/w.sum()/(sx*sz))


def field_stats(pop, k, t):
    sc, y, g = pop["sc"][k], pop["y"], pop["g"]; P, A, D, PI, NN = [], [], [], [], []
    for f in pop["fl"]:
        neg = pop["fl_neg"][k][f]; allf = sc[g == f]; Nf = len(neg); nf = len(allf)
        p = (Nf - np.searchsorted(neg, t, side="left"))/Nf; a = float((allf >= t).mean())
        P.append(p); A.append(a); D.append(p-a); PI.append(float((y[g == f] == 1).mean())); NN.append(Nf)
    P, A, D, PI, NN = map(np.asarray, (P, A, D, PI, NN))
    one = np.ones(len(P)); vp = P.var(); vd = D.var()
    out = dict(n_fields=int(len(P)), sd_p=float(P.std()), sd_a=float(A.std()), sd_d=float(D.std()), sd_pi=float(PI.std()),
               sd_p_w=wstats(P, NN)[1], sd_a_w=wstats(A, NN)[1], sd_d_w=wstats(D, NN)[1],
               rho_pa=wcorr(P, A, one), rho_pa_w=wcorr(P, A, NN),
               var_ratio_d_p=(None if vp == 0 else float(vd/vp)), var_ratio_d_p_w=(None if wstats(P, NN)[1] == 0 else float((wstats(D, NN)[1]/wstats(P, NN)[1])**2)),
               skew_d=fin(skew(D, bias=False)), exkurt_d=fin(kurtosis(D, fisher=True, bias=False)),
               k2_bound_over_var_p=(None if vp == 0 else float((A.std()+D.std())**2/vp)))
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); out = Path(a.out); assert not out.exists(), "output exists"
    inp = L.Inputs(a.a2link, a.ids_iter1, a.mhlh, a.rev1_bbbc, a.rev1_mh, a.rev1_lh, a.mh_rev1_py)
    cells = []
    for site in L.DEV:
        agg = json.load(open(inp.rev1[site]))["agg"]
        per = {(k, e): [] for k in CANDS for e, _ in L.SETTINGS}
        for sd in L.SEEDS:
            pop = L.population(inp, sd, site)
            for k in CANDS:
                for e, _ in L.SETTINGS:
                    ep = agg[L.STORED_KEY[site].format(rule=RULE_OF[k], eps=e)]["fpr_median"]; t = float(np.quantile(pop["sneg"][k], 1 - ep))
                    per[(k, e)].append(dict(seed=sd, eps_prime=ep, threshold=t, **field_stats(pop, k, t)))
        for (k, e), P in per.items():
            keys = [kk for kk in P[0] if kk not in ("seed",)]
            mean = {kk: (None if any(p[kk] is None for p in P) else float(np.mean([p[kk] for p in P]))) for kk in keys}
            cells.append(dict(site=site, cand=k, eps=e, mean_over_seeds=mean, per_seed=P))
    rep = dict(script_sha256=L.sha(__file__), lib_sha256=L.sha(L.__file__), inputs_sha256=dict(L.PIN), score_files_sha256=inp.used,
               nature="design-side description for K2 (problem card draft4 R2); no judgement, not confirmatory evidence", cells=cells)
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False)
    out.parent.mkdir(parents=False, exist_ok=False)
    with open(out, "x") as fh: fh.write(txt)
    print("written", out)


if __name__ == "__main__":
    main()

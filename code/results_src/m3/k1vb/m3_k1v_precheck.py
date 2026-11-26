# P5 K1-v design-side pre-check (problem card draft6 f859dd4e section 3 K1-v; addendum A c7571fdf: sign e + v, NIH validation files, Platt).  LOCAL, read-only.
# WRITTEN, NOT RUN, until 实验管理 approves; runs only after the K1vB evaluation has finished (one heavy job at a time).
# NATURE: description only, no judgement line.  Pre-stated condition for a next-round problem card (addendum A1): on old MP-IDB, Var(e + v)/Var(e) < 0.5.
# DATA:
#   development sites (bbbc041_all, mpidb_full rebuilt, lacuna_B rebuilt, old MP-IDB) via m3_stress.StressInputs (manifest / SHA256SUMS checked; old site one key);
#   NIH VALIDATION fold only, for the Platt calibration of q: a2ref_s{41..45}_nih_validation_logits.npz (ONLY key "logits") and ids_nih_validation.npz (ONLY key "y"),
#   sha pinned (filled by program).  No nih_locked_test key, no name/field key, no iteration over keys.  No M5, IML, PlasmoCount.
# CALIBRATION (addendum A3, fixed before the pre-check): per seed, q(s) = sigmoid(a s + b) on the raw score s = logit1 - logit0, maximum likelihood (scipy minimize,
#   no penalty); parameters frozen and reported.
# QUANTITIES per site, seed, candidate k in {c0, c1s}, eps in {0.10, 0.05}: t = quantile(population negatives of k, 1 - eps'), eps' = stored REV-1 fpr_median of the
#   single-candidate rule (B1u_c0 / B1u_c1s) (as R2).  Per field f of F: z_f = #neg >= t - eps n_f; w_f = #all >= t - eps c_f (scores of k); e_f = z_f - w_f;
#   v_f = sum_cells q(s_raw) (1{s_k >= t} - eps).  Reported: Var(e)/n, Var(e + v)/Var(e) (unweighted over fields), corr(e, -v), and the identity check
#   max_f |e_f + P_f (TPR_f(t) - eps)| (must be ~ 0; positives counted from labels, P_f = 0 fields give 0).
# Output: one new JSON in a new directory ("x", allow_nan=False; undefined -> None).  Prints only the output path.
# Usage: python -I -B m3_k1v_precheck.py <the m3_stress inputs: --a2link .. --ids_iter1 .. --mhlh .. --rev1_bbbc .. --rev1_mh .. --rev1_lh .. --mh_rev1_py ..
#        --old_s41 .. --old_dir ..> --nihval_dir <runs/c1a6_stage2_20261027T133543Z/a2ref> --out <dir/json>
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, log_expit

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "k2prep")); sys.path.insert(0, str(ROOT / "stress"))
import m3_replay_lib as L   # noqa: E402
import m3_stress as ST   # noqa: E402

NIHVAL_PIN = {41: "8e7fcd34d0dc119179dab406898ca9d08604e866885b6da83f5c262923b6cbd8", 42: "6d85e599802609f7a33dedc471a53f58d52bf4d1dbe52fb9b3ca94fa9fbab79f", 43: "71cd404368a600ad7a3a5470da93c0a7a83791e3173bdca14da8d06e6f3d886b", 44: "4b046baf7967608712d98f7fd30d3cdc1f6a0926c4014d6d8febdb43eaaa427f", 45: "7d71e4e47451fe32eb54530470176ddd05194330ac824df2eff458c77edaab66", "ids": "7780177f854dc31b0122d4bba45dbf79c2c735d5f243e54f3954abf457039d62"}
OLD = ST.OLD; SITES = ("bbbc041_all", "mpidb_full", "lacuna_B", OLD); CANDS = ("c0", "c1s"); RULE_OF = {"c0": "B1u_c0", "c1s": "B1u_c1s"}


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def fin(v):
    v = float(v); return v if np.isfinite(v) else None


def platt(s, y):
    s = np.asarray(s, float); y = np.asarray(y, float)
    nll = lambda p: -(y*log_expit(p[0]*s + p[1]) + (1 - y)*log_expit(-(p[0]*s + p[1]))).sum()
    r = minimize(nll, x0=np.array([1.0, 0.0]), method="BFGS"); assert r.success, ("platt fit", r.message)
    return float(r.x[0]), float(r.x[1])


def load_nihval(d):
    d = Path(d); assert "iml" not in str(d).lower() and "plasmo" not in str(d).lower(), "confirmatory-site path refused"
    fi = d / "ids_nih_validation.npz"; assert sha(fi) == NIHVAL_PIN["ids"], ("sha", fi); y = np.asarray(np.load(fi)["y"])   # key y only
    out = {}
    for sd in L.SEEDS:
        f = d / f"a2ref_s{sd}_nih_validation_logits.npz"; assert sha(f) == NIHVAL_PIN[sd], ("sha", f)
        z = np.load(f)["logits"]; assert len(z) == len(y), ("length", f); out[sd] = ((z[:, 1] - z[:, 0]).astype(float), y)        # key logits only
    return out


def field_stats(pop, k, t, eps, q):
    """q: per-cell calibrated positive probability on the raw score (same order as pop['s'])."""
    sk, g, y = pop["sc"][k], pop["g"], pop["y"]; E_, V_, ID = [], [], []
    for f in pop["fl"]:
        mk = g == f; neg = sk[mk & (y == 0)]; allf = sk[mk]; pos = sk[mk & (y == 1)]
        z = (neg >= t).sum() - eps*len(neg); w = (allf >= t).sum() - eps*len(allf); e = z - w
        v = float((q[mk]*((allf >= t).astype(float) - eps)).sum())
        P = len(pos); tpr = float((pos >= t).mean()) if P else 0.0
        E_.append(e); V_.append(v); ID.append(abs(e + P*(tpr - eps)))
    E_, V_ = np.asarray(E_, float), np.asarray(V_, float); ve = E_.var()
    c = np.corrcoef(E_, -V_)[0, 1] if E_.std() > 0 and V_.std() > 0 else np.nan
    return dict(n_fields=int(len(E_)), var_e=float(ve), var_ratio_e_plus_v=(None if ve == 0 else float((E_ + V_).var()/ve)), corr_e_minus_v=fin(c),
                identity_max_abs=float(max(ID)))


def main():
    ap = argparse.ArgumentParser()
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "old_s41", "old_dir", "nihval_dir", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); out = Path(a.out); assert not out.exists(), "output exists"
    inp = ST.StressInputs(a); assert sorted(inp.old) == [41, 42, 43, 44, 45], "old MP-IDB needs all 5 seeds"
    nv = load_nihval(a.nihval_dir); params = {sd: platt(*nv[sd]) for sd in L.SEEDS}
    cells = []
    for site in SITES:
        agg = json.load(open(inp.rev1[site]))["agg"]; per = {(k, e): [] for k in CANDS for e, _ in L.SETTINGS}
        for sd in L.SEEDS:
            pop = L.population(inp, sd, site); aa, bb = params[sd]; q = expit(aa*pop["s"] + bb)
            for k in CANDS:
                for e, _ in L.SETTINGS:
                    ep = agg[ST.PATTERN[site].format(m=22, rule=RULE_OF[k], eps=e)]["fpr_median"]; t = float(np.quantile(pop["sneg"][k], 1 - ep))
                    per[(k, e)].append(dict(seed=sd, eps_prime=ep, threshold=t, **field_stats(pop, k, t, e, q)))
        for (k, e), P in per.items():
            r = [p["var_ratio_e_plus_v"] for p in P]
            cells.append(dict(site=site, cand=k, eps=e, var_ratio_e_plus_v_mean=(None if any(x is None for x in r) else float(np.mean(r))),
                              identity_max_abs=float(max(p["identity_max_abs"] for p in P)), per_seed=P))
    old_r = [c["var_ratio_e_plus_v_mean"] for c in cells if c["site"] == OLD]
    rep = dict(script_sha256=sha(__file__), stress_sha256=sha(ST.__file__), lib_sha256=sha(L.__file__), inputs_sha256=dict(L.PIN), score_files_sha256=inp.used,
               nihval_sha256={str(k): v for k, v in NIHVAL_PIN.items()}, platt={str(sd): dict(a=p[0], b=p[1]) for sd, p in params.items()},
               nature="design-side description (draft6 K1-v pre-check, addendum A); no judgement",
               next_round_condition="old MP-IDB Var(e + v)/Var(e) < 0.5 (addendum A1)",
               old_mpidb_ratios=old_r, old_mpidb_condition_met=bool(old_r and all(x is not None and x < 0.5 for x in old_r)), cells=cells)
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False)
    out.parent.mkdir(parents=False, exist_ok=False)
    with open(out, "x") as fh: fh.write(txt)
    print("written", out)


if __name__ == "__main__":
    main()

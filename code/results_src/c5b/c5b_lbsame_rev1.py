# C5b EVAL (PREREG_P5_C5B_FRAMING_AUG_20261026.md + errata 1-3): derived from lb_same_paired_rev1.py (sha256 9641fa7e6adc...) by build_c5b_eval.py. ARMS = ['a2'] keeps ai = 0 (the frozen A2's random streams); the arm label LAB is read from env C5B_ARM (c5b-f / c5b-fs), logits from env C5B_LOGDIR/<LAB>_s<seed>_<site>_logits.npz; output file names carry LAB; see the build script for the full list of edits.
# SAMEPOP (PREREG_P5_SAMEPOP_DECOMP_20261026.md): rev1 with the paired stream (key8 = 8, sidx = 2) on the lb_same subset (old de-duplicated Lacuna-B restricted to the cells that the rebuilt set holds as positives); derived from lb_dedup_paired_rev1.py; only the site / data / output names differ.
# PAIRED (PREREG_P5_PAIRED_STREAM_SUPPLEMENT_20261026.md, sha256 2881a982..., erratum bea6b88b...): restores the ORIGINAL random stream of Lacuna-B in lb_dedup_rev1.py: seed [20261022, 8, ai, sd, m, 2] (key8 = 8, sidx = 2); per-draw output added (plus one extra list in rec()); nothing else differs from lb_dedup_rev1.py.
# LB_DEDUP adaptation (old Lacuna-B, positives de-duplicated, site name lb_dedup, site index 33) per PREREG_P5_LACUNAB_OLD_DEDUP_20261025.md, derived from the REVIEWED mh_rev1.py by build_lb_dedup.py: only data loading, site name (mh -> lb_dedup), site index (31 -> 33) and output names changed; the mh header and original docstring follow.
# (mh header) MH adaptation (rebuilt MP-IDB host-cell set, site index 31) per PREREG_P5_MPIDB_HOSTCELL_REEVAL_20261025.md section 3: only paths, manifest column mapping, site name and output names changed; original docstring follows.
"""R8 IML per PREREG_R8_IML_20261023.md; adapted from rev1.py (only site list and data loading changed; REV-1 definitions, 30a81969...). Post hoc re-analysis: new reference (oracle at FPR=eps), realized-FPR stats, R1-field and R1-all added. CPU only."""
import json, sys, numpy as np
from pathlib import Path
from multiprocessing import Pool
from scipy.special import betainc
from scipy.stats import beta as betad, f as fdist
sys.path.insert(0, str(Path(__file__).parent)); import s8_bounds as B8
INF = Path(sys.argv[1]); E7 = Path(sys.argv[2]); E8 = Path(sys.argv[3]); OUT = Path(sys.argv[4]); OUT.mkdir(parents=True, exist_ok=True)
import os as _os; LAB = _os.environ["C5B_ARM"]; LOGDIR = Path(_os.environ["C5B_LOGDIR"]); assert LAB in ("a2", "c5b-f", "c5b-fs"), LAB   # C5b: arm label of the files; ARMS keeps ONE entry so that ai = 0 (the A2 streams)
ARMS = ["a2"]; SEEDS = [41,42,43,44,45]; NDRAW = 200; SETTINGS = [(0.10, 0.10), (0.05, 0.10)]; G = 200
SITES = ["lb_same"]; CEN = ["c0", "c1s", "c2s"]; MS_DEV = [3, 5, 10, 22, 30, 50]
def b1u(x, n, delta):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X/N
    n0 = (N-(n**2).sum()/N)/(m-1); pi = x/n; msb = (n*(pi-ph)**2).sum()/(m-1); msw = (n*pi*(1-pi)).sum()/(N-m)
    if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0
    else:
        F = msb/msw; flo = fdist.ppf(delta, m-1, N-m); th = (F/flo-1)/n0; iccu = 0.0 if th <= 0 else th/(1+th)
    nbw = (n**2).sum()/N; neff = N/(1+(nbw-1)*iccu); xeff = X*neff/N
    return float(betad.ppf(1-delta, xeff+1, neff-xeff))
def R1(n, eps, delta):
    for k in range(1, int(n)+1):
        if betainc(n-k+1, k, eps) >= 1-delta: return k, True
    return int(n), False
def load_site(arm, sd, site):
    if site == "lb_same":
        lid = np.load(E8/"ids_lb_same.npz", allow_pickle=True); z = np.load(LOGDIR/f"{LAB}_s{sd}_lb_same_logits.npz")["logits"]; return (z[:, 1]-z[:, 0]).astype(float), lid["y"], lid["field"]
    if site == "lacuna_B":
        lid = np.load(E7/"ids_lacuna_B.npz", allow_pickle=True); z = np.load(E7/f"{arm}_s{sd}_lacuna_B_logits.npz")["logits"]; return (z[:, 1]-z[:, 0]).astype(float), lid["y"], lid["field"]
    ids = np.load(INF/"ids_iter1.npz", allow_pickle=True); z = np.load(INF/f"{arm}_s{sd}_scores.npz")[site+"_logits"]; return (z[:, 1]-z[:, 0]).astype(float), ids[site+"_y"], ids[site+"_group"]
def cstat(v, kind, gm):
    if kind == "c0": return 0.0
    w = len(v)/(len(v)+20.0)
    if kind == "c1s": return float(gm+w*(np.median(v)-gm))
    return float(gm+w*(np.median(np.sort(v)[:int(np.ceil(0.7*len(v)))])-gm))
def task(a):
    ai, sd, si = a; arm = ARMS[ai]; site = SITES[si]; key8 = 8; sidx = 2
    s, y, g = load_site(arm, sd, site); fl = list(np.unique(g)); gm = float(np.median(s)); cs = {f: {k: cstat(s[g == f], k, gm) for k in CEN} for f in fl}
    sc = {k: s-np.array([cs[f][k] for f in g]) for k in CEN}
    fl_neg = {k: {f: np.sort(sc[k][(g == f) & (y == 0)]) for f in fl} for k in CEN}; fl_posv = {k: {f: sc[k][(g == f) & (y == 1)] for f in fl} for k in CEN}
    fl = [f for f in fl if len(fl_neg["c0"][f]) > 0]
    sneg = {k: np.sort(sc[k][y == 0]) for k in CEN}; spos = {k: np.sort(sc[k][y == 1]) for k in CEN}; N = len(sneg["c0"])
    sens_at = lambda k, t: 1-np.searchsorted(spos[k], t, side="left")/len(spos[k]); fpr_at = lambda k, t: 1-np.searchsorted(sneg[k], t, side="left")/N
    orac_unc = {e: sens_at("c0", np.quantile(sneg["c0"], 1-e)) for e in (0.10, 0.05)}; orac_own = {(k, e): sens_at(k, np.quantile(sneg[k], 1-e)) for k in CEN for e in (0.10, 0.05)}; old05 = sens_at("c0", np.quantile(sneg["c0"], 0.95))
    # Proposition 1 conservativeness at the oracle_0.10 threshold (population level, uncentered)
    tau = np.quantile(sneg["c0"], 0.90); pbar = float(np.mean([ (fl_neg["c0"][f] >= tau).mean() for f in fl])); q = float(np.mean([ (fl_neg["c0"][f].max() >= tau) for f in fl]))
    res = {}
    def rec(name, m, eps, k, t, cert, clicks, nlab):
        e = res.setdefault((arm, site, m, name, eps), {"fpr": [], "sens": [], "dse": [], "dso": [], "dold": [], "cert": [], "clicks": [], "nlab": [], "kk": []}); sn = sens_at(k, t)
        e["fpr"].append(fpr_at(k, t)); e["sens"].append(sn); e["dse"].append(sn-orac_unc[eps]); e["dso"].append(sn-orac_own[(k, eps)]); e["dold"].append(sn-old05); e["cert"].append(cert); e["clicks"].append(clicks); e["nlab"].append(nlab); e["kk"].append(k)
    for m in [22]:
        if m > len(fl): continue
        rng = np.random.default_rng([20261022, key8, ai, sd, m, sidx])
        for _ in range(NDRAW):
            sel_ = rng.choice(len(fl), m, replace=False); F = [fl[i] for i in sel_]; clicks = sum(len(fl_posv["c0"][f]) for f in F); ni = np.array([len(fl_neg["c0"][f]) for f in F]); subs = {}
            for eps, dl in SETTINGS:
                r_ = np.random.default_rng(int(rng.integers(1 << 30))); subs[eps] = [r_.choice(n_, 2, replace=False) if n_ >= 2 else np.arange(n_) for n_ in ni]
            D = {}
            for k in CEN:
                negs = [fl_neg[k][f] for f in F]; pooled = np.concatenate(negs); grid = np.unique(np.quantile(pooled, 1-np.arange(G)/G))[::-1]
                X = np.stack([ni[i]-np.searchsorted(negs[i], grid, side="left") for i in range(m)]).astype(float)
                posl = np.concatenate([fl_posv[k][f] for f in F]) if clicks > 0 else np.array([]); D[k] = (negs, pooled, grid, X, posl)
            cache = {}
            def ub(rule, k, j):
                key = (rule, k, j)
                if key not in cache:
                    x = D[k][3][:, j]
                    cache[key] = B8.b1(x, ni, 0.10)[0] if rule == "B1" else b1u(x, ni, 0.10) if rule == "U1" else b1u(x, ni, 0.10/3)
                return cache[key]
            def tau_of(rule, k, eps):
                grid = D[k][2]; j = 0
                try:
                    while j < len(grid) and ub(rule, k, j) <= eps: j += 1
                except Exception: j = 0
                return (grid[j-1], True) if j > 0 else (grid[0], False)
            for eps, dl in SETTINGS:
                T = {("B1", "c0"): tau_of("B1", "c0", eps)}
                for k in CEN: T[("U1", k)] = tau_of("U1", k, eps); T[("U3", k)] = tau_of("U3", k, eps)
                nlab_all = len(D["c0"][1])
                k2, ok2 = R1(2*m, eps, dl); negs0 = D["c0"][0]; two = np.sort(np.concatenate([negs0[i][subs[eps][i]] for i in range(m)])); rec("R1_2", m, eps, "c0", two[min(k2, len(two))-1], ok2, clicks, 2*m)
                mx = np.sort([a_.max() for a_ in negs0]); kf, okf = R1(m, eps, dl); rec("R1_field", m, eps, "c0", mx[kf-1], okf, clicks, nlab_all)
                pl = np.sort(D["c0"][1]); ka, oka = R1(len(pl), eps, dl); rec("R1_all", m, eps, "c0", pl[ka-1], oka, clicks, nlab_all)
                t_, c_ = T[("B1", "c0")]; rec("B1", m, eps, "c0", t_, c_, clicks, nlab_all)
                for k in CEN: t_, c_ = T[("U1", k)]; rec("B1u_"+k, m, eps, k, t_, c_, clicks, nlab_all)
                cert_c = [k for k in CEN if T[("U3", k)][1]]
                if not cert_c: rec("S3", m, eps, "c0", T[("U3", "c0")][0], False, clicks, nlab_all)
                else:
                    if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
                    else:
                        sl = [(D[k][4] >= T[("U3", k)][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
                    rec("S3", m, eps, ch, T[("U3", ch)][0], True, clicks, nlab_all)
    print("done", arm, sd, site, flush=True)
    return res, {"pbar": pbar, "q": q}
if __name__ == "__main__":
    tasks = [(ai, sd, si) for ai in range(1) for sd in SEEDS for si in range(len(SITES))]
    with Pool(4) as p: parts = p.map(task, tasks, chunksize=1)
    merged = {}; prop1 = {}
    for (ai, sd, si), (r, pr) in zip(tasks, parts):
        prop1.setdefault(f"{ARMS[ai]}|{SITES[si]}", []).append(pr)
        for k, e in r.items():
            d = merged.setdefault(k, {kk: [] for kk in e})
            for kk in d: d[kk] += e[kk]
    agg = {}
    for (arm, site, m, name, eps), e in merged.items():
        fpr = np.array(e["fpr"])
        agg[f"{arm}|{site}|{m}|{name}|{eps}"] = {"cov": float((fpr <= eps).mean()), "cov_lb010": float((np.maximum(0, fpr-0.10) <= eps).mean()), "cannot_certify": float(1-np.mean(e["cert"])), "ds_eps": float(np.mean(e["dse"])), "ds_own": float(np.mean(e["dso"])), "ds_old05": float(np.mean(e["dold"])),
            "fpr_mean": float(fpr.mean()), "fpr_median": float(np.median(fpr)), "fpr_p90": float(np.quantile(fpr, 0.9)), "sens": float(np.mean(e["sens"])), "clicks": float(np.mean(e["clicks"])), "n_labeled_neg": float(np.mean(e["nlab"])), "n": len(fpr)}
    p1 = {k: {"pbar": float(np.mean([x["pbar"] for x in v])), "q": float(np.mean([x["q"] for x in v])), "ratio_q_over_pbar": float(np.mean([x["q"]/x["pbar"] for x in v if x["pbar"] > 0]))} for k, v in prop1.items()}
    json.dump({"agg": agg, "prop1": p1}, open(OUT/("c5b_lbsame_rev1_results_%s.json" % LAB), "w"), indent=1); print("saved all")
    # ---- PAIRED per-draw output (supplement section 4): writes the per-draw lists that rec() has already accumulated (cert, fpr, sens, clicks, selected candidate); nothing is computed here; one write after the aggregation and the summary JSON
    _stem = Path(__file__).stem
    for _arm in ARMS:
        _o = {}
        for _key, _e in merged.items():
            if _key[0] != _arm or _key[2] != 22: continue
            _tag = "|".join(str(_x) for _x in _key[3:])
            for _f in ("cert", "fpr", "sens", "clicks"): _o[_tag + "|" + _f] = np.array(_e[_f])
            if "kk" in _e and _tag.startswith("S3"): _o[_tag + "|sel"] = np.array([_k if _c else "none" for _k, _c in zip(_e["kk"], _e["cert"])])
            if "sel" in _e: _o[_tag + "|sel"] = np.array([CEN[_i] if _i >= 0 else "none" for _i in _e["sel"]])
        np.savez(OUT / (_stem + "_perdraw_" + LAB + "_m22.npz"), **_o)

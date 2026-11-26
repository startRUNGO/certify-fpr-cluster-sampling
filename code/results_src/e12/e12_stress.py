# P5 E12: synthetic stress tests of S3, S3c, SPLIT6 (spec P5_规格_E11E12_水平修正与压力测试_开发站与合成_20261108_draft2.md, frozen in PRECHECK_LEDGER).  LOCAL, no project data.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# NATURE: a supplementary stress test added after the internal blind review; no rule and no constant is changed.  Every score is generated here by a copy of the E10 generator extended with four stress families (make_site2, identical to
#   e10_synth.make_site when every stress parameter is neutral).  Real data, M5, IML, PlasmoCount, NIH, images, classifiers: never read.  The only file read besides the frozen libraries is the E10 result JSON (anchor).
# Cells (24, 5 realisations each; seeds 41..45 are realisation ids; base = tau 0.5, r 100, pi 0.05, normal tail, F 300, SEP 2.5):
#   0 base (E10 streams, so its S3 / SPLIT6 coverage must equal E10 cell 0, anchor d);
#   A rare high-exceedance fields: rho in {0.01, 0.03, 0.10} x delta in {2, 4}: round(rho F) fields (at least one) get +delta on all their NEGATIVE scores (6 cells);
#   B heavy-tailed field sizes: n_f = max(5, round(r exp(sigma Z - sigma^2/2))), sigma in {0.5, 1.0, 1.5, 2.0} (4 cells);
#   C slide -> field -> cell nesting (CONTROL: fields still drawn uniformly over the whole site): S in {3, 6, 12} slides, omega in {0.5, 0.8}: mu_f = a_s + b_f, a ~ N(0, omega tau^2), b ~ N(0, (1 - omega) tau^2) (6 cells);
#   C2 two-stage sampling (STRESS): S = 12, omega in {0.5, 0.8}; first s in {2, 4} slides without replacement, then 22 fields without replacement among the fields of those slides; truth is the whole site (4 cells);
#   D small F: F in {60, 100, 150} (3 cells).
# Arms per cell (same draws): S3, S3c (every candidate b1u at 0.10/6), SPLIT6, B1u_c0 / B1u_c1s from the frozen rules.  m = 22, 200 draws per realisation, eps in {0.10, 0.05}, delta = 0.10.
# The field sampler is a hook of the copy of m3_split_dev_eval.run_site_seed below; with the uniform sampler the copy equals the library row by row (anchor c); SPLIT6 only (K1, K1vB, SPLIT4/8 are not used here and use separate streams, so S3/SPLIT6 rows are unchanged).
# Gate (exit 4, before anything is read): any argument containing plasmo / iml / r8 / m5 (any case), or an already-imported module named m5_* / e8_* / r8_* / e7_*.
# Anchors (exit 3, no output directory): pins of the frozen libraries and e10_synth; (a) generator checks of E10 (ICC at tau 0 and 1, iid type-I of B1u); (b) make_site2 with neutral parameters equals e10_synth.make_site bit for bit;
#   (c) uniform sampler: copy of run_site_seed equals SP.run_site_seed on a base-cell realisation (5 draws) for S3, B1u_c0, B1u_c1s, SPLIT6; (d) base cell S3 and SPLIT6 coverage per realisation equal E10 cell 0 within 1e-12 (only when
#   --reps 5 --ndraw 200 and cell 0 is run; otherwise recorded as not_applicable); (e) per draw: S3c certified => S3 certified.
# Reading data per cell, rule, eps: coverage (mean over realisations), per-realisation coverage, seed-cluster 90% interval (ebatch_lib.seed_cluster_ci, as E10), band, near_boundary, cannot_certify, ds_own; paired S3c - S3; realised ICC per realisation;
#   cells A: share of draws containing a rare field and S3/S3c coverage with / without one.  band: >= 0.90 hold_1-delta; 0.80..0.90 below_1-delta_within_1-2delta; < 0.80 fail; near_boundary: within 0.01 of 0.90 or 0.80.  No cell is excluded.
# Output: one NEW directory: e12_summary.json, e12_rows.json (per draw, without ni and fields), SHA256SUMS ('x' mode, allow_nan=False).  Prints the directory only.
# Usage: python -I -B e12_stress.py --e10_json <ebatch_e10.json> --out <new dir> [--cells 0,5] [--reps 5] [--ndraw 200] [--check_draws 2000]
import os
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ.setdefault(_k, "8")
import argparse, hashlib, json, re, sys, time
from pathlib import Path
import numpy as np
from scipy.special import expit

HERE = Path(__file__).resolve().parent; R = HERE.parent
FORBID_ARG = re.compile(r"plasmo|iml|r8|m5", re.I); FORBID_MOD = re.compile(r"^(m5_|e8_|r8_|e7_)")
PINS = {"m3_replay_lib": "ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de", "m3_split_dev_eval": "cbf1e7f4f33dd1c74cbfd0f3ac4d34a397e8a7ef44b82c9365ca117733bf6457",
        "s8_bounds": "54a2af8d3704c16b0994dcdc584d4f407cd9564cc07d49ad74cdf7c4cad5ae7f", "e10_synth": "00be97b1f15a0fd5ac7c03b06d41f6adc0ac44a2c11b165e54748fbe68f7e08e",
        "ebatch_lib": "359b6bce2dfdbf50a2acb8b909d631ab40d445447275915e2318e70981ef007f"}
EPSS = (0.10, 0.05); M = 22; LEV = {"S3": 0.10 / 3, "S3c": 0.10 / 6}


class AnchorError(Exception): pass


def gate(argv):
    bad = [a for a in argv if FORBID_ARG.search(a)]
    if bad: print("GATE (exit 4): confirmatory-site string in argument:", bad[:3], file=sys.stderr); sys.exit(4)
    mods = [m for m in sys.modules if FORBID_MOD.match(m.split(".")[-1])]
    if mods: print("GATE (exit 4): confirmatory-site module loaded:", mods[:3], file=sys.stderr); sys.exit(4)


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def cells():
    base = dict(name="base", tau=0.5, r=100, pi=0.05, tail="normal", F=300, size_sigma=0.0, rare=None, slides=None, two_stage=None); out = [base]
    for rho in (0.01, 0.03, 0.10):
        for dl in (2.0, 4.0): out.append(dict(base, name=f"A_rho{rho}_d{dl}", rare=(rho, dl)))
    for sg in (0.5, 1.0, 1.5, 2.0): out.append(dict(base, name=f"B_sig{sg}", size_sigma=sg))
    for S in (3, 6, 12):
        for om in (0.5, 0.8): out.append(dict(base, name=f"C_S{S}_om{om}", slides=(S, om)))
    for om in (0.5, 0.8):
        for s in (2, 4): out.append(dict(base, name=f"C2_S12_om{om}_s{s}", slides=(12, om), two_stage=s))
    for F in (60, 100, 150): out.append(dict(base, name=f"D_F{F}", F=F))
    return out


def stream(ci, sd, purpose):
    """Cell 0 uses the E10 streams (so that it reproduces E10 cell 0); the stress cells use their own."""
    return [20261103, 10, 0, sd, purpose] if ci == 0 else [20261108, 12, ci, sd, purpose]


def make_site2(E10, rng, tau, r, pi, tail, F, size_sigma=0.0, rare=None, slides=None, sep=None):
    """Scores, labels, field ids, slide-of-field (or None) and the rare field ids.  With size_sigma = 0, rare = None, slides = None this is e10_synth.make_site, bit for bit."""
    sep = E10.SEP if sep is None else sep
    if size_sigma > 0:
        z = rng.normal(size=F); nf = np.maximum(5, np.round(r * np.exp(size_sigma * z - size_sigma ** 2 / 2)).astype(int))
    else: nf = np.maximum(5, rng.poisson(r, F))
    npos = rng.poisson(pi / (1 - pi) * nf); slide_of = None
    if slides is None: mu = rng.normal(0, tau, F) if tau > 0 else np.zeros(F)
    else:
        S, om = slides; a = rng.normal(0, np.sqrt(om) * tau, S); b = rng.normal(0, np.sqrt(1 - om) * tau, F)
        slide_of = np.minimum(np.arange(F) // (F // S), S - 1); mu = a[slide_of] + b
    shift = np.zeros(F); rare_fields = np.array([], int)
    if rare is not None:
        rho, dl = rare; k = max(1, int(round(rho * F))); rare_fields = np.sort(rng.choice(F, k, replace=False)); shift[rare_fields] = dl
    s, y, g = [], [], []
    for f in range(F):
        s.append(mu[f] + shift[f] + E10.noise(rng, nf[f], tail)); y.append(np.zeros(nf[f], int)); g.append(np.full(nf[f], f))
        if npos[f]: s.append(mu[f] + sep + E10.noise(rng, npos[f], tail)); y.append(np.ones(npos[f], int)); g.append(np.full(npos[f], f))
    return np.concatenate(s), np.concatenate(y), np.concatenate(g), slide_of, rare_fields


def uniform_sampler(rng, fl, m): return rng.choice(len(fl), m, replace=False)


def two_stage_sampler(slide_of, S, s):
    def smp(rng, fl, m):
        chosen = set(int(x) for x in rng.choice(S, s, replace=False)); cand = [i for i, f in enumerate(fl) if int(slide_of[int(f)]) in chosen]
        pick = rng.choice(len(cand), m, replace=False); return [cand[i] for i in pick]
    return smp


def run_site_seed_e12(SP, inp, sd, site, platt, sampler):
    """Copy of m3_split_dev_eval.run_site_seed restricted to S3, B1u_c0/c1s and SPLIT6 (ms = 6), with a field-sampler hook and the added arm S3c.  K1, K1vB, SPLIT4/8 are dropped: they use their own streams (bootstrap weights, perm) and do not
    consume the main stream, so the kept rows equal the library rows (anchor c)."""
    L, E, VB = SP.L, SP.E, SP.VB; A3 = SP.A3
    pop = L.population(inp, sd, site); fl, fl_neg, fl_posv = pop["fl"], pop["fl_neg"], pop["fl_posv"]; sens_at, fpr_at = pop["sens_at"], pop["fpr_at"]
    AC = {k: E.AllCells(pop, k) for k in L.CEN}; q = expit(platt[0] * pop["s"] + platt[1]); QC = {k: SP.QCells(pop["sc"][k], q, pop["g"], fl) for k in L.CEN}
    key8, sidx = L.SITECFG[site]; ai = L.AI; rows = []; m = L.M; si = L.DEV.index(site)
    rng = np.random.default_rng([20261022, key8, ai, sd, m, sidx])
    for d_ in range(L.NDRAW):
        sel_ = sampler(rng, fl, m); F = [fl[i] for i in sel_]; clicks = sum(len(fl_posv["c0"][f]) for f in F); ni = np.array([len(fl_neg["c0"][f]) for f in F])
        for eps, dl in L.SETTINGS:   # kept: consumes the stream exactly as REV-1 (R1-2 subsets)
            r_ = np.random.default_rng(int(rng.integers(1 << 30))); [r_.choice(n_, 2, replace=False) if n_ >= 2 else np.arange(n_) for n_ in ni]
        D = {}
        for k in L.CEN:
            negs = [fl_neg[k][f] for f in F]; pooled = np.concatenate(negs); grid = np.unique(np.quantile(pooled, 1 - np.arange(L.G) / L.G))[::-1]
            X = np.stack([ni[i] - np.searchsorted(negs[i], grid, side="left") for i in range(m)]).astype(float)
            posl = np.concatenate([fl_posv[k][f] for f in F]) if clicks > 0 else np.array([]); D[k] = (negs, pooled, grid, X, posl)
        perm = np.random.default_rng([SP.SPLIT_SEED, si, sd, d_]).permutation(m)
        cache = {}

        def ubf(level, k, j, idx=None):
            key = (level, k, j, None if idx is None else tuple(idx))
            if key not in cache: cache[key] = L.b1u(D[k][3][:, j] if idx is None else D[k][3][idx, j], ni if idx is None else ni[idx], level)
            return cache[key]

        def rec(name, eps, k, t, cert, extra=None):
            sn = sens_at(k, t)
            rows.append(dict(site=site, seed=sd, draw=d_, eps=eps, rule=name, cand=k, cert=bool(cert), fpr=float(fpr_at(k, t)), sens=float(sn),
                             dse=float(sn - pop["orac_unc"][eps]), dso=float(sn - pop["orac_own"][(k, eps)]), ni=ni.tolist(), fields=[str(f) for f in F], **(extra or {})))

        def choose(T, eps, name, extra=None):
            cert_c = [k for k in L.CEN if T[k][1]]
            if not cert_c: rec(name, eps, "c0", T["c0"][0], False, (extra or {}).get("c0")); return
            if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
            else:
                sl = [(D[k][4] >= T[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
            rec(name, eps, ch, T[ch][0], True, (extra or {}).get(ch))
        for eps, dl in L.SETTINGS:
            choose({k: L.fixed_sequence(D[k][2], lambda j, k=k: ubf(0.10 / 3, k, j), eps) for k in L.CEN}, eps, "S3")
            choose({k: L.fixed_sequence(D[k][2], lambda j, k=k: ubf(0.10 / 6, k, j), eps) for k in L.CEN}, eps, "S3c")
            for k in ("c0", "c1s"):
                t_, c_, _ = L.fixed_sequence(D[k][2], lambda j, k=k: ubf(0.10, k, j), eps); rec("B1u_" + k, eps, k, t_, c_)
            cl = {k: AC[k].c_lab(F) for k in L.CEN}
            for ms in (6,):
                arm = f"SPLIT{ms}"; Ls = sorted(perm[:ms].tolist()); Lc = sorted(perm[ms:].tolist()); mc = m - ms; Fs = [F[i] for i in Ls]; Fc = [F[i] for i in Lc]
                Ts, extra = {}, {}
                for k in L.CEN:
                    RS = AC[k].F - ms
                    pooled_s = np.concatenate([D[k][0][i] for i in Ls]); ts = float(np.quantile(pooled_s, 1 - eps))
                    xs_t = np.array([ni[i] - np.searchsorted(D[k][0][i], ts, side="left") for i in Ls], float); ns_ = ni[Ls].astype(float)
                    ws_t = AC[k].g_lab(Fs, ts) - eps * cl[k][Ls]; es_t = (xs_t - eps * ns_) - ws_t; vs_t = QC[k].v_fields(Fs, ts, eps)
                    b, Wd = SP.select_bound(es_t, vs_t, xs_t, ns_, mc, RS); SP.CHOSEN[(site, eps, arm, k)][b] += 1
                    Wc = np.random.default_rng([E.BOOT_SEED, si, sd, d_, int(round(eps * 100)), E.CAND_CODE[k], 1000 + ms]).multinomial(mc, np.full(mc, 1 / mc), size=E.B).astype(float)

                    def split_conditions(j, k=k, b=b, Ls=Ls, Lc=Lc, Fs=Fs, Fc=Fc, RS=RS, mc=mc, Wc=Wc, arm=arm):
                        t = D[k][2][j]; xs = D[k][3][Ls, j]; c2 = float(xs.sum() / ni[Ls].sum()) <= eps
                        if b == "b1u":
                            c1 = ubf(A3, k, j, Lc) <= eps
                        else:
                            Sw, _ = AC[k].w_stats(t, eps); w_all = AC[k].g_lab(F, t) - eps * cl[k]; z_all = D[k][3][:, j] - eps * ni
                            SwR = Sw - float(w_all[Ls].sum()); e_c = z_all[Lc] - w_all[Lc]
                            if b == "K1": U = VB.k1_ub_alpha(e_c, SwR, RS, mc, Wc, A3, (site, eps, arm))
                            else:
                                v_all = QC[k].v_fields(F, t, eps); SvR = QC[k].sum_v(t, eps) - float(v_all[Ls].sum())
                                U = VB.k1_ub_alpha(e_c + v_all[Lc], SwR - SvR, RS, mc, Wc, A3, (site, eps, arm))
                            c1 = (U is not None) and np.isfinite(U) and U <= 0
                        return bool(c1), bool(c2)

                    memo = {}

                    def seq_fn(j, sc=split_conditions, key=(site, eps, arm, k), memo=memo):
                        c1, c2 = sc(j); memo[j] = (c1, c2)
                        if c1 and not c2: SP.C2STOP[key] += 1
                        return SP.split_to_seq(c1, c2)
                    t_, c_, jj = L.fixed_sequence(D[k][2], E.counted(arm, site, eps, seq_fn), eps); Ts[k] = (t_, c_, jj)
                    ps_at = float(D[k][3][Ls, jj].sum() / ni[Ls].sum()) if jj >= 0 else None
                    extra[k] = dict(split_bound=b, split_ps_at=ps_at, split_c1_at=(memo[jj][0] if jj >= 0 else None), split_widths={kk: float(v) for kk, v in Wd.items()})
                choose({k: Ts[k][:2] for k in L.CEN}, eps, arm, extra)
    return rows


def band(c): return "hold_1-delta" if c >= 0.90 else ("below_1-delta_within_1-2delta" if c >= 0.80 else "fail")
def near(c): return bool(abs(c - 0.90) <= 0.01 or abs(c - 0.80) <= 0.01)


def summarise(EB, rows_by, rare_by, eps_list=EPSS):
    out = {}; sds = list(rows_by)
    for eps in eps_list:
        for rule in ("S3", "S3c", "SPLIT6", "B1u_c1s"):
            sel = {sd: [r for r in rows_by[sd] if r["rule"] == rule and abs(r["eps"] - eps) < 1e-12] for sd in sds}
            cov = [float(np.mean([r["fpr"] <= eps for r in sel[sd]])) for sd in sds]; cc = [float(1 - np.mean([r["cert"] for r in sel[sd]])) for sd in sds]; dso = [float(np.mean([r["dso"] for r in sel[sd]])) for sd in sds]
            m = float(np.mean(cov)); out[f"{rule}|{eps}"] = dict(coverage=m, coverage_per_realisation=cov, coverage_ci90_seedcluster=EB.seed_cluster_ci(cov), band=band(m), near_boundary=near(m), cannot_certify=float(np.mean(cc)), ds_own=float(np.mean(dso)))
        a = [np.mean([r["fpr"] <= eps for r in rows_by[sd] if r["rule"] == "S3c" and abs(r["eps"] - eps) < 1e-12]) - np.mean([r["fpr"] <= eps for r in rows_by[sd] if r["rule"] == "S3" and abs(r["eps"] - eps) < 1e-12]) for sd in sds]
        b = [np.mean([r["dso"] for r in rows_by[sd] if r["rule"] == "S3c" and abs(r["eps"] - eps) < 1e-12]) - np.mean([r["dso"] for r in rows_by[sd] if r["rule"] == "S3" and abs(r["eps"] - eps) < 1e-12]) for sd in sds]
        out[f"S3c-S3|{eps}"] = dict(delta_coverage=float(np.mean(a)), delta_ds_own=float(np.mean(b)), delta_coverage_per_realisation=[float(x) for x in a], delta_ds_own_per_realisation=[float(x) for x in b])
        if any(len(v) for v in rare_by.values()):
            for rule in ("S3", "S3c"):
                has, hasnot = [], []
                for sd in sds:
                    rf = set(int(x) for x in rare_by[sd])
                    for r in rows_by[sd]:
                        if r["rule"] == rule and abs(r["eps"] - eps) < 1e-12: (has if any(int(f) in rf for f in r["fields"]) else hasnot).append(r["fpr"] <= eps)
                out[f"rare_split|{rule}|{eps}"] = dict(share_draws_with_rare_field=float(len(has) / max(1, len(has) + len(hasnot))), coverage_with_rare=(float(np.mean(has)) if has else None), coverage_without_rare=(float(np.mean(hasnot)) if hasnot else None), n_with=len(has), n_without=len(hasnot))
    return out


def check_s3c_s3(rows):
    cert = {(r["seed"], r["draw"], r["eps"], r["rule"]): r["cert"] for r in rows}
    return [(sd, d, e) for (sd, d, e, rule), c in cert.items() if rule == "S3c" and c and not cert[(sd, d, e, "S3")]]


def main(argv=None, _hooks=None):
    """_hooks: test seam used only by e12_smoke_local.py ({'skip_pins': True}); None in every real run."""
    argv = sys.argv[1:] if argv is None else argv; gate(argv)
    ap = argparse.ArgumentParser(); ap.add_argument("--e10_json", required=True); ap.add_argument("--out", required=True); ap.add_argument("--cells", default=None); ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--ndraw", type=int, default=200); ap.add_argument("--check_draws", type=int, default=2000); a = ap.parse_args(argv); out = Path(a.out)
    if out.exists(): raise SystemExit("output directory exists: use a new directory")
    sys.path[:0] = [str(R / "e10"), str(R / "e_batch"), str(R / "m3" / "split"), str(R / "m3" / "k2prep"), str(R / "c5b")]
    import e10_synth as E10, ebatch_lib as EB, m3_replay_lib as L, m3_split_dev_eval as SP, m3_k1v_precheck as KV, s8_bounds as B8
    gate(argv)
    if not (_hooks and _hooks.get("skip_pins")):
        for name, mod in (("m3_replay_lib", L), ("m3_split_dev_eval", SP), ("s8_bounds", B8), ("e10_synth", E10), ("ebatch_lib", EB)):
            if sha(mod.__file__) != PINS[name]: raise AnchorError((name, "differs from its pin"))
    assert SP.L is L; L.M = M; L.NDRAW = a.ndraw; site = "bbbc041_all"; assert site in L.SITECFG and site in L.DEV
    t0 = time.time(); grid = cells(); sel = [int(x) for x in a.cells.split(",")] if a.cells else list(range(len(grid)))
    class Synth:
        used = {}
        def __init__(self, sites): self.sites = sites
        def load_site(self, sd, st): return self.sites[sd]
    # (a) generator checks of E10
    chk = E10.generator_checks(np.random.default_rng([20261108, 12, 9999, 0, 0]), L, B8, ndraw=a.check_draws)
    # (b) neutral make_site2 == make_site, bit for bit
    for ci_, kw in ((0, dict(tau=0.5, r=100, pi=0.05, tail="normal", F=300)), (0, dict(tau=0.0, r=20, pi=0.2, tail="t5", F=60))):
        x1 = E10.make_site(np.random.default_rng([1, 2, 3]), kw["tau"], kw["r"], kw["pi"], kw["tail"], F=kw["F"]); x2 = make_site2(E10, np.random.default_rng([1, 2, 3]), **kw)[:3]
        if not all(np.array_equal(u, v) for u, v in zip(x1, x2)): raise AnchorError(("make_site2 with neutral parameters differs from make_site", kw))
    # (c) uniform sampler: copy of run_site_seed equals the library
    s0, y0, g0 = E10.make_site(np.random.default_rng(stream(0, 41, 0)), tau=0.5, r=100, pi=0.05, tail="normal"); ss0, ys0, _ = E10.source_set(np.random.default_rng(stream(0, 41, 1)), 100, 0.05, "normal")
    pl0, _ = E10.platt_synth(KV, ss0, ys0); keep = L.NDRAW; L.NDRAW = 5
    try:
        lib = SP.run_site_seed(Synth({41: (s0, y0, g0)}), 41, site, pl0); mine = run_site_seed_e12(SP, Synth({41: (s0, y0, g0)}), 41, site, pl0, uniform_sampler)
    finally: L.NDRAW = keep
    core = ("cert", "cand", "fpr", "sens", "dse", "dso"); kk = lambda r: (r["draw"], r["eps"], r["rule"])
    A = {kk(r): r for r in lib if r["rule"] in ("S3", "B1u_c0", "B1u_c1s", "SPLIT6")}; B = {kk(r): r for r in mine if r["rule"] in ("S3", "B1u_c0", "B1u_c1s", "SPLIT6")}
    if set(A) != set(B) or not all(all(A[k][f] == B[k][f] for f in core) for k in A): raise AnchorError(("copy of run_site_seed differs from the library", len(A), len(B)))
    res_cells = {}; rows_out = []; platt_info = {}
    for ci in sel:
        cell = grid[ci]; sites, platt, icc, rare_by, rows_by = {}, {}, {}, {}, {}
        for rid in range(a.reps):
            sd = 41 + rid
            s, y, g, slide_of, rare_f = make_site2(E10, np.random.default_rng(stream(ci, sd, 0)), cell["tau"], cell["r"], cell["pi"], cell["tail"], cell["F"], cell["size_sigma"], cell["rare"], cell["slides"])
            sites[sd] = (s, y, g, slide_of); rare_by[sd] = rare_f; icc[sd] = E10.realised_icc(s, y, g)
            ss, ys, _ = E10.source_set(np.random.default_rng(stream(ci, sd, 1)), cell["r"], cell["pi"], cell["tail"]); platt[sd], pinfo = E10.platt_synth(KV, ss, ys); platt_info[(ci, sd)] = pinfo
        SP.E.EXC.clear(); SP.E.DEG.clear(); SP.E.DROP.clear(); SP.CHOSEN.clear(); SP.C2STOP.clear()
        for sd in sites:
            smp = uniform_sampler if cell["two_stage"] is None else two_stage_sampler(sites[sd][3], cell["slides"][0], cell["two_stage"])
            rows_by[sd] = run_site_seed_e12(SP, Synth({sd: sites[sd][:3]}), sd, site, platt[sd], smp)
            bad = check_s3c_s3(rows_by[sd])
            if bad: raise AnchorError(("e: S3c certified a draw that S3 did not", ci, sd, bad[:3]))
        summ = summarise(EB, rows_by, rare_by)
        res_cells[str(ci)] = dict(cell={k: v for k, v in cell.items()}, realised_icc_per_realisation=[icc[sd] for sd in sites], realised_icc_mean=float(np.mean(list(icc.values()))), results=summ)
        for sd in sites:
            for r in rows_by[sd]: rows_out.append({k: v for k, v in r.items() if k not in ("ni", "fields") and not k.startswith("split_")} | dict(cell=ci))
    # (d) base cell equals E10 cell 0
    anchor_d = "not_applicable (cell 0 not run or non-default size)"
    if 0 in sel and a.reps == 5 and a.ndraw == 200:
        e10 = json.load(open(a.e10_json)); assert e10["reps"] == 5 and e10["ndraw"] == 200, "E10 size"; lr = e10["cells"]["0"]["library_rules"]
        for eps in EPSS:
            for rule in ("S3", "SPLIT6"):
                mine_c = res_cells["0"]["results"][f"{rule}|{eps}"]["coverage_per_realisation"]; ref_c = lr[f"{rule}|{eps}"]["coverage_per_realisation"]
                if len(mine_c) != len(ref_c) or any(abs(u - v) > 1e-12 for u, v in zip(mine_c, ref_c)): raise AnchorError(("d: base cell differs from E10 cell 0", rule, eps, mine_c, ref_c))
        anchor_d = "passed (S3 and SPLIT6 per-realisation coverage at eps 0.10 and 0.05 equal E10 cell 0 within 1e-12)"
    flags = [dict(cell=int(ci), name=c["cell"]["name"], key=k, coverage=v["coverage"], band=v["band"], near_boundary=v["near_boundary"]) for ci, c in res_cells.items() for k, v in c["results"].items()
             if k.split("|")[0] in ("S3", "S3c", "SPLIT6", "B1u_c1s") and isinstance(v, dict) and "band" in v and v["band"] != "hold_1-delta"]
    rep = dict(version="e12_r1", spec="P5_规格_E11E12_水平修正与压力测试_开发站与合成_20261108_draft2", script_sha256=sha(__file__), pins=PINS, reps=a.reps, ndraw=a.ndraw, m=M, levels=LEV,
               anchors=dict(a_generator_checks=chk, b="passed (make_site2 neutral == make_site bit for bit)", c="passed (copy of run_site_seed == library for S3, B1u_c0, B1u_c1s, SPLIT6, 5 draws)", d=anchor_d, e="passed (S3c certified => S3 certified, every draw)"),
               cells=res_cells, cells_below_hold_1_delta=flags, platt_fallback=[[ci, sd] for (ci, sd), v in platt_info.items() if v["fallback"]], seconds=time.time() - t0, ci_definition=EB.CI_DEFINITION,
               note="synthetic; no project data; supplementary stress test; no rule or constant changed; every cell is reported")
    out.mkdir(parents=False, exist_ok=False)
    files = {"e12_summary.json": json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False), "e12_rows.json": json.dumps(rows_out, allow_nan=False)}
    for n, txt in files.items():
        with open(out / n, "x") as fh: fh.write(txt)
    with open(out / "SHA256SUMS", "x") as fh: fh.write("".join(f"{sha(out / n)}  {n}\n" for n in files))
    print(out)


if __name__ == "__main__":
    try: main()
    except Exception as e:   # AnchorError of this script and EB.AnchorError of ebatch_lib (raised by E10.generator_checks) are both anchor mismatches
        if type(e).__name__ == "AnchorError": print("ANCHOR MISMATCH (exit 3, no output):", e, file=sys.stderr); sys.exit(3)
        raise

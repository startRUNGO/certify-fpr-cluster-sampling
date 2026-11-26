# P5 SPLIT development-site evaluation (problem card draft7.1 d2c5a182; theory draft2 86a4c9e4 Proposition 2' + addendum A a876bf8f).  LOCAL, read-only.
# WRITTEN, NOT RUN, until 实验管理 approves.
# NATURE: development-site mechanism selection on 4 sites (bbbc041_all, mpidb_full, lacuna_B, old MP-IDB); not confirmatory evidence; no M5, IML, PlasmoCount, NIH test.
# ARMS: S3, B1u_c1s (= c1s), K1 (draft5), K1vB (draft6) as controls; SPLIT6 (main); SPLIT4, SPLIT8 (sensitivity, report only).
# SPLIT_m_s, per draw:
#   split: perm = default_rng([20261102, site index, seed, draw]).permutation(22) of the 22 drawn fields; L_s = first m_s, L_c = rest (nested across m_s 4/6/8);
#          this stream is separate from, and does not consume, the REV-1 stream.
#   choose (per candidate k, per eps), on L_s labels only: t_s = quantile(pooled L_s negatives of k, 1 - eps); at t_s
#          W_K1 = z* sd(e_s)/(sqrt(m_c) nbar_s) fpc_R;  W_Kv = z* sd(e_s + v_s)/(sqrt(m_c) nbar_s) fpc_R;
#          W_b1u = z* sqrt(p_s (1 - p_s)(1 + (nbar_s - 1) ICC_s)/(m_c nbar_s)), ICC_s = one-way ANOVA point estimate on L_s clipped to [0, 0.999];
#          z* = Phi^-1(1 - delta/3), fpc_R = sqrt(1 - m_c/|R|), R = F minus L_s, sd with ddof 1, nbar_s = mean negatives per L_s field;
#          argmin, ties b1u > K1 > K1-v.
#   certify (Proposition 2'): at each grid point t of the S3 grid (quantiles of the pooled negatives of all 22 fields; addendum A2), split_conditions(j) returns
#          (C1, C2): C1 = the chosen bound on L_c with population R certifies p_R(t) <= eps  [K1: k1_ub_alpha(e_c, sum_R w, |R|, m_c, W_c, delta/3) <= 0;
#          K1-v: same with residual e + v and sum_R (w - v); b1u: b1u(x_c, n_c, delta/3) <= eps];  C2 = p_s(t) = sum_{L_s} x / sum_{L_s} n <= eps.
#          split_to_seq(C1, C2) = 0.0 iff both hold, else 1.0 (written fixed, as to_seq); fixed-sequence test stops at the first point where either fails.
#          Union over c0, c1s, c2s with delta/3 each; candidate choice as S3 (labelled positives of all 22 fields).
#   sums over R are exact: sum_R w = sum_F w (AllCells) - sum_{L_s} w_i; sum_R v = sum_F v (prefix sums of q over all cells) - sum_{L_s} v_i.
#   q: frozen Platt per seed on the NIH validation fold (m3_k1v_precheck.load_nihval / platt; keys logits, y only), applied to the raw score s.
#   bootstrap weights for L_c: default_rng([20261101, site index, seed, draw, eps code, candidate code, 1000 + m_s]).multinomial(m_c, ..., B = 999).
# CHECKS before writing: replay self-check (S3, B1u_c0, B1u_c1s of this loop vs stored REV-1, 4 sites, 1e-12).  Rows record, at the chosen threshold, p_s and the
#   chosen bound's statistic so that the smoke can verify C1 and C2 for every certified draw.
# CRITERIA (draft6 section 4, as fixed in draft7.1 section 4) for SPLIT6 against S3: 8-cell gate (coverage >= 0.90, 90% seed-cluster lower >= 0.88, cannot_certify
#   <= S3 + 0.02, no flag); mean dSens_own improvement over the 3 rebuilt sites at eps 0.10 >= +0.02 with none of those 6 cells worse than S3 by > 0.01; old MP-IDB both
#   eps not worse than S3 by > 0.01.  Fail -> root-cause analysis; criteria and m_s unchanged.
# Output: one new JSON in a new directory ("x", allow_nan=False).  Prints only the output path.
# Usage: python -I -B m3_split_dev_eval.py <m3_stress inputs> --nihval_dir <runs/c1a6_stage2_20261027T133543Z/a2ref> --out <dir/json>
import argparse, hashlib, json, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from scipy.special import expit
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
for sub in ("k2prep", "stress", "k1vb"): sys.path.insert(0, str(ROOT / sub))
import m3_replay_lib as L   # noqa: E402
import m3_k12_dev_eval as E   # noqa: E402
import m3_stress as ST   # noqa: E402
import m3_k1vb_dev_eval as VB   # noqa: E402
import m3_k1v_precheck as KV   # noqa: E402

DELTA = 0.10; A3 = DELTA/3; A6 = DELTA/6; OLD = ST.OLD; REBUILT = ("bbbc041_all", "mpidb_full", "lacuna_B"); SITES = REBUILT + (OLD,)
SPLITS = (4, 6, 8); MAIN = "SPLIT6"; ARMS = ("S3", "B1u_c1s", "K1", "K1vB") + tuple(f"SPLIT{ms}" for ms in SPLITS); ORDER = ("b1u", "K1", "K1-v")
SPLIT_SEED = 20261102; GATE_COV, GATE_LO, GATE_CC, MAIN_GAIN, MAIN_WORSE = 0.90, 0.88, 0.02, 0.02, 0.01
CHOSEN = defaultdict(Counter); C2STOP = defaultdict(int)


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def split_to_seq(c1, c2):
    """Fixed interface (实验管理 2026-10-30): 0.0 iff C1 and C2 both hold at this grid point, else 1.0."""
    return 0.0 if (bool(c1) and bool(c2)) else 1.0


def anova_icc(x, n):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum()
    if m < 2 or N <= m: return 0.0
    ph = x.sum()/N; pi = x/n; n0 = (N - (n**2).sum()/N)/(m - 1); msb = (n*(pi - ph)**2).sum()/(m - 1); msw = (n*pi*(1 - pi)).sum()/(N - m)
    if msw <= 0: return 0.0 if msb <= 0 else 0.999
    th = (msb/msw - 1)/n0; return float(min(max(th/(1 + th), 0.0), 0.999)) if th > 0 else 0.0


def select_bound(e_s, v_s, x_s, n_s, m_c, R_size):
    zs = norm.ppf(1 - A3); nb = float(np.mean(n_s)); fpc = np.sqrt(max(1 - m_c/R_size, 0.0)); den = np.sqrt(m_c)*nb
    W = {"K1": zs*np.std(e_s, ddof=1)/den*fpc, "K1-v": zs*np.std(e_s + v_s, ddof=1)/den*fpc}
    p = float(np.sum(x_s)/np.sum(n_s)); icc = anova_icc(x_s, n_s); W["b1u"] = zs*np.sqrt(max(p*(1 - p), 0.0)*(1 + (nb - 1)*icc)/(m_c*nb))
    return min(ORDER, key=lambda b: (W[b], ORDER.index(b))), W


class QCells:
    """Exact sum over all cells of the site of q (1{s_k >= t} - eps), and per-field v_i, for one candidate."""
    def __init__(self, sk, q, g, fl):
        keep = np.isin(g, np.asarray(fl)); s = sk[keep]; qq = q[keep]; o = np.argsort(-s, kind="stable")
        self.desc = s[o]; self.cq = np.concatenate([[0.0], np.cumsum(qq[o])]); self.Q = float(qq.sum()); self.asc = np.sort(s); self.n = len(s)
        self.field = {f: (np.sort(sk[g == f])[::-1], np.concatenate([[0.0], np.cumsum(q[g == f][np.argsort(-sk[g == f], kind="stable")])]), float(q[g == f].sum())) for f in fl}

    def sum_v(self, t, eps):
        K = self.n - np.searchsorted(self.asc, t, side="left"); return float(self.cq[K] - eps*self.Q)

    def v_fields(self, fields, t, eps):
        out = []
        for f in fields:
            sd, cq, Qf = self.field[f]; K = int(np.searchsorted(-sd, -t, side="right")); out.append(cq[K] - eps*Qf)
        return np.array(out)


def run_site_seed(inp, sd, site, platt):
    pop = L.population(inp, sd, site); fl, fl_neg, fl_posv = pop["fl"], pop["fl_neg"], pop["fl_posv"]; sens_at, fpr_at = pop["sens_at"], pop["fpr_at"]
    AC = {k: E.AllCells(pop, k) for k in L.CEN}; q = expit(platt[0]*pop["s"] + platt[1]); QC = {k: QCells(pop["sc"][k], q, pop["g"], fl) for k in L.CEN}
    key8, sidx = L.SITECFG[site]; ai = L.AI; rows = []; m = L.M; si = L.DEV.index(site)
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
        perm = np.random.default_rng([SPLIT_SEED, si, sd, d_]).permutation(m)
        cache = {}

        def ubf(level, k, j, idx=None):
            key = (level, k, j, None if idx is None else tuple(idx))
            if key not in cache:
                cache[key] = L.b1u(D[k][3][:, j] if idx is None else D[k][3][idx, j], ni if idx is None else ni[idx], level)
            return cache[key]

        def rec(name, eps, k, t, cert, extra=None):
            sn = sens_at(k, t)
            rows.append(dict(site=site, seed=sd, draw=d_, eps=eps, rule=name, cand=k, cert=bool(cert), fpr=float(fpr_at(k, t)), sens=float(sn),
                             dse=float(sn-pop["orac_unc"][eps]), dso=float(sn-pop["orac_own"][(k, eps)]), ni=ni.tolist(), fields=[str(f) for f in F], **(extra or {})))

        def choose(T, eps, name, extra=None):
            cert_c = [k for k in L.CEN if T[k][1]]
            if not cert_c: rec(name, eps, "c0", T["c0"][0], False, (extra or {}).get("c0")); return
            if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
            else:
                sl = [(D[k][4] >= T[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
            rec(name, eps, ch, T[ch][0], True, (extra or {}).get(ch))
        for eps, dl in L.SETTINGS:
            choose({k: L.fixed_sequence(D[k][2], lambda j, k=k: ubf(0.10/3, k, j), eps) for k in L.CEN}, eps, "S3")
            for k in ("c0", "c1s"):
                t_, c_, _ = L.fixed_sequence(D[k][2], lambda j, k=k: ubf(0.10, k, j), eps); rec("B1u_" + k, eps, k, t_, c_)
            Wk = {k: np.random.default_rng([E.BOOT_SEED, si, sd, d_, int(round(eps*100)), E.CAND_CODE[k]]).multinomial(m, np.full(m, 1/m), size=E.B).astype(float) for k in L.CEN}
            cl = {k: AC[k].c_lab(F) for k in L.CEN}

            def k1_U(j, k, alpha, arm, eps=eps):
                t = D[k][2][j]; z = D[k][3][:, j] - eps*ni; w = AC[k].g_lab(F, t) - eps*cl[k]; Sw, _ = AC[k].w_stats(t, eps)
                return VB.k1_ub_alpha(z - w, Sw, AC[k].F, m, Wk[k], alpha, (site, eps, arm))

            def ub_vb(j, k, eps=eps):
                if ubf(A6, k, j) <= eps: return 0.0
                return E.to_seq(k1_U(j, k, A6, "K1vB"))
            choose({k: L.fixed_sequence(D[k][2], E.counted("K1", site, eps, lambda j, k=k: E.to_seq(k1_U(j, k, A3, "K1"))), eps) for k in L.CEN}, eps, "K1")
            choose({k: L.fixed_sequence(D[k][2], E.counted("K1vB", site, eps, lambda j, k=k: ub_vb(j, k)), eps) for k in L.CEN}, eps, "K1vB")
            for ms in SPLITS:
                arm = f"SPLIT{ms}"; Ls = sorted(perm[:ms].tolist()); Lc = sorted(perm[ms:].tolist()); mc = m - ms; Fs = [F[i] for i in Ls]; Fc = [F[i] for i in Lc]
                Ts, extra = {}, {}
                for k in L.CEN:
                    RS = AC[k].F - ms
                    pooled_s = np.concatenate([D[k][0][i] for i in Ls]); ts = float(np.quantile(pooled_s, 1 - eps))
                    xs_t = np.array([ni[i] - np.searchsorted(D[k][0][i], ts, side="left") for i in Ls], float); ns_ = ni[Ls].astype(float)
                    ws_t = AC[k].g_lab(Fs, ts) - eps*cl[k][Ls]; es_t = (xs_t - eps*ns_) - ws_t; vs_t = QC[k].v_fields(Fs, ts, eps)
                    b, Wd = select_bound(es_t, vs_t, xs_t, ns_, mc, RS); CHOSEN[(site, eps, arm, k)][b] += 1
                    Wc = np.random.default_rng([E.BOOT_SEED, si, sd, d_, int(round(eps*100)), E.CAND_CODE[k], 1000 + ms]).multinomial(mc, np.full(mc, 1/mc), size=E.B).astype(float)

                    def split_conditions(j, k=k, b=b, Ls=Ls, Lc=Lc, Fs=Fs, Fc=Fc, RS=RS, mc=mc, Wc=Wc, arm=arm):
                        t = D[k][2][j]; xs = D[k][3][Ls, j]; c2 = float(xs.sum()/ni[Ls].sum()) <= eps
                        if b == "b1u":
                            c1 = ubf(A3, k, j, Lc) <= eps
                        else:
                            Sw, _ = AC[k].w_stats(t, eps); w_all = AC[k].g_lab(F, t) - eps*cl[k]; z_all = D[k][3][:, j] - eps*ni
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
                        if c1 and not c2: C2STOP[key] += 1
                        return split_to_seq(c1, c2)
                    t_, c_, jj = L.fixed_sequence(D[k][2], E.counted(arm, site, eps, seq_fn), eps); Ts[k] = (t_, c_, jj)
                    ps_at = float(D[k][3][Ls, jj].sum()/ni[Ls].sum()) if jj >= 0 else None
                    extra[k] = dict(split_bound=b, split_ps_at=ps_at, split_c1_at=(memo[jj][0] if jj >= 0 else None), split_widths={kk: float(v) for kk, v in Wd.items()})
                choose({k: Ts[k][:2] for k in L.CEN}, eps, arm, extra)
    return rows


def judge(cells, arm):
    get = lambda a, s, e: [c for c in cells if c["arm"] == a and c["site"] == s and abs(c["eps"] - e) < 1e-12][0]
    gate, diffs = [], []
    for s in SITES:
        for e, _ in L.SETTINGS:
            a, b = get(arm, s, e), get("S3", s, e)
            gate.append(dict(site=s, eps=e, cov_ok=a["coverage"] >= GATE_COV, lo_ok=a["coverage_ci90_seedcluster"][0] >= GATE_LO,
                             cc_ok=a["cannot_certify"] <= b["cannot_certify"] + GATE_CC + 1e-12, flagged=a["flag"] is not None))
            diffs.append(dict(site=s, eps=e, diff=a["ds_own"] - b["ds_own"]))
    g_ok = all(x["cov_ok"] and x["lo_ok"] and x["cc_ok"] and not x["flagged"] for x in gate)
    mg = float(np.mean([d["diff"] for d in diffs if d["site"] in REBUILT and abs(d["eps"] - 0.10) < 1e-12]))
    r_ok = all(d["diff"] >= -MAIN_WORSE for d in diffs if d["site"] in REBUILT); o_ok = all(d["diff"] >= -MAIN_WORSE for d in diffs if d["site"] == OLD)
    m_ok = mg >= MAIN_GAIN and r_ok and o_ok
    return dict(arm=arm, gate=gate, gate_pass=g_ok, mean_gain_rebuilt_eps010=mg, per_cell_diff=diffs, rebuilt_no_cell_worse=r_ok, old_not_worse=o_ok,
                main_pass=m_ok, passes=bool(g_ok and m_ok))


def main():
    ap = argparse.ArgumentParser()
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "old_s41", "old_dir", "nihval_dir", "out"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args(); out = Path(a.out); assert not out.exists(), "output exists"
    inp = ST.StressInputs(a); assert sorted(inp.old) == [41, 42, 43, 44, 45], "old MP-IDB needs all 5 seeds"
    nv = KV.load_nihval(a.nihval_dir); platt = {sd: KV.platt(*nv[sd]) for sd in L.SEEDS}
    L.SITECFG[OLD] = (8, 1); L.DEV = tuple(L.DEV) + (OLD,); L.M = 22
    E.EXC.clear(); E.DEG.clear(); E.DROP.clear(); CHOSEN.clear(); C2STOP.clear(); rows = []
    for site in SITES:
        for sd in L.SEEDS: rows += run_site_seed(inp, sd, site, platt[sd])
    status, bad = ST.check(inp, rows, SITES, 22)
    if bad or not all(v.startswith("self-check passed") for v in status.values()):
        raise SystemExit("self-check failed or not possible, nothing written: " + repr((bad[:5], status)))
    rng = np.random.default_rng(E.CI_SEED)
    cells = [dict(arm=arm, site=s, eps=e, **E.summarise([r for r in rows if r["site"] == s and abs(r["eps"] - e) < 1e-12 and r["rule"] == arm], e, s, arm, rng))
             for arm in ARMS for s in SITES for e, _ in L.SETTINGS]
    for c in cells:
        if c["arm"].startswith("SPLIT"):
            c["bound_chosen"] = {k: dict(CHOSEN[(c["site"], c["eps"], c["arm"], k)]) for k in L.CEN}
            c["c2_bottleneck_points"] = {k: int(C2STOP[(c["site"], c["eps"], c["arm"], k)]) for k in L.CEN}
    rep = dict(script_sha256=sha(__file__), lib_sha256=sha(L.__file__), eval_sha256=sha(E.__file__), stress_sha256=sha(ST.__file__), vb_sha256=sha(VB.__file__),
               precheck_sha256=sha(KV.__file__), inputs_sha256=dict(L.PIN), score_files_sha256=inp.used, platt={str(k): list(v) for k, v in platt.items()},
               problem_card="draft7.1 d2c5a182; theory draft2 86a4c9e4 + addendum A a876bf8f", nature="development-site mechanism selection; not confirmatory evidence",
               settings=dict(delta=DELTA, alpha=A3, split_seed=SPLIT_SEED, m_s_main=6, m_s_sensitivity=[4, 8], B=E.B, boot_seed=E.BOOT_SEED, ci_B=E.CI_B, ci_seed=E.CI_SEED),
               self_check=status, cells=cells, judgement=judge(cells, MAIN), sensitivity_judgement_report_only=[judge(cells, f"SPLIT{ms}") for ms in (4, 8)])
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False)
    out.parent.mkdir(parents=False, exist_ok=False)
    with open(out, "x") as fh: fh.write(txt)
    print("written", out)


if __name__ == "__main__":
    main()

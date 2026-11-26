# P5 M3/K2 preparation: shared replay library for R2 (m3_r2_aux.py), R3 (m3_r3_bound.py) and W3 (m3_w3_baselines.py).  LOCAL, read-only.  WRITTEN, NOT RUN,
# until 实验管理 approves (problem card draft4 fcf90632).
# What it replays: the REV-1 m = 22 draws of arm a2 on the three DEVELOPMENT sites, with the same score files, random streams, centring, grid and b1u as the runs
#   behind the manuscript.  The setup and draw loop are copied from results_src/diag_icc0/d_icc0.py (8e32f299, itself a verbatim copy of the REV-1 task body, whose
#   S3 aggregates matched the stored REV-1 results to 1e-12); b1u is asserted line for line against results_src/mh/mh_rev1.py (47f79f1d).
# Data boundary: development sites only (bbbc041_all, mpidb_full, lacuna_B).  Score files are checked against their manifests before use (as m3_rc.py 57993952):
#   bbbc041_all <a2link>/a2_s{sd}_bbbc041_all_logits.npz vs MANIFEST.tsv; labels and fields from ids_iter1.npz, keys bbbc041_all_y and bbbc041_all_group ONLY;
#   mpidb_full / lacuna_B <mhlh>/a2_s{sd}_{mh,lh}_logits.npz and ids_{mh,lh}.npz vs SHA256SUMS.  No IML, no PlasmoCount, no NIH key of ids_iter1.npz.
# Self-check (callers must run it before writing anything): for S3 (rule U3, delta/3, union, labelled choice), B1u_c0 and B1u_c1s (rule U1, delta, one candidate):
#   cov, fpr_mean, cannot_certify, ds_eps, ds_own and n recomputed here equal the stored REV-1 aggregates within 1e-12, 3 sites x 2 eps.
import hashlib, json
from pathlib import Path
import numpy as np
from scipy.stats import beta as betad, f as fdist

PIN = {"manifest": "94df12d55a1bacc048ff8f7ae5b279bee36f7ffc8941958c1527794d2c36bdaf", "sums": "7fb60c41f83d8bec2edb2c4710f1b21ec8e3ef2a71678a304870d9e9e12e0f49", "ids_iter1": "75b46e6ff2d58c1d4ae3b4f4e1436691398acedf4f8dc13e64c40f530fa60b26", "rev1_bbbc": "bd4f45fbf224a3975ed04cdd7f6e7dc851ea2a057c1a051a8c85e627ac8ced41", "rev1_mh": "e87ca813427c3c0d178c3908bab5737fb0eec64f724c23a5f1d9e83d79090243",
       "rev1_lh": "396b91ab957dcbbbfe1b4c9f717310213c9b81d41fd1d79163f4b4ad0ab6ccf8", "mh_rev1_py": "47f79f1df632e5efd1be5ea9485826d082eeaac38f89ee6f60a6aed2292116c5"}
SEEDS = [41, 42, 43, 44, 45]; NDRAW = 200; SETTINGS = [(0.10, 0.10), (0.05, 0.10)]; G = 200; CEN = ["c0", "c1s", "c2s"]; M = 22; AI = 0; ARM = "a2"
DEV = ("bbbc041_all", "mpidb_full", "lacuna_B")
SITECFG = {"bbbc041_all": (8, 0), "mpidb_full": (31, 0), "lacuna_B": (8, 2)}
STORED_KEY = {"bbbc041_all": "a2|bbbc041_all|22|{rule}|{eps}", "mpidb_full": "a2|mh|22|{rule}|{eps}", "lacuna_B": "a2|lh|22|{rule}|{eps}"}
SHORT = {"mpidb_full": "mh", "lacuna_B": "lh"}


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


# b1u below is verbatim from mh_rev1.py; check_b1u_copy() asserts the copy line for line
def b1u(x, n, delta):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X/N
    n0 = (N-(n**2).sum()/N)/(m-1); pi = x/n; msb = (n*(pi-ph)**2).sum()/(m-1); msw = (n*pi*(1-pi)).sum()/(N-m)
    if msw <= 0: iccu = 0.0 if msb <= 0 else 1.0
    else:
        F = msb/msw; flo = fdist.ppf(delta, m-1, N-m); th = (F/flo-1)/n0; iccu = 0.0 if th <= 0 else th/(1+th)
    nbw = (n**2).sum()/N; neff = N/(1+(nbw-1)*iccu); xeff = X*neff/N
    return float(betad.ppf(1-delta, xeff+1, neff-xeff))


def check_b1u_copy(mh_rev1_py):
    grab = lambda L: L[next(i for i, l in enumerate(L) if l.startswith("def b1u(")):][:8]
    assert grab(Path(mh_rev1_py).read_text().splitlines()) == grab(Path(__file__).read_text().splitlines()), "b1u copy differs"


def cstat(v, kind, gm):
    if kind == "c0": return 0.0
    w = len(v)/(len(v)+20.0)
    if kind == "c1s": return float(gm+w*(np.median(v)-gm))
    return float(gm+w*(np.median(np.sort(v)[:int(np.ceil(0.7*len(v)))])-gm))


class Inputs:
    """Pinned inputs; every score file is checked against its manifest when loaded."""
    def __init__(self, a2link, ids_iter1, mhlh, rev1_bbbc, rev1_mh, rev1_lh, mh_rev1_py):
        self.a2link, self.ids_iter1, self.mhlh = Path(a2link), Path(ids_iter1), Path(mhlh)
        self.rev1 = {"bbbc041_all": Path(rev1_bbbc), "mpidb_full": Path(rev1_mh), "lacuna_B": Path(rev1_lh)}
        for p in (a2link, ids_iter1, mhlh, rev1_bbbc, rev1_mh, rev1_lh):
            assert "iml" not in str(p).lower() and "plasmo" not in str(p).lower(), "confirmatory-site path refused"
        for k, p in (("manifest", self.a2link/"MANIFEST.tsv"), ("sums", self.mhlh/"SHA256SUMS"), ("ids_iter1", self.ids_iter1), ("rev1_bbbc", rev1_bbbc),
                     ("rev1_mh", rev1_mh), ("rev1_lh", rev1_lh), ("mh_rev1_py", mh_rev1_py)):
            assert sha(p) == PIN[k], ("pinned input differs", k)
        check_b1u_copy(mh_rev1_py); self.used = {}

    def _manifest_sha(self, sd):
        rows = [l.split("\t") for l in (self.a2link/"MANIFEST.tsv").read_text().splitlines()[1:] if l.strip()]
        hit = [r for r in rows if r[0].endswith(f"a2_s{sd}_scores.npz[bbbc041_all_logits]")]; assert len(hit) == 1, ("manifest row", sd); return hit[0][3]

    def _sums_sha(self, name):
        hit = [l.split()[0] for l in (self.mhlh/"SHA256SUMS").read_text().splitlines() if l.strip() and l.split()[1].lstrip("./").lstrip("*") == name]
        assert len(hit) == 1, ("SHA256SUMS line", name); return hit[0]

    def load_site(self, sd, site):
        assert site in DEV, ("site outside the development sites", site)
        if site == "bbbc041_all":
            f = self.a2link/f"a2_s{sd}_bbbc041_all_logits.npz"; h = sha(f); assert h == self._manifest_sha(sd), ("sha", f); self.used[str(f)] = h
            z = np.load(f)["logits"]; ids = np.load(self.ids_iter1, allow_pickle=True); y, g = np.asarray(ids["bbbc041_all_y"]), np.asarray(ids["bbbc041_all_group"])
        else:
            sh = SHORT[site]; f = self.mhlh/f"a2_s{sd}_{sh}_logits.npz"; fi = self.mhlh/f"ids_{sh}.npz"
            for p in (f, fi):
                h = sha(p); assert h == self._sums_sha(p.name), ("sha", p); self.used[str(p)] = h
            z = np.load(f)["logits"]; ids = np.load(fi, allow_pickle=True); y, g = np.asarray(ids["y"]), np.asarray(ids["field"])
        s = (z[:, 1]-z[:, 0]).astype(float); assert len(s) == len(y) == len(g), ("length", site, sd)
        return s, y, g


def population(inp, sd, site):
    """Site population for one seed: centred scores per candidate, sorted negatives / positives per candidate, oracles (as d_icc0.py / mh_rev1.py)."""
    s, y, g = inp.load_site(sd, site); fl = list(np.unique(g)); gm = float(np.median(s)); cs = {f: {k: cstat(s[g == f], k, gm) for k in CEN} for f in fl}
    sc = {k: s-np.array([cs[f][k] for f in g]) for k in CEN}
    fl_neg = {k: {f: np.sort(sc[k][(g == f) & (y == 0)]) for f in fl} for k in CEN}; fl_posv = {k: {f: sc[k][(g == f) & (y == 1)] for f in fl} for k in CEN}
    fl = [f for f in fl if len(fl_neg["c0"][f]) > 0]
    sneg = {k: np.sort(sc[k][y == 0]) for k in CEN}; spos = {k: np.sort(sc[k][y == 1]) for k in CEN}; N = len(sneg["c0"])
    sens_at = lambda k, t: 1-np.searchsorted(spos[k], t, side="left")/len(spos[k]); fpr_at = lambda k, t: 1-np.searchsorted(sneg[k], t, side="left")/N
    orac_unc = {e: sens_at("c0", np.quantile(sneg["c0"], 1-e)) for e in (0.10, 0.05)}
    orac_own = {(k, e): sens_at(k, np.quantile(sneg[k], 1-e)) for k in CEN for e in (0.10, 0.05)}
    return dict(s=s, y=y, g=g, sc=sc, fl=fl, fl_neg=fl_neg, fl_posv=fl_posv, sneg=sneg, spos=spos, N=N, sens_at=sens_at, fpr_at=fpr_at,
                orac_unc=orac_unc, orac_own=orac_own)


def fixed_sequence(grid, ub, eps):
    """REV-1 tau_of: descending grid, step down while the bound is <= eps; exception -> j = 0 (not certified)."""
    j = 0
    try:
        while j < len(grid) and ub(j) <= eps: j += 1
    except Exception: j = 0
    return (grid[j-1], True, j-1) if j > 0 else (grid[0], False, -1)


def replay(inp, sd, site, extra=None):
    """Replay the REV-1 draws for one seed and site.  Returns (rows, pop).  rows: one per draw x eps x rule with fpr, sens, cand, cert, dse, dso.
    Rules always computed: S3 (U3 = b1u at delta/3 on c0, c1s, c2s, union, labelled choice), B1u_c0 and B1u_c1s (U1 = b1u at delta 0.10, one candidate).
    extra: optional {name: (candidate, bound_fn(x, ni, j_context) -> upper bound)} evaluated with the same fixed-sequence test (used by W3)."""
    pop = population(inp, sd, site); fl, fl_neg, fl_posv = pop["fl"], pop["fl_neg"], pop["fl_posv"]; sens_at, fpr_at = pop["sens_at"], pop["fpr_at"]
    key8, sidx = SITECFG[site]; ai = AI; rows = []; m = M
    rng = np.random.default_rng([20261022, key8, ai, sd, m, sidx])
    for d_ in range(NDRAW):
        sel_ = rng.choice(len(fl), m, replace=False); F = [fl[i] for i in sel_]; clicks = sum(len(fl_posv["c0"][f]) for f in F); ni = np.array([len(fl_neg["c0"][f]) for f in F])
        for eps, dl in SETTINGS:   # kept: consumes the stream exactly as REV-1 (R1-2 subsets)
            r_ = np.random.default_rng(int(rng.integers(1 << 30))); [r_.choice(n_, 2, replace=False) if n_ >= 2 else np.arange(n_) for n_ in ni]
        D = {}
        for k in CEN:
            negs = [fl_neg[k][f] for f in F]; pooled = np.concatenate(negs); grid = np.unique(np.quantile(pooled, 1-np.arange(G)/G))[::-1]
            X = np.stack([ni[i]-np.searchsorted(negs[i], grid, side="left") for i in range(m)]).astype(float)
            posl = np.concatenate([fl_posv[k][f] for f in F]) if clicks > 0 else np.array([]); D[k] = (negs, pooled, grid, X, posl)
        cache = {}

        def ubf(rule, k, j):
            if (rule, k, j) not in cache: cache[(rule, k, j)] = b1u(D[k][3][:, j], ni, 0.10 if rule == "U1" else 0.10/3)
            return cache[(rule, k, j)]

        def rec(name, eps, k, t, cert):
            sn = sens_at(k, t)
            rows.append(dict(site=site, seed=sd, draw=d_, eps=eps, rule=name, cand=k, cert=bool(cert), fpr=float(fpr_at(k, t)), sens=float(sn),
                             dse=float(sn-pop["orac_unc"][eps]), dso=float(sn-pop["orac_own"][(k, eps)]), ni=ni.tolist(), fields=[str(f) for f in F]))
        for eps, dl in SETTINGS:
            T3 = {k: fixed_sequence(D[k][2], lambda j, k=k: ubf("U3", k, j), eps) for k in CEN}
            cert_c = [k for k in CEN if T3[k][1]]
            if not cert_c: rec("S3", eps, "c0", T3["c0"][0], False)
            else:
                if len(cert_c) == 1 or clicks == 0: ch = cert_c[0]
                else:
                    sl = [(D[k][4] >= T3[k][0]).mean() for k in cert_c]; mxs = max(sl); ch = cert_c[[i for i, v in enumerate(sl) if v == mxs][0]]
                rec("S3", eps, ch, T3[ch][0], True)
            for k in ("c0", "c1s"):
                t_, c_, _ = fixed_sequence(D[k][2], lambda j, k=k: ubf("U1", k, j), eps); rec("B1u_" + k, eps, k, t_, c_)
            for name, (k, fn) in (extra or {}).items():
                t_, c_, _ = fixed_sequence(D[k][2], lambda j, k=k, fn=fn: fn(D[k][3][:, j], ni, (site, sd, d_, eps, j)), eps); rec(name, eps, k, t_, c_)
    return rows, pop


def self_check(inp, rows):
    """S3, B1u_c0, B1u_c1s vs the stored REV-1 aggregates (1e-12).  Returns the list of mismatches (empty = passed)."""
    bad = []
    for site in DEV:
        agg = json.load(open(inp.rev1[site]))["agg"]
        for eps, _ in SETTINGS:
            for rule in ("S3", "B1u_c0", "B1u_c1s"):
                R = [r for r in rows if r["site"] == site and r["eps"] == eps and r["rule"] == rule]; st = agg[STORED_KEY[site].format(rule=rule, eps=eps)]
                fpr = np.array([r["fpr"] for r in R])
                mine = {"cov": float((fpr <= eps).mean()), "fpr_mean": float(fpr.mean()), "cannot_certify": float(1-np.mean([r["cert"] for r in R])),
                        "ds_eps": float(np.mean([r["dse"] for r in R])), "ds_own": float(np.mean([r["dso"] for r in R])), "n": len(R)}
                for kk, v in mine.items():
                    if abs(v-st[kk]) > 1e-12: bad.append((site, eps, rule, kk, v, st[kk]))
    return bad

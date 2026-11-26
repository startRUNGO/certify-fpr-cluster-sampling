# P5 stress tests of K1+S3 before the W1 pre-registration (实验管理 2026-10-30).  LOCAL, read-only.  WRITTEN, NOT RUN, until approved.
# NATURE: report only, no judgement line.  Any cell with coverage < 0.90 is listed in "risk_cells" and must go into the risk section of the W1 pre-registration.
# Uses the APPROVED functions unchanged: m3_replay_lib (ed9d1f17) and m3_k12_dev_eval.run_site_seed / summarise (0091a333), with K_ARMS set to ("K1",) in this
#   process only.  Arms reported: S3, B1u_c1s (= single-candidate c1s), K1 (= K1+S3).
# (b) m sensitivity: the three development sites, m in {10, 15, 22}.  REV-1 seeds its stream with m, so m = 10 and 15 are different random streams from m = 22.
#     Self-check (1e-12, before anything is written) wherever a stored REV-1 aggregate exists for (site, m): bbbc041_all at m = 10 and 22 (rev1_results.json),
#     mpidb_full and lacuna_B at m = 22; other (site, m) are labelled "new random stream, no stored aggregate".
# (a) old MP-IDB (development-site earlier version; S3 chose c0 in 60.3% of draws there): m = 22, REV-1 stream of rev1.py (site index 1: key8 = 8, sidx = 1).
#     Scores: seed 41 from the local iter1 file (sha must equal the a2_link MANIFEST source_sha256 of a2_s41_scores.npz), ONLY key "mpidb_full_logits" accessed;
#     seeds 42-45 from the VM extraction directory (each file's sha must match its SHA256SUMS line; key "logits").  Labels and fields from ids_iter1.npz keys
#     mpidb_full_y and mpidb_full_group only.  Self-check against rev1_results.json "a2|mpidb_full|22|..." (S3, B1u_c0, B1u_c1s) when all 5 seeds are present; with
#     seed 41 only, the check is not possible and the output says so (reference only).
# No IML, no PlasmoCount, no NIH key.  Output: one new JSON in a new directory ("x", allow_nan=False).  Prints only the output path.
# Usage: python -I -B m3_stress.py --a2link .. --ids_iter1 .. --mhlh .. --rev1_bbbc .. --rev1_mh .. --rev1_lh .. --mh_rev1_py .. --old_s41 <iter1 a2_s41_scores.npz>
#        [--old_dir <runs/a2_iter1_mpidb_old_<UTC>>] --out <dir/json>
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

K2 = Path(__file__).resolve().parents[1] / "k2prep"; sys.path.insert(0, str(K2))
import m3_replay_lib as L   # noqa: E402
import m3_k12_dev_eval as E   # noqa: E402

OLD = "mpidb_old"; MS = (10, 15, 22); ARMS = ("S3", "B1u_c1s", "K1")
PATTERN = {"bbbc041_all": "a2|bbbc041_all|{m}|{rule}|{eps}", "mpidb_full": "a2|mh|{m}|{rule}|{eps}", "lacuna_B": "a2|lh|{m}|{rule}|{eps}",
           OLD: "a2|mpidb_full|{m}|{rule}|{eps}"}


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


class StressInputs(L.Inputs):
    def __init__(self, a):
        super().__init__(a.a2link, a.ids_iter1, a.mhlh, a.rev1_bbbc, a.rev1_mh, a.rev1_lh, a.mh_rev1_py)
        self.rev1[OLD] = Path(a.rev1_bbbc); self.old = {}
        for p in [a.old_s41] + ([a.old_dir] if a.old_dir else []):
            assert "iml" not in str(p).lower() and "plasmo" not in str(p).lower(), "confirmatory-site path refused"
        h41 = sha(a.old_s41); assert h41 == self._manifest_src_sha(41), ("sha", a.old_s41); self.old[41] = (Path(a.old_s41), "mpidb_full_logits"); self.used[str(a.old_s41)] = h41
        if a.old_dir:
            d = Path(a.old_dir); lines = {l.split()[1].lstrip("./"): l.split()[0] for l in (d / "SHA256SUMS").read_text().splitlines() if l.strip()}
            for sd in (42, 43, 44, 45):
                f = d / f"a2_s{sd}_mpidb_full_logits.npz"; h = sha(f); assert h == lines.get(f.name), ("sha", f); self.old[sd] = (f, "logits"); self.used[str(f)] = h

    def _manifest_src_sha(self, sd):
        rows = [l.split("\t") for l in (self.a2link/"MANIFEST.tsv").read_text().splitlines()[1:] if l.strip()]
        hit = [r for r in rows if r[0].endswith(f"a2_s{sd}_scores.npz[bbbc041_all_logits]")]; assert len(hit) == 1, ("manifest row", sd); return hit[0][1]

    def load_site(self, sd, site):
        if site != OLD: return super().load_site(sd, site)
        f, key = self.old[sd]; z = np.load(f)[key]                           # one key only
        ids = np.load(self.ids_iter1, allow_pickle=True); y, g = np.asarray(ids["mpidb_full_y"]), np.asarray(ids["mpidb_full_group"])
        s = (z[:, 1]-z[:, 0]).astype(float); assert len(s) == len(y) == len(g), ("length", site, sd); return s, y, g


def check(inp, rows, sites, m):
    status, bad = {}, []
    for site in sites:
        agg = json.load(open(inp.rev1[site]))["agg"]; keys = {(eps, rule): PATTERN[site].format(m=m, rule=rule, eps=eps) for eps, _ in L.SETTINGS for rule in ("S3", "B1u_c0", "B1u_c1s")}
        if not all(k in agg for k in keys.values()) or len(L.SEEDS) < 5:
            status[site] = "new random stream: no stored aggregate to check (reference only)" if len(L.SEEDS) == 5 else f"fewer than 5 seeds ({len(L.SEEDS)}): self-check not possible (reference only)"
            continue
        for (eps, rule), k in keys.items():
            R = [r for r in rows if r["site"] == site and abs(r["eps"] - eps) < 1e-12 and r["rule"] == rule]; st = agg[k]; fpr = np.array([r["fpr"] for r in R])
            mine = {"cov": float((fpr <= eps).mean()), "fpr_mean": float(fpr.mean()), "cannot_certify": float(1-np.mean([r["cert"] for r in R])),
                    "ds_eps": float(np.mean([r["dse"] for r in R])), "ds_own": float(np.mean([r["dso"] for r in R])), "n": len(R)}
            bad += [(site, m, eps, rule, kk, v, st[kk]) for kk, v in mine.items() if abs(v - st[kk]) > 1e-12]
        status[site] = "self-check passed (REV-1 stream, 1e-12)"
    return status, bad


def run_block(inp, sites, m, test):
    L.M = m; E.EXC.clear(); E.DEG.clear(); E.DROP.clear(); rows = []
    for site in sites:
        for sd in L.SEEDS: rows += E.run_site_seed(inp, sd, site)
    status, bad = check(inp, rows, sites, m)
    if bad: raise SystemExit("self-check failed, nothing written: " + repr(bad[:5]))
    rng = np.random.default_rng(E.CI_SEED)
    cells = [dict(test=test, site=s, m=m, arm=a, eps=e, stream=status[s], **E.summarise([r for r in rows if r["site"] == s and abs(r["eps"] - e) < 1e-12 and r["rule"] == a], e, s, a, rng))
             for s in sites for a in ARMS for e, _ in L.SETTINGS]
    return cells


def main():
    ap = argparse.ArgumentParser()
    for k in ("a2link", "ids_iter1", "mhlh", "rev1_bbbc", "rev1_mh", "rev1_lh", "mh_rev1_py", "old_s41", "out"): ap.add_argument("--" + k, required=True)
    ap.add_argument("--old_dir", default=None)
    a = ap.parse_args(); out = Path(a.out); assert not out.exists(), "output exists"
    inp = StressInputs(a); E.K_ARMS = ("K1",); seeds_all = list(L.SEEDS); cells = []
    for m in MS: cells += run_block(inp, L.DEV, m, "m_sensitivity")
    L.SITECFG[OLD] = (8, 1); L.DEV = tuple(L.DEV) + (OLD,); L.SEEDS = [s for s in seeds_all if s in inp.old]
    cells += run_block(inp, (OLD,), 22, "old_mpidb"); L.SEEDS = seeds_all
    risk = [dict(test=c["test"], site=c["site"], m=c["m"], arm=c["arm"], eps=c["eps"], coverage=c["coverage"]) for c in cells if c["coverage"] < 0.90]
    rep = dict(script_sha256=sha(__file__), lib_sha256=sha(L.__file__), eval_sha256=sha(E.__file__), inputs_sha256=dict(L.PIN), score_files_sha256=inp.used,
               nature="stress tests, report only; cells with coverage < 0.90 go into the W1 pre-registration risk section",
               old_mpidb_seeds=sorted(inp.old), cells=cells, risk_cells=risk)
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False)
    out.parent.mkdir(parents=False, exist_ok=False)
    with open(out, "x") as fh: fh.write(txt)
    print("written", out)


if __name__ == "__main__":
    main()

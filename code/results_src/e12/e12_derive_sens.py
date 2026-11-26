# P5a E12 derive: sensitivity, Delta Sens_eps, cannot-certify and their Monte Carlo SE for S3 and SPLIT6 in the 24 E12 cells, from the stored per-draw rows (实验管理, 2026-11-24).
# WRITTEN, NOT RUN, until 实验管理 reviews it.  READ-ONLY on its inputs; no simulation is run, no rule or constant is changed; every cell is reported.
# Inputs (pinned, exit 3 if a sha differs): e12_rows.json (ece7c768, runs/e11_e12_20261108T121449Z/out_e12) and e12_summary.json (f3c40a94, same directory).
# Per cell (0..23), rule (S3, SPLIT6), eps (0.10, 0.05) and realisation (seeds 41..45): n, coverage (fpr <= eps), cannot_certify (1 - mean cert), sens, ds_eps (dse), ds_own (dso), fpr_mean, as means over the
#   200 draws of the realisation; then the mean over realisations, the 5 per-realisation values and MCSE = SD (ddof 1) / sqrt(5).  Draws without a certified threshold are evaluated at the fallback
#   threshold, as in the stored rows and the manuscript (cannot-certify is reported beside them).
# Anchor (exit 3, nothing written): coverage_per_realisation, cannot_certify and ds_own recomputed here equal e12_summary.json within 1e-12 for every cell x rule x eps; every (cell, rule, eps, seed) has 200 rows.
# Output: one NEW directory e12derive_<UTC yyyymmddThhmmssZ>: e12_derived_sens.json and SHA256SUMS ('x' mode, allow_nan=False), for the add-only master builder.
# Usage: python3 -I -B e12_derive_sens.py --rows <e12_rows.json> --summary <e12_summary.json> --out <parent>/e12derive_<UTC>
import argparse, hashlib, json, re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

PINS = {"rows": "ece7c768005f0a8b79e0e4a883a1a50dc9f6f70e5b70995e4ece4d9949d54a20", "summary": "f3c40a947e24a92bdf360b640f921c6fc57c9f1955d704f0a7a498806a313941"}
RULES = ("S3", "SPLIT6"); EPSS = (0.1, 0.05); SEEDS = (41, 42, 43, 44, 45); NDRAW = 200; N_CELLS = 24
METRICS = ("coverage", "cannot_certify", "sens", "ds_eps", "ds_own", "fpr_mean"); OUT_RE = re.compile(r"^e12derive_\d{8}T\d{6}Z$")


class AnchorError(Exception): pass


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--rows", required=True); ap.add_argument("--summary", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(argv); out = Path(a.out)
    if out.exists() or not OUT_RE.match(out.name) or not out.parent.is_dir(): raise SystemExit("--out must be a new directory e12derive_<UTC yyyymmddThhmmssZ> in an existing parent")
    for k, p in (("rows", a.rows), ("summary", a.summary)):
        if sha(p) != PINS[k]: raise AnchorError(("input differs from its pin", k))
    with open(a.summary) as fh: summ = json.load(fh)
    with open(a.rows) as fh: rows = json.load(fh)
    by = defaultdict(list)
    for r in rows:
        if r["rule"] in RULES: by[(r["cell"], r["rule"], round(r["eps"], 6), r["seed"])].append(r)
    cells = {}
    for ci in range(N_CELLS):
        res = {}
        for eps in EPSS:
            for rule in RULES:
                vals = {mn: [] for mn in METRICS}
                for sd in SEEDS:
                    R_ = by.get((ci, rule, round(eps, 6), sd), [])
                    if len(R_) != NDRAW: raise AnchorError(("row count", ci, rule, eps, sd, len(R_)))
                    fpr = np.array([r["fpr"] for r in R_])
                    vals["coverage"].append(float(np.mean([r["fpr"] <= eps for r in R_]))); vals["cannot_certify"].append(float(1 - np.mean([r["cert"] for r in R_])))
                    vals["sens"].append(float(np.mean([r["sens"] for r in R_]))); vals["ds_eps"].append(float(np.mean([r["dse"] for r in R_])))
                    vals["ds_own"].append(float(np.mean([r["dso"] for r in R_]))); vals["fpr_mean"].append(float(fpr.mean()))
                ref = summ["cells"][str(ci)]["results"][f"{rule}|{eps}"]
                ok = all(abs(u - v) <= 1e-12 for u, v in zip(vals["coverage"], ref["coverage_per_realisation"])) and len(vals["coverage"]) == len(ref["coverage_per_realisation"])
                ok = ok and abs(float(np.mean(vals["cannot_certify"])) - ref["cannot_certify"]) <= 1e-12 and abs(float(np.mean(vals["ds_own"])) - ref["ds_own"]) <= 1e-12
                if not ok: raise AnchorError(("recomputed values differ from e12_summary.json", ci, rule, eps))
                res[f"{rule}|{eps}"] = dict(n=NDRAW * len(SEEDS), **{mn: dict(mean=float(np.mean(v)), per_realisation=v, mcse=float(np.std(v, ddof=1) / np.sqrt(len(v)))) for mn, v in vals.items()})
        cells[str(ci)] = dict(name=summ["cells"][str(ci)]["cell"]["name"], results=res)
    rep = dict(version="e12derive_r1", request="实验管理 2026-11-24: E12 sens, Delta Sens_eps, cannot-certify and MCSE for S3 and SPLIT6", script_sha256=sha(__file__), input_sha256=PINS,
               anchor="passed (coverage per realisation, cannot_certify and ds_own equal e12_summary.json within 1e-12, 24 cells x 2 rules x 2 eps; 200 rows per cell x rule x eps x seed)",
               mcse_definition="SD (ddof 1) of the 5 per-realisation means / sqrt(5)", cells=cells, note="derived from stored rows; nothing simulated; fallback-threshold convention as the stored rows")
    txt = json.dumps(rep, indent=1, ensure_ascii=False, allow_nan=False)
    out.mkdir(parents=False, exist_ok=False)
    with open(out / "e12_derived_sens.json", "x") as fh: fh.write(txt)
    with open(out / "SHA256SUMS", "x") as fh: fh.write(f"{sha(out / 'e12_derived_sens.json')}  e12_derived_sens.json\n")
    print(out); return 0


if __name__ == "__main__":
    try: sys.exit(main())
    except AnchorError as e: print("ANCHOR MISMATCH (exit 3, nothing written):", e, file=sys.stderr); sys.exit(3)

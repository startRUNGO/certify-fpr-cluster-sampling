# P5a derivation D-c (实验管理 2026-11-25, F-c option (i)): draw-level realized FPR of the simulation study pooled over the 50 synthetic sites, for the box plots of F-c,
#   and the count of uncertified draws with realized FPR above eps for every simulation cell (decides the cannot-certify clause of the F-b / F-c captions).
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  LOCAL, synthetic data only, no VM.  Reads pinned inputs, writes one new directory ('x' mode).
# Settings of F-c, fixed by design before any number is read (实验管理 2026-11-25):
#   e10|0 base; e10|3 between-field shift xi = 1.0 (high); e10|9 heavy tail t5; e12|19 two-stage slides, S = 12, omega = 0.8, s = 2, chosen as the stronger slide stress
#   by design (omega = 0.8 puts the larger share of between-field variance between slides; e12|17 has omega = 0.5).  Rules: SPLIT6, S3, CB_c1s, GEE_c1s, BB_c1s; eps 0.1 and 0.05.
# Row sources (as sim50_extend.py ed74e23d, sources_orig5 of every sim50 summary):
#   seeds 46..90, every cell and rule: sim50_rows.npz of the three batches (pinned);
#   seeds 41..45, E10 S3 / SPLIT6: the E10 replay rows in the batch-1 npz (anchor b);
#   seeds 41..45, E12 S3 / SPLIT6: e12_rows.json (ece7c768, pulled 20261125T190947Z, rows with rule S3 / SPLIT6, field 'cell');
#   seeds 41..45, CB_c1s / GEE_c1s / BB_c1s: simw3_rows.json (bed0438a, simw3_20261124T194046Z).
# Draw semantics (code, quoted to 实验管理 2026-11-25): an uncertified draw carries the realized site FPR at grid[0] (tau_1, the largest grid threshold) of its candidate
#   (m3_replay_lib.fixed_sequence, line 105; e12_stress.py line 130; simw3_baselines.py lines 141 / 148 / 150).  Coverage = mean(fpr <= eps) over all draws.
# Gates (exit 3, nothing written):
#   (1) sha256 of every input == PINS (three summaries, three npz, e12_rows.json, simw3_rows.json, results_master.json 9b04913e);
#   (2) for every (design, cell, rule, eps) of the 340: rows for exactly the seeds 41..90, each with exactly DRAWS draws, draw ids 0..DRAWS-1 without repeats;
#   (3) per-seed coverage mean(fpr <= eps) and cannot-certify 1 - mean(cert) equal the all50 per_realisation values of the summary that holds the result key
#       (seed order = its all50 'seeds' list) within 1e-12, and the coverage values equal master sim50/<rk>|all50|coverage|per_realisation|<i> within 1e-12;
#   (4) every fpr finite and in [0, 1]; cert boolean.
# Output per (design, cell, rule, eps), all 340: n, n_cert, n_uncert, n_uncert_fpr_gt_eps.
#   For the F-c settings only (4 x 5 x 2 = 40): among certified draws q05, q25, q50, q75, q95 (numpy.percentile, method 'linear'), mean, min, max; and when
#   n_cert < STRIP_N (100) the sorted certified values themselves (the figure then draws a strip with n printed, decided by n_cert, not by eye).
# Output (new directory, 'x' mode): d_c_simbox.json (allow_nan=False), SHA256SUMS.
# Usage: python3 -E -B d_c_simbox.py --scratch <scratchpad dir holding the input dirs> --master <tree/results_master.json> --out <new dir>
import argparse, hashlib, json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

BATCHES = ["fifty_20261124T214602Z", "fifty_20261125T094940Z", "fifty_20261125T141231Z"]
PINS = {
    "fifty_20261124T214602Z/sim50_summary.json": "7bb098ca3561a3216a39274bbce1dc88b8a0e96847fc5a1d898168bcbc2c29ab",
    "fifty_20261125T094940Z/sim50_summary.json": "fa88e4f0102b1d90062a96164772c3b0cdcec43cf614dd982df8ac0e5b549d0f",
    "fifty_20261125T141231Z/sim50_summary.json": "115968f1a0ce0b30c28cfa2eadc31d28c8c5cf5eeac9622c59ccfe647dd2f66e",
    "fifty_20261124T214602Z/sim50_rows.npz": "88baad0d6bc01063548bbe1c996dbce6434548cd45f47fe9705afbcf81fdb95b",
    "fifty_20261125T094940Z/sim50_rows.npz": "6a04d4b9ad59f018fe758b581784f8466dfff3bb886593543935ee7f54e56306",
    "fifty_20261125T141231Z/sim50_rows.npz": "c82330c224438bc041d6fac94c5a582aaf7525f023af2b5d716d33a1ed4910cf",
    "e12rows_pull_20261125T190947Z/e12_rows.json": "ece7c768005f0a8b79e0e4a883a1a50dc9f6f70e5b70995e4ece4d9949d54a20",
    "simw3_20261124T194046Z/simw3_rows.json": "bed0438a1f0717d98482502a0baa93d2fddfefb88e3179048e44ea8bff2d3794"}
MASTER_SHA = "9b04913eba0b2044e12fce07b935de17a7e5065fb616ee636dad9c54fb6fec4a"
RULES = ["SPLIT6", "S3", "CB_c1s", "GEE_c1s", "BB_c1s"]
EPS = [0.1, 0.05]
CELLS = [("e10", i) for i in range(10)] + [("e12", i) for i in range(24)]
FC_SETTINGS = [("e10", 0), ("e10", 3), ("e10", 9), ("e12", 19)]
SEEDS = list(range(41, 91))
DRAWS = 200
STRIP_N = 100
QS = (5, 25, 50, 75, 95)
TOL = 1e-12


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def es(e): return "0.1" if abs(e - 0.1) < 1e-12 else "0.05" if abs(e - 0.05) < 1e-12 else None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--scratch", required=True); ap.add_argument("--master", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); S = Path(a.scratch); out = Path(a.out)
    if out.exists(): return fail("output directory exists", out)
    for rel, want in PINS.items():
        if not (S / rel).is_file(): return fail("missing input", rel)
        if sha(S / rel) != want: return fail("input sha differs from the pin", rel)
    if sha(a.master) != MASTER_SHA: return fail("master sha differs from the pin")
    V = json.load(open(a.master))["values"]
    # rows[(design, cell, rule, eps_str)][seed] -> list of (draw, cert, fpr)
    rows = defaultdict(lambda: defaultdict(list))

    def add(d, c, r, e, sd, dr, ct, fp):
        k = es(float(e))
        if k is None or r not in RULES: return
        rows[(d, int(c), r, k)][int(sd)].append((int(dr), bool(ct), float(fp)))
    for b in BATCHES:
        z = np.load(S / b / "sim50_rows.npz"); dn, rn = list(z["design_names"]), list(z["rule_names"])
        for d_, c_, r_, e_, sd, dr, ct, fp in zip(z["design"], z["cell"], z["rule"], z["eps"], z["seed"], z["draw"], z["cert"], z["fpr"]):
            d, r = str(dn[d_]), str(rn[r_])
            # orig seeds from the npz only for the E10 S3 / SPLIT6 replay (anchor b); every other orig row comes from its own source below
            if sd <= 45 and not (d == "e10" and r in ("S3", "SPLIT6")): return fail("unexpected seed 41..45 row in", b, d, r)
            add(d, c_, r, e_, sd, dr, ct, fp)
    for x in json.load(open(S / "e12rows_pull_20261125T190947Z/e12_rows.json")):
        if x["rule"] in ("S3", "SPLIT6"): add("e12", x["cell"], x["rule"], x["eps"], x["seed"], x["draw"], x["cert"], x["fpr"])
    for x in json.load(open(S / "simw3_20261124T194046Z/simw3_rows.json")):
        if x["rule"] in ("CB_c1s", "GEE_c1s", "BB_c1s"): add(x["design"], x["cell"], x["rule"], x["eps"], x["seed"], x["draw"], x["cert"], x["fpr"])
    summ = {}
    for b in BATCHES:
        for rk, v in json.load(open(S / b / "sim50_summary.json"))["results"].items():
            if rk in summ: return fail("result key in two summaries", rk)
            summ[rk] = v["all50"]
    res = {}; fc = {}
    for d, c in CELLS:
        for r in RULES:
            for e in EPS:
                k = (d, c, r, es(e)); rk = f"{d}|{c}|{r}|{es(e)}"
                if rk not in summ: return fail("result key missing from the summaries", rk)
                per = rows.get(k)
                if per is None or sorted(per) != SEEDS: return fail("seeds differ from 41..90", rk, sorted(per or []))
                A = summ[rk]
                if A["seeds"] != SEEDS: return fail("summary seed order differs from 41..90", rk)
                cert_all, fpr_all = [], []
                for i, sd in enumerate(SEEDS):
                    rr = per[sd]
                    if sorted(x[0] for x in rr) != list(range(DRAWS)): return fail("draw ids differ from 0..DRAWS-1", rk, sd)
                    ct = np.array([x[1] for x in rr], bool); fp = np.array([x[2] for x in rr], float)
                    if not np.all(np.isfinite(fp)) or fp.min() < 0 or fp.max() > 1: return fail("fpr not finite or outside [0, 1]", rk, sd)
                    cov, cc = float(np.mean(fp <= e)), float(1 - np.mean(ct))
                    if abs(cov - A["coverage"]["per_realisation"][i]) > TOL: return fail("per-seed coverage differs from the summary", rk, sd)
                    if abs(cc - A["cannot_certify"]["per_realisation"][i]) > TOL: return fail("per-seed cannot-certify differs from the summary", rk, sd)
                    mk = f"sim50/{rk}|all50|coverage|per_realisation|{i}"
                    if mk not in V or abs(cov - V[mk]) > TOL: return fail("per-seed coverage differs from the master", mk)
                    cert_all.append(ct); fpr_all.append(fp)
                ct, fp = np.concatenate(cert_all), np.concatenate(fpr_all)
                res[rk] = dict(n=int(ct.size), n_cert=int(ct.sum()), n_uncert=int((~ct).sum()), n_uncert_fpr_gt_eps=int(((~ct) & (fp > e)).sum()))
                if (d, c) in FC_SETTINGS:
                    v = np.sort(fp[ct]); o = dict(n=int(ct.size), n_cert=int(v.size))
                    if v.size:
                        o.update({f"q{q:02d}": float(np.percentile(v, q, method="linear")) for q in QS}, mean=float(v.mean()), min=float(v[0]), max=float(v[-1]))
                    if v.size < STRIP_N: o["strip_values"] = [float(x) for x in v]
                    fc[rk] = o
    tot = {e: sum(x["n_uncert_fpr_gt_eps"] for rk, x in res.items() if rk.endswith("|" + e)) for e in ("0.1", "0.05")}
    for o in list(res.values()) + list(fc.values()):
        for x in o.values():
            if isinstance(x, float) and not math.isfinite(x): return fail("non-finite output value")
    payload = json.dumps(dict(script_sha256=sha(__file__), inputs_sha256=PINS, master_sha256=MASTER_SHA, seeds=SEEDS, draws_per_realisation=DRAWS,
                              fc_settings=[f"{d}|{c}" for d, c in FC_SETTINGS], strip_n=STRIP_N, percentile_method="numpy.percentile linear, among certified draws",
                              uncert_semantics="an uncertified draw carries the realized site FPR at the largest grid threshold of its candidate (m3_replay_lib.fixed_sequence)",
                              counts=res, n_uncert_fpr_gt_eps_total=tot, fc=fc), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    with open(out / "d_c_simbox.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh: fh.write(f"{sha(out / 'd_c_simbox.json')}  d_c_simbox.json\n")
    print("written", out); return 0


if __name__ == "__main__": sys.exit(main())

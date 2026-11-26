# P5a ADD-ONLY results_master builder: derived counts, minima and Monte Carlo margins of the 50-site simulation, sim50_derived/* (实验管理 2026-11-25, OLD->NEW list B and D).
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  No input file: every new key is computed from sim50/* keys already in the master (add_sim50_keys.py 71af83c7).
# Same mechanics as add_simbox_keys.py (c19380c2): write targets checked absent first; indent=1; every old key, source and files_sha256 entry unchanged; key count checked;
#   every new key recomputed from the written file and compared.
# Coverage = sim50/<d>|<ci>|<rule>|<eps>|all50|coverage|mean, MCSE = ...|all50|coverage|mcse (SD of the 50 per-site coverages / sqrt(50)); delta = protocol/delta.
# Keys (rule in S3, SPLIT6):
#   sim50_derived/e12|<rule>|{n_settings, min_coverage, argmin, min_mcse, n_hold, n_between, n_fail, n_below_by_2mcse}
#     over the 24 stress cells x 2 eps; n_hold: coverage >= 1 - delta; n_between: 1 - 2 delta <= coverage < 1 - delta; n_fail: coverage < 1 - 2 delta;
#     n_below_by_2mcse: coverage + 2 MCSE < 1 - delta; argmin is the string '<d>|<ci>|<eps>'.
#   sim50_derived/e10|<rule>|<eps>|{n_cells, min_coverage, argmin, min_mcse, n_hold} over the 10 single-factor cells.
#   sim50_derived/shortfall|<d>|<ci>|<rule>|<eps>|{coverage, mcse, shortfall_in_mcse, below_by_2mcse} for every (d, ci, rule, eps) of S3 / SPLIT6 with coverage < 1 - delta;
#     shortfall_in_mcse = (1 - delta - coverage) / MCSE; below_by_2mcse = 1 if coverage + 2 MCSE < 1 - delta else 0 (stored as int).
# Sources: {'derived': <rule in words>, 'from_keys': <the sim50/* keys used>} (as sim50/platt_fallback|*).
# Gates (exit 3, nothing written): master sha dba06e5c; protocol/delta == 0.1; no existing key starts with sim50_derived/; every input key present and finite with
#   n_realisations == 50; n_hold + n_between + n_fail == n_settings == 48; every new key new; all write targets absent.
# Backup: results_master_before_sim50_derived_dba06e5c.json ('x').
# Usage: python3 -E -B add_sim50_derived_keys.py --tree <tree_bj_v1_...>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "dba06e5cf5c9f73b432eb2717bc542b990bb473ed16370f869207c1e6b0e12d1"
PRE = "sim50_derived"
RULES = ("S3", "SPLIT6")
EPS = ("0.1", "0.05")
E12, E10 = range(24), range(10)


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def derive(V):
    """Returns {key: (value, source)} computed from V, or raises ValueError."""
    delta = V["protocol/delta"]
    if delta != 0.1: raise ValueError("protocol/delta differs from 0.1")
    lo1, lo2 = 1 - delta, 1 - 2 * delta

    def cov(d, ci, r, e):
        b = f"sim50/{d}|{ci}|{r}|{e}|all50|"
        ks = [b + "coverage|mean", b + "coverage|mcse", b + "n_realisations"]
        for k in ks:
            if k not in V or not isinstance(V[k], (int, float)) or not math.isfinite(V[k]): raise ValueError(f"missing or non-finite {k}")
        if V[ks[2]] != 50: raise ValueError(f"n_realisations differs from 50: {b}")
        return V[ks[0]], V[ks[1]], ks[:2]
    new = {}

    def put(k, v, words, used):
        k = f"{PRE}/{k}"
        if k in new: raise ValueError(f"duplicate {k}")
        new[k] = (v, dict(derived=words, from_keys=sorted(set(used))))
    for r in RULES:
        rows = [(cov("e12", ci, r, e), ci, e) for ci in E12 for e in EPS]
        used = [k for (c, m, ks), _, _ in rows for k in ks]
        (cmin, mmin, _), ci0, e0 = min(rows, key=lambda x: x[0][0])
        nh = sum(c >= lo1 for (c, _, _), _, _ in rows); nb = sum(lo2 <= c < lo1 for (c, _, _), _, _ in rows); nf = sum(c < lo2 for (c, _, _), _, _ in rows)
        n2 = sum(c + 2 * m < lo1 for (c, m, _), _, _ in rows)
        if nh + nb + nf != len(rows) or len(rows) != 48: raise ValueError(f"e12 {r}: counts do not add up")
        w = f"over the 24 E12 stress cells x eps 0.1, 0.05 of all50 coverage of rule {r}, delta = protocol/delta"
        put(f"e12|{r}|n_settings", len(rows), "number of settings " + w, used)
        put(f"e12|{r}|min_coverage", cmin, "minimum " + w, used)
        put(f"e12|{r}|argmin", f"e12|{ci0}|{e0}", "setting of the minimum " + w, used)
        put(f"e12|{r}|min_mcse", mmin, "all50 coverage MCSE at the minimum " + w, used)
        put(f"e12|{r}|n_hold", nh, "settings with coverage >= 1 - delta " + w, used)
        put(f"e12|{r}|n_between", nb, "settings with 1 - 2 delta <= coverage < 1 - delta " + w, used)
        put(f"e12|{r}|n_fail", nf, "settings with coverage < 1 - 2 delta " + w, used)
        put(f"e12|{r}|n_below_by_2mcse", n2, "settings with coverage + 2 MCSE < 1 - delta " + w, used)
        for e in EPS:
            rows10 = [(cov("e10", ci, r, e), ci) for ci in E10]; used10 = [k for (c, m, ks), _ in rows10 for k in ks]
            (c10, m10, _), ci1 = min(rows10, key=lambda x: x[0][0])
            w10 = f"over the 10 E10 single-factor cells of all50 coverage of rule {r} at eps {e}"
            put(f"e10|{r}|{e}|n_cells", len(rows10), "number of cells " + w10, used10)
            put(f"e10|{r}|{e}|min_coverage", c10, "minimum " + w10, used10)
            put(f"e10|{r}|{e}|argmin", f"e10|{ci1}|{e}", "cell of the minimum " + w10, used10)
            put(f"e10|{r}|{e}|min_mcse", m10, "all50 coverage MCSE at the minimum " + w10, used10)
            put(f"e10|{r}|{e}|n_hold", sum(c >= lo1 for (c, _, _), _ in rows10), "cells with coverage >= 1 - delta " + w10, used10)
        for d, cells in (("e12", E12), ("e10", E10)):
            for ci in cells:
                for e in EPS:
                    c, m, ks = cov(d, ci, r, e)
                    if c >= lo1: continue
                    if not m > 0: raise ValueError(f"MCSE not positive at a shortfall: {d}|{ci}|{r}|{e}")
                    b = f"shortfall|{d}|{ci}|{r}|{e}"; ww = f"all50 coverage of {d}|{ci}|{r}|{e} below 1 - delta"
                    put(b + "|coverage", c, ww + ": coverage", ks)
                    put(b + "|mcse", m, ww + ": MCSE", ks)
                    put(b + "|shortfall_in_mcse", (lo1 - c) / m, ww + ": (1 - delta - coverage) / MCSE", ks)
                    put(b + "|below_by_2mcse", int(c + 2 * m < lo1), ww + ": 1 if coverage + 2 MCSE < 1 - delta else 0", ks)
    return new


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); a = ap.parse_args()
    tree = Path(a.tree); rm = tree / "results_master.json"
    try:
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not dba06e5c")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        if any(k.startswith(PRE + "/") for k in vals): raise ValueError(f"existing keys already use prefix {PRE}/")
        new = derive(vals)
        for k in new:
            if k in vals: raise ValueError(f"key exists: {k}")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    bak = tree / "results_master_before_sim50_derived_dba06e5c.json"; tmp = tree / "results_master.json.sim50_derived_tmp"
    existing = [str(x) for x in (bak, tmp) if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}; old_fsha = dict(fsha)
    with open(rm, "rb") as s_, open(bak, "xb") as d_: d_.write(s_.read())
    for k, (v, s) in new.items(): vals[k] = v; srcs[k] = s
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp))
    try: again = derive({k: v for k, v in chk["values"].items() if not k.startswith(PRE + "/")})
    except ValueError as e: print("GATE (exit 3): recompute failed", e, file=sys.stderr); return 3
    for k, (v, _) in new.items():
        if chk["values"].get(k) != v or again[k][0] != v: print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    for k, v in old_fsha.items():
        if chk["files_sha256"].get(k) != v: print("GATE (exit 3): an old files_sha256 entry changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new): print("GATE (exit 3): key count mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    print(json.dumps(dict(added=len(new), n_values=d["n_values"], new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

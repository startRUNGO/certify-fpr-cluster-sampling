# P5a ADD-ONLY results_master builder r2: further derived keys of the 50-site simulation for the §4 ADEMP prose (实验管理 2026-11-25, outline decision 3).
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  No input file: every new key is computed from sim50/* and simbox/* keys already in the master.
# Same mechanics as add_sim50_derived_keys.py (350511e8): write targets checked absent first; indent=1; every old key, source and files_sha256 entry unchanged;
#   key count checked; every new key recomputed from the written file and compared.
# Coverage = sim50/<d>|<ci>|<rule>|<eps>|all50|coverage|mean, MCSE = ...|coverage|mcse; delta = protocol/delta.
# Keys (prefix sim50_derived/, none may exist before the run):
#   e12|<rule>|{n_settings, min_coverage, argmin, min_mcse, n_hold, n_between, n_fail, n_below_by_2mcse} and
#   e10|<rule>|<eps>|{n_cells, min_coverage, argmin, min_mcse, n_hold} for rule in CB_c1s, GEE_c1s, BB_c1s (definitions as in r1 for S3 / SPLIT6);
#   max_mcse|value, max_mcse|argmax ('<d>|<ci>|<rule>|<eps>') over the 340 settings (34 cells x 5 rules x 2 eps) of all50 coverage MCSE;
#   range|<d>|<rule>|<eps>|<stat>|{min, max, argmin, argmax} for d in e10, e12, the 5 rules, both eps, stat in cannot_certify, ds_eps (all50 ...|<stat>|mean over the
#     cells of the design; argmin / argmax is the cell '<d>|<ci>');
#   fc_q95_gt_eps|<d>|<ci>|<rule>|<eps> = 1 if simbox/fc|<d>|<ci>|<rule>|<eps>|q95 > eps else 0, for the 40 F-c entries.
# Sources: {'derived': <rule in words>, 'from_keys': <keys used>}.
# Gates (exit 3, nothing written): master sha effb684e; protocol/delta == 0.1; every new key absent; every input key present and finite, n_realisations == 50;
#   n_hold + n_between + n_fail == n_settings == 48; 340 MCSE values; 40 F-c entries with n_cert == n; all write targets absent.
# Backup: results_master_before_sim50_derived2_effb684e.json ('x').
# Usage: python3 -E -B add_sim50_derived2_keys.py --tree <tree_bj_v1_...>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "effb684e0b83f85dd52f45e0c50d77ceb9505ec9ba9364459fdd67af7dc62ce4"
PRE = "sim50_derived"
RULES_ALL = ("S3", "SPLIT6", "CB_c1s", "GEE_c1s", "BB_c1s")
RULES_NEW = ("CB_c1s", "GEE_c1s", "BB_c1s")
EPS = ("0.1", "0.05")
CELLS = {"e10": range(10), "e12": range(24)}
FC = (("e10", 0), ("e10", 3), ("e10", 9), ("e12", 19))


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def derive(V):
    """Returns {key: (value, source)} computed from V, or raises ValueError."""
    delta = V["protocol/delta"]
    if delta != 0.1: raise ValueError("protocol/delta differs from 0.1")
    lo1, lo2 = 1 - delta, 1 - 2 * delta

    def num(k):
        if k not in V or not isinstance(V[k], (int, float)) or isinstance(V[k], bool) or not math.isfinite(V[k]): raise ValueError(f"missing or non-finite {k}")
        return V[k]

    def cov(d, ci, r, e):
        b = f"sim50/{d}|{ci}|{r}|{e}|all50|"
        if num(b + "n_realisations") != 50: raise ValueError(f"n_realisations differs from 50: {b}")
        return num(b + "coverage|mean"), num(b + "coverage|mcse"), [b + "coverage|mean", b + "coverage|mcse"]
    new = {}

    def put(k, v, words, used):
        k = f"{PRE}/{k}"
        if k in new: raise ValueError(f"duplicate {k}")
        new[k] = (v, dict(derived=words, from_keys=sorted(set(used))))
    for r in RULES_NEW:
        rows = [(cov("e12", ci, r, e), ci, e) for ci in CELLS["e12"] for e in EPS]
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
            rows10 = [(cov("e10", ci, r, e), ci) for ci in CELLS["e10"]]; used10 = [k for (c, m, ks), _ in rows10 for k in ks]
            (c10, m10, _), ci1 = min(rows10, key=lambda x: x[0][0])
            w10 = f"over the 10 E10 single-factor cells of all50 coverage of rule {r} at eps {e}"
            put(f"e10|{r}|{e}|n_cells", len(rows10), "number of cells " + w10, used10)
            put(f"e10|{r}|{e}|min_coverage", c10, "minimum " + w10, used10)
            put(f"e10|{r}|{e}|argmin", f"e10|{ci1}|{e}", "cell of the minimum " + w10, used10)
            put(f"e10|{r}|{e}|min_mcse", m10, "all50 coverage MCSE at the minimum " + w10, used10)
            put(f"e10|{r}|{e}|n_hold", sum(c >= lo1 for (c, _, _), _ in rows10), "cells with coverage >= 1 - delta " + w10, used10)
    allm = [(cov(d, ci, r, e)[1], f"{d}|{ci}|{r}|{e}", cov(d, ci, r, e)[2][1]) for d, cs in CELLS.items() for ci in cs for r in RULES_ALL for e in EPS]
    if len(allm) != 340: raise ValueError(f"{len(allm)} MCSE values, expected 340")
    mx = max(allm, key=lambda x: x[0]); wm = "over the 340 settings (34 cells x 5 rules x 2 eps) of the all50 coverage MCSE"
    put("max_mcse|value", mx[0], "maximum " + wm, [k for _, _, k in allm])
    put("max_mcse|argmax", mx[1], "setting of the maximum " + wm, [k for _, _, k in allm])
    for d, cs in CELLS.items():
        for r in RULES_ALL:
            for e in EPS:
                for st in ("cannot_certify", "ds_eps"):
                    vals = [(num(f"sim50/{d}|{ci}|{r}|{e}|all50|{st}|mean"), ci) for ci in cs]; used = [f"sim50/{d}|{ci}|{r}|{e}|all50|{st}|mean" for ci in cs]
                    lo_, hi_ = min(vals), max(vals); w = f"over the {len(vals)} {d.upper()} cells of all50 {st} of rule {r} at eps {e}"
                    b = f"range|{d}|{r}|{e}|{st}|"
                    put(b + "min", lo_[0], "minimum " + w, used); put(b + "max", hi_[0], "maximum " + w, used)
                    put(b + "argmin", f"{d}|{lo_[1]}", "cell of the minimum " + w, used); put(b + "argmax", f"{d}|{hi_[1]}", "cell of the maximum " + w, used)
    for d, ci in FC:
        for r in RULES_ALL:
            for e in EPS:
                b = f"simbox/fc|{d}|{ci}|{r}|{e}|"
                if num(b + "n_cert") != num(b + "n"): raise ValueError(f"not every draw certified: {b}")
                put(f"fc_q95_gt_eps|{d}|{ci}|{r}|{e}", int(num(b + "q95") > float(e)), f"1 if the 95th percentile of realized FPR among certified draws exceeds eps, {d}|{ci}|{r}|{e}", [b + "q95"])
    return new


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); a = ap.parse_args()
    tree = Path(a.tree); rm = tree / "results_master.json"
    try:
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not effb684e")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        new = derive(vals)
        clash = [k for k in new if k in vals]
        if clash: raise ValueError(f"keys exist already: {clash[:5]}")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    bak = tree / "results_master_before_sim50_derived2_effb684e.json"; tmp = tree / "results_master.json.sim50_derived2_tmp"
    existing = [str(x) for x in (bak, tmp) if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}; old_fsha = dict(fsha)
    with open(rm, "rb") as s_, open(bak, "xb") as d_: d_.write(s_.read())
    for k, (v, s) in new.items(): vals[k] = v; srcs[k] = s
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp))
    try: again = derive({k: v for k, v in chk["values"].items() if k not in new})
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

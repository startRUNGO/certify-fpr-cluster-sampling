# P5a ADD-ONLY results_master builder: sim50 (50 synthetic sites per cell; sim50_extend.py ea3ea087, VM runs fifty_<UTC>, three batches).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Built from add_ddintv_keys.py (0a5fd9ed, reviewed): same write mechanics (write targets checked absent first; backup 'x'; tmp master re-read; every new key re-read from
#   its source path; every old key, source and files_sha256 entry unchanged; key count checked; indent=1).
# Input: the three batch directories pulled via the scratchpad (each: sim50_summary.json + SHA256SUMS, checked here).  SRC_PINS holds the summary sha of each batch;
#   an empty pin refuses to run (batches 2 and 3 are filled through a reviewed diff after their runs).  All three batches are required: the cell set must be complete.
# v1.1 (实验管理 2026-11-25): (1) size: per-realisation values are registered only for coverage and fpr_mean of all50 (F-b / F-c use only these); seeds lists are
#   not registered; the other stats keep mean and mcse (coverage also mcse_binomial and ci90).  The mean-vs-realisations gate still checks all six stats from the file.
#   (2) orig5 coverage mean must equal the published master key of the same cell, rule and eps (ORIG5_KEY, 1e-12); a cell without such a key is printed, not failed.
#   (3) Platt denominator: sim50/platt_fallback|n_tasks_total = number of task_seconds entries over the three batches (E10: 5 replayed + 45 new realisations per cell;
#   E12: the 45 new realisations; the E12 orig5 realisations reuse stored rows and are not re-fitted, so they are outside the Platt count), gated == 1580; every
#   platt_fallback entry must be a task of its batch.  (4) anchor b of a batch without E10 cells reads "not_applicable (no E10 cell)" in sim50_extend.py (underscore).
# Keys (prefix sim50/, results key as in the file, e.g. e10|0|S3|0.1):
#   all50 (the 50 realisations): sim50/<cell>|all50|<stat>|<field>[|i], filtered as in v1.1 (1) above (all50_leaf);
#   orig5 / new45: summary leaves only (n, n_realisations, <stat>|mean, <stat>|mcse, coverage|mcse_binomial, coverage|ci90_seedcluster|0/1); their per-realisation
#     values are the 5 / 45 entries of all50 and are not registered twice;
#   Platt fallback (derived, source = the list platt_fallback of each batch): sim50/platt_fallback|n_tasks (all batches), sim50/platt_fallback|<family>|<cell>|n_tasks
#     and the denominator sim50/platt_fallback|n_tasks_total (v1.1 (3)).
# Gates (exit 3, nothing written): master sha MASTER_SHA; SHA256SUMS of each batch dir; summary sha == SRC_PINS; version sim50_r1; script_sha256 == SCRIPT_SHA in every batch;
#   anchors a and c "passed", anchor b "passed" or "not_applicable" (only in a batch without E10 cells); exceptions empty; cells of the batches disjoint and their union
#   == EXPECTED_CELLS x RULES x EPS (34 cells x 5 rules x 2 eps = 340 result keys); per result: all50 seeds == orig5 seeds + new45 seeds, n_realisations 5 / 45 / 50,
#   |mean - mean(per_realisation)| <= 1e-12 for every stat of all50, ci90 lo <= mean <= hi for coverage; no existing key starts with sim50/; no non-finite value.
# Backup: results_master_before_sim50_<first 8 of MASTER_SHA>.json ('x').  Copies: <tree>/results_src/sim50/<batch dir name>/sim50_summary.json (+ files_sha256).
# Usage: python3 -E -B add_sim50_keys.py --tree <tree_bj_v1_...> --src <batch1 dir> <batch2 dir> <batch3 dir>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "09e168643816c024fccf77fc7306e3b7cb6cd88ea3cb118c3c4d16bd3f9ebeb5"
SCRIPT_SHA = "ea3ea087620e8bcda8d08a95a80db6897fd26622ede8c5842262699eb8b602a6"
# batch dir name -> sha of its sim50_summary.json; batches 2 and 3 are filled through a reviewed diff after their runs
SRC_PINS = {"fifty_20261124T214602Z": "7bb098ca3561a3216a39274bbce1dc88b8a0e96847fc5a1d898168bcbc2c29ab",
            "fifty_20261125T094940Z": "fa88e4f0102b1d90062a96164772c3b0cdcec43cf614dd982df8ac0e5b549d0f",
            "fifty_20261125T141231Z": "115968f1a0ce0b30c28cfa2eadc31d28c8c5cf5eeac9622c59ccfe647dd2f66e"}
EXPECTED_CELLS = [f"e10|{i}" for i in range(10)] + [f"e12|{i}" for i in range(24)]
RULES, EPS = ["S3", "SPLIT6", "CB_c1s", "GEE_c1s", "BB_c1s"], ["0.1", "0.05"]
STATS = ["coverage", "cannot_certify", "sens", "ds_eps", "ds_own", "fpr_mean"]
PRE, FNAME, DST = "sim50", "sim50_summary.json", "results_src/sim50"
TOL = 1e-12


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def leaves(o, path):
    if isinstance(o, dict):
        for k, v in o.items(): yield from leaves(v, path + [str(k)])
    elif isinstance(o, list):
        for i, v in enumerate(o): yield from leaves(v, path + [str(i)])
    else: yield path, o


def get(o, path):
    for p in path: o = o[int(p)] if isinstance(o, list) else o[p]
    return o


# v1.1: published master key of the orig5 coverage mean of a result key (design|cell|rule|eps)
def ORIG5_KEY(rk):
    d, ci, rule, eps = rk.split("|")
    if rule in ("S3", "SPLIT6"):
        return f"ebatch_e10/cells|{ci}|library_rules|{rule}|{eps}|coverage" if d == "e10" else f"e12/cells|{ci}|results|{rule}|{eps}|coverage"
    return f"simw3/{d}|{ci}|{rule}|{eps}|coverage|mean"


# v1.1: all50 keeps the summary leaves plus the per-realisation values of coverage and fpr_mean
def all50_leaf(rel): return summary_leaf(rel) or (len(rel) == 3 and rel[0] in ("coverage", "fpr_mean") and rel[1] == "per_realisation")


N_TASKS_TOTAL = 10 * 50 + 24 * 45


def summary_leaf(rel):
    """orig5 / new45: keep n, n_realisations, <stat>|mean, <stat>|mcse, coverage|mcse_binomial, coverage|ci90_seedcluster|i."""
    if rel in (["n"], ["n_realisations"]): return True
    if len(rel) == 2 and rel[0] in STATS and rel[1] in ("mean", "mcse"): return True
    return rel[:1] == ["coverage"] and (rel[1:] == ["mcse_binomial"] or rel[1:2] == ["ci90_seedcluster"])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); ap.add_argument("--src", required=True, nargs="+"); a = ap.parse_args()
    tree = Path(a.tree); rm = tree / "results_master.json"; dirs = [Path(s) for s in a.src]
    try:
        if len(dirs) != 3 or sorted(d.name for d in dirs) != sorted(SRC_PINS): raise ValueError(f"--src must be the three batch dirs {sorted(SRC_PINS)}")
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master differs from MASTER_SHA")
        objs = {}
        for d in dirs:
            if not SRC_PINS[d.name]: raise ValueError(f"SRC_PINS[{d.name}] is empty: fill it through a reviewed diff after the run")
            rows = [l.split("  ", 1) for l in (d / "SHA256SUMS").read_text().split("\n") if l.strip()]
            if any(len(r) != 2 or sha(d / r[1]) != r[0] for r in rows) or FNAME not in {r[1] for r in rows}: raise ValueError(f"SHA256SUMS check failed in {d}")
            if sha(d / FNAME) != SRC_PINS[d.name]: raise ValueError(f"summary sha differs from the pin in {d}")
            o = json.load(open(d / FNAME))
            has_e10 = any(k.startswith("e10|") for k in o["results"])
            if o["version"] != "sim50_r1" or o["script_sha256"] != SCRIPT_SHA or o["exceptions"]: raise ValueError(f"version / script / exceptions in {d}")
            an = o["anchors"]
            if not (an["a"].startswith("passed") and an["c"].startswith("passed") and (an["b"].startswith("passed") or (not has_e10 and an["b"].startswith("not_applicable")))):
                raise ValueError(f"anchors in {d}: {an}")
            objs[d.name] = o
        seen, want = {}, {f"{c}|{r}|{e}" for c in EXPECTED_CELLS for r in RULES for e in EPS}
        for bn, o in objs.items():
            for rk in o["results"]:
                if rk in seen: raise ValueError(f"result {rk} in two batches")
                seen[rk] = bn
        if set(seen) != want: raise ValueError(f"result keys differ from the expected 340: missing {sorted(want - set(seen))[:5]}, extra {sorted(set(seen) - want)[:5]}")
        for rk, bn in seen.items():
            g = objs[bn]["results"][rk]
            if [g["orig5"]["n_realisations"], g["new45"]["n_realisations"], g["all50"]["n_realisations"]] != [5, 45, 50]: raise ValueError(f"realisation counts {rk}")
            if g["all50"]["seeds"] != g["orig5"]["seeds"] + g["new45"]["seeds"]: raise ValueError(f"all50 seeds differ from orig5 + new45 {rk}")
            for st in STATS:
                pr = g["all50"][st]["per_realisation"]
                if len(pr) != 50 or abs(g["all50"][st]["mean"] - sum(pr) / 50) > TOL: raise ValueError(f"all50 mean differs from its realisations {rk} {st}")
            lo, hi = g["all50"]["coverage"]["ci90_seedcluster"]
            if not lo <= g["all50"]["coverage"]["mean"] <= hi: raise ValueError(f"coverage outside its interval {rk}")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        if any(k.startswith(PRE + "/") for k in vals): raise ValueError(f"existing keys already use prefix {PRE}/")
        # v1.1: orig5 coverage mean == published master key (missing keys are reported, not failed)
        no_key = []
        for rk, bn in seen.items():
            k = ORIG5_KEY(rk)
            if k not in vals: no_key.append(rk); continue
            if abs(vals[k] - objs[bn]["results"][rk]["orig5"]["coverage"]["mean"]) > TOL: raise ValueError(f"orig5 coverage of {rk} differs from {k}")
        print(f"orig5 coverage vs master: {len(seen) - len(no_key)} equal, {len(no_key)} without a master key", no_key, file=sys.stderr)
        # v1.1: Platt denominator and scope
        n_tasks = sum(len(o["task_seconds"]) for o in objs.values())
        if n_tasks != N_TASKS_TOTAL: raise ValueError(f"task count {n_tasks}, expected {N_TASKS_TOTAL}")
        for bn, o in objs.items():
            tk = {tuple(k.split("|")[1:]) for k in o["task_seconds"]}
            if any((f, str(c), str(s)) not in tk for f, c, s in o["platt_fallback"]): raise ValueError(f"a platt_fallback entry of {bn} is not a task of that batch")
    except (ValueError, KeyError, OSError, TypeError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    new, derived = {}, {}
    for rk, bn in seen.items():
        rel_file = f"{DST}/{bn}/{FNAME}"
        for grp in ("orig5", "new45", "all50"):
            for path, v in leaves(objs[bn]["results"][rk][grp], ["results", rk, grp]):
                rel = path[3:]
                if not (all50_leaf(rel) if grp == "all50" else summary_leaf(rel)): continue
                if isinstance(v, bool) or not isinstance(v, (int, float)): continue
                if isinstance(v, float) and not math.isfinite(v): print("GATE (exit 3): non-finite", path, file=sys.stderr); return 3
                key = PRE + "/" + "|".join([rk, grp] + rel)
                if key in vals or key in new: print("GATE (exit 3): duplicate key", key, file=sys.stderr); return 3
                new[key] = (v, rel_file, path)
    per_cell = {}
    for bn, o in objs.items():
        for fam, cell, *_ in o["platt_fallback"]: per_cell[f"{fam}|{cell}"] = per_cell.get(f"{fam}|{cell}", 0) + 1
    derived[f"{PRE}/platt_fallback|n_tasks"] = (sum(per_cell.values()), "number of entries of platt_fallback summed over the three batches")
    derived[f"{PRE}/platt_fallback|n_tasks_total"] = (n_tasks, "number of entries of task_seconds summed over the three batches (E10: 5 replayed + 45 new realisations per cell; "
                                                     "E12: the 45 new realisations; E12 orig5 realisations reuse stored rows and are not re-fitted)")
    for c, n in sorted(per_cell.items()): derived[f"{PRE}/platt_fallback|{c}|n_tasks"] = (n, f"number of entries of platt_fallback with family|cell {c}, summed over the three batches")
    tag = MASTER_SHA[:8]
    bak = tree / f"results_master_before_sim50_{tag}.json"; tmp = tree / "results_master.json.sim50_tmp"
    dsts = [tree / DST / bn / FNAME for bn in objs]
    existing = [str(x) for x in [bak, tmp] + dsts if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}; old_fsha = dict(fsha)
    with open(rm, "rb") as s_, open(bak, "xb") as d_: d_.write(s_.read())
    for d_dir in dirs:
        (tree / DST / d_dir.name).mkdir(parents=True, exist_ok=True)
        with open(d_dir / FNAME, "rb") as s_, open(tree / DST / d_dir.name / FNAME, "xb") as d_: d_.write(s_.read())
        if sha(tree / DST / d_dir.name / FNAME) != SRC_PINS[d_dir.name]: print("GATE (exit 3): copy sha differs", d_dir.name, file=sys.stderr); return 3
        fsha[f"{DST}/{d_dir.name}/{FNAME}"] = SRC_PINS[d_dir.name]
    for k, (v, rel_file, path) in new.items(): vals[k] = v; srcs[k] = dict(file=rel_file, path=path)
    for k, (v, how) in derived.items():
        if k in vals: print("GATE (exit 3): duplicate derived key", k, file=sys.stderr); return 3
        vals[k] = v; srcs[k] = dict(files=sorted(f"{DST}/{bn}/{FNAME}" for bn in objs), derived=how)
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp))
    src_objs = {bn: json.load(open(tree / DST / bn / FNAME)) for bn in objs}
    for k, (v, rel_file, path) in new.items():
        bn = rel_file.split("/")[2]
        if get(src_objs[bn], path) != chk["values"][k]: print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    for k, v in old_fsha.items():
        if chk["files_sha256"].get(k) != v: print("GATE (exit 3): an old files_sha256 entry changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new) + len(derived): print("GATE (exit 3): key count mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    print(json.dumps(dict(added=len(new) + len(derived), from_results=len(new), derived=len(derived), n_values=d["n_values"], new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

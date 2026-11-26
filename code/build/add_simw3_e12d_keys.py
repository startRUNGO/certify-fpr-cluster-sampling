# P5a ADD-ONLY results_master builder: simulation baselines simw3/* (simw3_baselines.py 45f1ebcd) and E12 derived sensitivity e12d/* (e12_derive_sens.py 6b1f3c65).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Same mechanics as add_figmech_v39_keys.py (52ea5b64): write targets checked absent first; indent=1; every new key re-read from its source path; every old key, source and files_sha256 entry unchanged; key count checked.
# Inputs:
#   simw3_summary.json  (full run <DATA>/p5_simw3/simw3_20261124T194046Z, pulled via the scratchpad, SHA256SUMS checked)  -> keys simw3/<p1>|...|<pn> for every int / float / bool leaf under "results" and "settings"
#       (strings, None, task_seconds, seconds, cells_spec are not registered).  Its pin SIMW3_SHA is EMPTY until the full run is done; it is filled only through a reviewed diff, and the script refuses to run while it is empty.
#       Gate: anchors b, c, d and e of the summary all start with "passed" (e must pass at full scale).
#   e12_derived_sens.json  ea9a05d4 (scratchpad e12derive_20261124T194128Z, anchor passed)  -> keys e12d/<p1>|...|<pn> for every int / float / bool leaf under "cells".
#   simw3_rows.json is NOT copied into the tree (about 80 MB); its sha is recorded in files_sha256 under its run path for traceability only if SIMW3_ROWS_SHA is filled (same rule as SIMW3_SHA).
# Both summaries are copied to <tree>/results_src/simw3/ and <tree>/results_src/e12derive/ and added to files_sha256.  No existing key is touched (add-only).
# Gates (exit 3, nothing written): master sha 1f68e729; input shas; no existing key starts with simw3/ or e12d/; every new key new; no non-finite value; all write targets absent.
# Backup: results_master_before_simw3_e12d_1f68e729.json ('x').
# Usage: python3 -B add_simw3_e12d_keys.py --tree <tree_bj_v1_...> --simw3 <dir with simw3_summary.json> --e12d <dir with e12_derived_sens.json>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "1f68e729ba8c2f499021363ea2bc59f6d6776bc64483559d0a0a5f37dc6e13bc"
SIMW3_SHA = "495a9b29f6b871af74237d39e956dab73fd3a50c7b0919ac66e254768b850b13"
SIMW3_ROWS_SHA = "bed0438a1f0717d98482502a0baa93d2fddfefb88e3179048e44ea8bff2d3794"
E12D_SHA = "ea9a05d487e47de2797e007d758ae695ff9e3182e865e9a710bb52b2e12bcd86"
INPUTS = (("simw3", "simw3_summary.json", "results_src/simw3", ("results", "settings")), ("e12d", "e12_derived_sens.json", "results_src/e12derive", ("cells",)))


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


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); ap.add_argument("--simw3", required=True); ap.add_argument("--e12d", required=True); a = ap.parse_args()
    tree = Path(a.tree); rm = tree / "results_master.json"; src_dir = {"simw3": Path(a.simw3), "e12d": Path(a.e12d)}; pin = {"simw3": SIMW3_SHA, "e12d": E12D_SHA}
    try:
        if not SIMW3_SHA: raise ValueError("SIMW3_SHA is empty: fill it through a reviewed diff after the full run")
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not 1f68e729")
        for pre, f, _, _ in INPUTS:
            if sha(src_dir[pre] / f) != pin[pre]: raise ValueError(f"input sha differs: {f}")
        summ = json.load(open(src_dir["simw3"] / INPUTS[0][1]))
        bad = [k for k in ("b", "c", "d", "e") if not str(summ["anchors"].get(k, "")).startswith("passed")]
        if bad: raise ValueError(f"simw3 anchors not passed: {bad}")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        for pre, *_ in INPUTS:
            if any(k.startswith(pre + "/") for k in vals): raise ValueError(f"existing keys already use prefix {pre}/")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    new = {}
    for pre, f, dst, tops in INPUTS:
        obj = json.load(open(src_dir[pre] / f)); rel = f"{dst}/{f}"
        for top in tops:
            for path, v in leaves(obj[top], [top]):
                if isinstance(v, bool) or isinstance(v, (int, float)):
                    if isinstance(v, float) and not math.isfinite(v): print("GATE (exit 3): non-finite", pre, path, file=sys.stderr); return 3
                    key = pre + "/" + "|".join(path[1:])
                    if key in vals or key in new: print("GATE (exit 3): duplicate key", key, file=sys.stderr); return 3
                    new[key] = (v, rel, path)
    bak = tree / "results_master_before_simw3_e12d_1f68e729.json"; tmp = tree / "results_master.json.simw3_e12d_tmp"
    targets = [bak, tmp] + [tree / dst / f for _, f, dst, _ in INPUTS]
    existing = [str(x) for x in targets if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}; old_fsha = dict(fsha)
    with open(rm, "rb") as src_fh, open(bak, "xb") as dst_fh: dst_fh.write(src_fh.read())
    for pre, f, dst, _ in INPUTS:
        (tree / dst).mkdir(parents=True, exist_ok=True)
        with open(src_dir[pre] / f, "rb") as src_fh, open(tree / dst / f, "xb") as dst_fh: dst_fh.write(src_fh.read())
        if sha(tree / dst / f) != pin[pre]: print("GATE (exit 3): copy sha differs", f, file=sys.stderr); return 3
        fsha[f"{dst}/{f}"] = pin[pre]
    if SIMW3_ROWS_SHA: fsha["runs/simw3_20261124T194046Z/simw3_rows.json"] = SIMW3_ROWS_SHA
    for k, (v, rel_, path) in new.items(): vals[k] = v; srcs[k] = dict(file=rel_, path=path)
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp)); src_obj = {f"{dst}/{f}": json.load(open(tree / dst / f)) for _, f, dst, _ in INPUTS}
    for k, (v, rel_, path) in new.items():
        if get(src_obj[rel_], path) != chk["values"][k]: print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    for k, v in old_fsha.items():
        if chk["files_sha256"].get(k) != v: print("GATE (exit 3): an old files_sha256 entry changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new): print("GATE (exit 3): key count mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    print(json.dumps(dict(added=len(new), added_simw3=sum(k.startswith("simw3/") for k in new), added_e12d=sum(k.startswith("e12d/") for k in new), n_values=d["n_values"], new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

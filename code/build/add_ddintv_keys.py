# P5a ADD-ONLY results_master builder: D-d development-site intervals d_d/* (d_d_dev_intervals.py fa9902c5, VM run ddintv_20261124T232236Z).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Built from add_dafunnel_keys.py (f54f8bcd, reviewed); the constants, the key naming (readable, below) and the consistency gate differ.
# Same mechanics as add_simw3_e12d_keys.py (bd7fb28f): write targets checked absent first; indent=1; every new key re-read from its source path;
#   every old key, source and files_sha256 entry unchanged; key count checked.
# Input: d_d_dev_intervals.json (pulled via the scratchpad, SHA256SUMS checked; SRC_SHA is EMPTY until the run is done and is filled through a reviewed diff; the builder refuses
#   to run while it is empty) -> for every cell c of 'cells' and every int / float / bool leaf f of c: key d_d/<rule>|<site>|<eps>|<f>[|<i>] (eps as in the file, e.g. 0.1),
#   source path ['cells', index, f(, i)].  The file is copied to <tree>/results_src/ddintv/ and added to files_sha256.
# Gates (exit 3, nothing written): master sha fccda464; input sha; no existing key starts with d_a/; every new key new; no non-finite value; all write targets absent;
#   d_d values: 72 cells (12 rules x 3 sites x 2 eps), each n == 1000, 5 per-seed values, lo <= point <= hi for coverage and Delta Sens_eps.
# Backup: results_master_before_ddintv_e459b90f.json ('x').
# Usage: python3 -E -B add_ddintv_keys.py --tree <tree_bj_v1_...> --src <dir with d_d_dev_intervals.json>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "e459b90f5dd3fb398e9acef2e8f7d2815a498d1475f82d1ca896b4592971830a"
# filled through a reviewed diff after the run
SRC_SHA = "8653b50bca20949b8b6e86d796ab52d2b554577b21e8b1b34045387b363551df"
PRE, FNAME, DST = "d_d", "d_d_dev_intervals.json", "results_src/ddintv"


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
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); ap.add_argument("--src", required=True); a = ap.parse_args()
    tree = Path(a.tree); rm = tree / "results_master.json"; src = Path(a.src) / FNAME
    try:
        if not SRC_SHA: raise ValueError("SRC_SHA is empty: fill it through a reviewed diff after the run")
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not e459b90f")
        if sha(src) != SRC_SHA: raise ValueError("input sha differs")
        obj = json.load(open(src))
        c = obj["cells"]
        ok = lambda x: x["n"] == 1000 and len(x["per_seed_coverage"]) == 5 and len(x["per_seed_ds_eps"]) == 5 and x["coverage_ci90_seedcluster"][0] <= x["coverage"] <= x["coverage_ci90_seedcluster"][1] \
            and x["ds_eps_ci90_seedcluster"][0] <= x["ds_eps"] <= x["ds_eps_ci90_seedcluster"][1]
        if len(c) != 72 or len({(x["rule"], x["site"], x["eps"]) for x in c}) != 72 or not all(ok(x) for x in c): raise ValueError("d_d cells inconsistent")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        if any(k.startswith(PRE + "/") for k in vals): raise ValueError(f"existing keys already use prefix {PRE}/")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    rel = f"{DST}/{FNAME}"; new = {}
    for i, cell in enumerate(obj["cells"]):
        for path, v in leaves(cell, ["cells", str(i)]):
            if isinstance(v, bool) or isinstance(v, (int, float)):
                if isinstance(v, float) and not math.isfinite(v): print("GATE (exit 3): non-finite", path, file=sys.stderr); return 3
                key = PRE + "/" + "|".join([cell["rule"], cell["site"], str(cell["eps"])] + path[2:])
                if key in vals or key in new: print("GATE (exit 3): duplicate key", key, file=sys.stderr); return 3
                new[key] = (v, path)
    bak = tree / "results_master_before_ddintv_e459b90f.json"; tmp = tree / "results_master.json.ddintv_tmp"; dst = tree / DST / FNAME
    existing = [str(x) for x in (bak, tmp, dst) if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}; old_fsha = dict(fsha)
    with open(rm, "rb") as s_, open(bak, "xb") as d_: d_.write(s_.read())
    (tree / DST).mkdir(parents=True, exist_ok=True)
    with open(src, "rb") as s_, open(dst, "xb") as d_: d_.write(s_.read())
    if sha(dst) != SRC_SHA: print("GATE (exit 3): copy sha differs", file=sys.stderr); return 3
    fsha[rel] = SRC_SHA
    for k, (v, path) in new.items(): vals[k] = v; srcs[k] = dict(file=rel, path=path)
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp)); src_obj = json.load(open(dst))
    for k, (v, path) in new.items():
        if get(src_obj, path) != chk["values"][k]: print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    for k, v in old_fsha.items():
        if chk["files_sha256"].get(k) != v: print("GATE (exit 3): an old files_sha256 entry changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new): print("GATE (exit 3): key count mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    print(json.dumps(dict(added=len(new), n_values=d["n_values"], new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

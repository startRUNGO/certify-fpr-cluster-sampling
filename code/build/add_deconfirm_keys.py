# P5a ADD-ONLY results_master builder: D-e confirmation-site intervals d_e/* (d_e_confirm_intervals.py a9f40100, local run deconfirm_<UTC>).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Built from add_dafunnel_keys.py (f54f8bcd, reviewed); only the constants, the input description and the consistency gate differ.
# Same mechanics as add_simw3_e12d_keys.py (bd7fb28f): write targets checked absent first; indent=1; every new key re-read from its source path;
#   every old key, source and files_sha256 entry unchanged; key count checked.
# Input: d_e_confirm_intervals.json (local run in the scratchpad, SHA256SUMS checked) -> keys d_e/cells|<i>|<field>... for every int / float / bool leaf of 'cells'
#   (point, per_seed|0..4, ci90_seedcluster|0/1, eps; strings rule / site / stat / master_key are not registered; the readable index is in the file).
#   The file is copied to <tree>/results_src/deconfirm/ and added to files_sha256.
# Gates (exit 3, nothing written): master sha fccda464; input sha; no existing key starts with d_a/; every new key new; no non-finite value; all write targets absent;
#   d_e values: 8 cells (IML; S3, R1_2; eps 0.1, 0.05; cov, dsens), each with 5 per-seed values and lo <= point <= hi.
# Backup: results_master_before_deconfirm_21496cc6.json ('x').
# Usage: python3 -E -B add_deconfirm_keys.py --tree <tree_bj_v1_...> --src <dir with d_e_confirm_intervals.json>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "21496cc647faa98ef17ee0261b7da65c766f5dbfe34aae64290396bf2e35124c"
SRC_SHA = "21dd87da35a552a950dedb34ff866ac5e87ea3bb811e6c04b1167e69c1c293a5"
PRE, FNAME, DST = "d_e", "d_e_confirm_intervals.json", "results_src/deconfirm"
SKIP_TOP = {"pins", "nature", "script_sha256", "no_interval", "stored_interval"}


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
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not 21496cc6")
        if sha(src) != SRC_SHA: raise ValueError("input sha differs")
        obj = json.load(open(src))
        c = obj["cells"]
        if len(c) != 8 or any(len(x["per_seed"]) != 5 or not (x["ci90_seedcluster"][0] <= x["point"] <= x["ci90_seedcluster"][1]) for x in c): raise ValueError("d_e cells inconsistent")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        if any(k.startswith(PRE + "/") for k in vals): raise ValueError(f"existing keys already use prefix {PRE}/")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    rel = f"{DST}/{FNAME}"; new = {}
    for top, v0 in obj.items():
        if top in SKIP_TOP: continue
        for path, v in leaves(v0, [top]):
            if isinstance(v, bool) or isinstance(v, (int, float)):
                if isinstance(v, float) and not math.isfinite(v): print("GATE (exit 3): non-finite", path, file=sys.stderr); return 3
                key = PRE + "/" + "|".join(path)
                if key in vals or key in new: print("GATE (exit 3): duplicate key", key, file=sys.stderr); return 3
                new[key] = (v, path)
    bak = tree / "results_master_before_deconfirm_21496cc6.json"; tmp = tree / "results_master.json.deconfirm_tmp"; dst = tree / DST / FNAME
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

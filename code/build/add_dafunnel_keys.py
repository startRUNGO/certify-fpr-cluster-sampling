# P5a ADD-ONLY results_master builder: D-a funnel export d_a/* (d_a_funnel_bbbc.py 0ff2f656, VM run dafunnel_20261124T230844Z).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Same mechanics as add_simw3_e12d_keys.py (bd7fb28f): write targets checked absent first; indent=1; every new key re-read from its source path;
#   every old key, source and files_sha256 entry unchanged; key count checked.
# Input: d_a_funnel_bbbc.json 79ffe00f (pulled via the scratchpad, SHA256SUMS checked) -> keys d_a/<p1>|...|<pn> for every int / float / bool leaf of the file
#   (top-level scalars and the lists n_f, X_f, frac_f, funnel_n, funnel_q05_frac, funnel_q95_frac; strings and the nested dicts input_sha256 / anchor strings are not registered;
#   anchor numbers are registered).  The file is copied to <tree>/results_src/dafunnel/ and added to files_sha256.
# Gates (exit 3, nothing written): master sha fccda464; input sha; no existing key starts with d_a/; every new key new; no non-finite value; all write targets absent;
#   d_a values: n_fields == len(n_f) == len(X_f); sum(n_f) == n_neg.
# Backup: results_master_before_dafunnel_fccda464.json ('x').
# Usage: python3 -E -B add_dafunnel_keys.py --tree <tree_bj_v1_...> --src <dir with d_a_funnel_bbbc.json>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "fccda46437a84c34783f0e355d88d2e94487e8effa91cdc85eaf381a28517642"
SRC_SHA = "79ffe00f62502ecd61e12f226ec7c25e1685e13e5dd3972c90691aee39ac6648"
PRE, FNAME, DST = "d_a", "d_a_funnel_bbbc.json", "results_src/dafunnel"
SKIP_TOP = {"input_sha256"}


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
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not fccda464")
        if sha(src) != SRC_SHA: raise ValueError("input sha differs")
        obj = json.load(open(src))
        if not (obj["n_fields"] == len(obj["n_f"]) == len(obj["X_f"])) or sum(obj["n_f"]) != obj["n_neg"]: raise ValueError("d_a field counts inconsistent")
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
    bak = tree / "results_master_before_dafunnel_fccda464.json"; tmp = tree / "results_master.json.dafunnel_tmp"; dst = tree / DST / FNAME
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

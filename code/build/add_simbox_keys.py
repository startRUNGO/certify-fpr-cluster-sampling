# P5a ADD-ONLY results_master builder: D-c simulation box statistics and uncertified-draw counts simbox/* (d_c_simbox.py 17223230, local run d_c_simbox_20261125T191523Z).
#   WRITTEN, NOT RUN, until 实验管理 reviews it.
# Built from add_ddintv_keys.py (reviewed); the constants, the key naming (below) and the consistency gate differ.  Same mechanics: write targets checked absent first;
#   indent=1; every new key re-read from its source path; every old key, source and files_sha256 entry unchanged; key count checked.
# Input: d_c_simbox.json (SHA256SUMS checked after the run; SRC_SHA filled from that run, shown in the review diff) -> keys
#   simbox/fc|<d>|<ci>|<rule>|<eps>|<f>[|<i>]   for every int / float leaf f of fc[<d>|<ci>|<rule>|<eps>] (n, n_cert, q05 .. q95, mean, min, max, strip_values|i if present),
#     source path ['fc', '<d>|<ci>|<rule>|<eps>', f(, i)];
#   simbox/uncert|<d>|<ci>|<rule>|<eps>|<f>    for every leaf f of counts[...] (n, n_cert, n_uncert, n_uncert_fpr_gt_eps), source path ['counts', key, f];
#   simbox/uncert_total|<eps>                    = n_uncert_fpr_gt_eps_total[eps], source path ['n_uncert_fpr_gt_eps_total', eps].
#   The file is copied to <tree>/results_src/simbox/ and added to files_sha256.
# Gates (exit 3, nothing written): master sha 9b04913e; input sha; script_sha256 in the file == DC_SHA; no existing key starts with simbox/; every new key new;
#   no non-finite value; all write targets absent; consistency: 340 counts with n == 10000 and n_cert + n_uncert == n and 0 <= n_uncert_fpr_gt_eps <= n_uncert;
#   the totals equal the sums over the counts; 40 fc entries (4 settings x 5 rules x 2 eps), each with n == 10000, n_cert equal to its counts entry,
#   min <= q05 <= q25 <= q50 <= q75 <= q95 <= max when n_cert > 0, strip_values present iff n_cert < 100 and then of length n_cert.
# Backup: results_master_before_simbox_9b04913e.json ('x').
# Usage: python3 -E -B add_simbox_keys.py --tree <tree_bj_v1_...> --src <dir with d_c_simbox.json>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "9b04913eba0b2044e12fce07b935de17a7e5065fb616ee636dad9c54fb6fec4a"
# from the run d_c_simbox_20261125T191523Z (SHA256SUMS checked)
SRC_SHA = "0d6b3b55b690ca69aba1712f47b9936e317be5001a91ac0ebfcb6f6d59b27ab9"
DC_SHA = "17223230b7fc8a1721e8cd52002ce89b4af3a78084e21b1622121e3858a81970"
PRE, FNAME, DST = "simbox", "d_c_simbox.json", "results_src/simbox"
Q = ("min", "q05", "q25", "q50", "q75", "q95", "max")


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
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not 9b04913e")
        if sha(src) != SRC_SHA: raise ValueError("input sha differs")
        obj = json.load(open(src))
        if obj["script_sha256"] != DC_SHA: raise ValueError("d_c_simbox.json was not written by the reviewed d_c_simbox.py")
        C, F, T = obj["counts"], obj["fc"], obj["n_uncert_fpr_gt_eps_total"]
        okc = lambda x: x["n"] == 10000 and x["n_cert"] + x["n_uncert"] == x["n"] and 0 <= x["n_uncert_fpr_gt_eps"] <= x["n_uncert"]
        if len(C) != 340 or not all(okc(x) for x in C.values()): raise ValueError("counts inconsistent")
        if T != {e: sum(x["n_uncert_fpr_gt_eps"] for k, x in C.items() if k.endswith("|" + e)) for e in ("0.1", "0.05")}: raise ValueError("totals differ from the counts")
        def okf(k, x):
            if x["n"] != 10000 or x["n_cert"] != C[k]["n_cert"]: return False
            if x["n_cert"] > 0 and not all(x[a] <= x[b] for a, b in zip(Q, Q[1:])): return False
            return ("strip_values" in x) == (x["n_cert"] < 100) and ("strip_values" not in x or len(x["strip_values"]) == x["n_cert"])
        if len(F) != 40 or not all(k in C and okf(k, x) for k, x in F.items()): raise ValueError("fc entries inconsistent")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        if any(k.startswith(PRE + "/") for k in vals): raise ValueError(f"existing keys already use prefix {PRE}/")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    rel = f"{DST}/{FNAME}"; new = {}
    items = [(f"{PRE}/fc|{k}", ["fc", k], x) for k, x in F.items()] + [(f"{PRE}/uncert|{k}", ["counts", k], x) for k, x in C.items()]
    items += [(f"{PRE}/uncert_total|{e}", ["n_uncert_fpr_gt_eps_total", e], T[e]) for e in ("0.1", "0.05")]
    for base, p0, o in items:
        for path, v in (leaves(o, p0) if isinstance(o, dict) else [(p0, o)]):
            if isinstance(v, bool) or not isinstance(v, (int, float)): print("GATE (exit 3): non-numeric leaf", path, file=sys.stderr); return 3
            if isinstance(v, float) and not math.isfinite(v): print("GATE (exit 3): non-finite", path, file=sys.stderr); return 3
            key = "|".join([base] + path[len(p0):])
            if key in vals or key in new: print("GATE (exit 3): duplicate key", key, file=sys.stderr); return 3
            new[key] = (v, path)
    bak = tree / "results_master_before_simbox_9b04913e.json"; tmp = tree / "results_master.json.simbox_tmp"; dst = tree / DST / FNAME
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

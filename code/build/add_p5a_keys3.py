# P5a (CertBound) ADD-ONLY results_master builder, round 3: the preregistered M5-1 thresholds of M2 and M3 (for sections 5-7, [[TBD:M2 threshold]] / [[TBD:M3 threshold]]).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Same mechanics as add_p5a_keys2.py v1.2 (7e78c5c8): write targets checked absent first; indent=1; every new key re-read from its source path; every old key/source unchanged; key count checked.
# Input (P5a_CertBound/inputs_20261124/):  prefix m51gate/   file m51gate/m51_gates.json   sha GATE_SHA (filled after make_m51_gates.py has run and been reported; empty -> exit 3)   all leaves.
# Keys: m51gate/M2|hw_ratio_max (0.80), m51gate/M2|min_sites (2), m51gate/M3|sens_diff_min (0.01), m51gate/M3|min_sites (2), m51gate/n_dev_sites (3).
# Gates (exit 3, nothing written): master sha 0f5217e6 (round-2 result); input sha; no existing key starts with m51gate/; every new key new.
# Backup: results_master_before_p5a3_0f5217e6.json ('x').  Input copied to <tree>/results_src/p5a3/<file> and added to files_sha256.
# Usage: python3 -B add_p5a_keys3.py --tree <tree_bj_v1_...> --inputs <inputs_20261124>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "0f5217e66aec2df1af44afd4fdec5dd424f49ac73e23ef293103eeb04ec8ccb7"
GATE_SHA = "9d32fd25d41d40f3094576845e1cbda60d7724790c05e2c1947ba43e78179213"   # v1.1: output of make_m51_gates.py cdefe9b0 (rc 0)
SPEC = [("m51gate/", "m51gate/m51_gates.json", GATE_SHA, None)]   # None = all leaves


def sha(p):
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
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); ap.add_argument("--inputs", required=True); a = ap.parse_args(); tree, inp = Path(a.tree), Path(a.inputs)
    rm = tree / "results_master.json"
    try:
        if GATE_SHA == "": raise ValueError("GATE_SHA is empty (make_m51_gates.py not yet run / reported)")
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not the round-2 master 0f5217e6")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        for pre, f, h, flds in SPEC:
            if sha(inp / f) != h: raise ValueError(f"input sha differs: {f}")
            if any(k.startswith(pre) for k in vals): raise ValueError(f"existing keys already use prefix {pre}")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    new = {}; per_prefix = {}
    for pre, f, h, flds in SPEC:
        obj = json.load(open(inp / f)); rel = f"results_src/p5a3/{f}"; n = 0
        if flds is not None:
            miss = [x for x in flds if x not in obj]
            if miss: print("GATE (exit 3): listed field missing", f, miss, file=sys.stderr); return 3
            items = [(p, v) for x in flds for p, v in leaves(obj[x], [x])]
        else: items = list(leaves(obj, []))
        for path, v in items:
            if isinstance(v, bool) or isinstance(v, (int, float)):
                if isinstance(v, float) and not math.isfinite(v): print("GATE (exit 3): non-finite", f, path, file=sys.stderr); return 3
                key = pre + "|".join(path)
                if key in vals or key in new: print("GATE (exit 3): duplicate key", key, file=sys.stderr); return 3
                new[key] = (v, rel, path); n += 1
        per_prefix[pre] = n
    bak = tree / "results_master_before_p5a3_0f5217e6.json"; tmp = tree / "results_master.json.p5a3_tmp"
    targets = [bak, tmp] + [tree / "results_src" / "p5a3" / f for pre, f, h, flds in SPEC]
    existing = [str(x) for x in targets if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}
    with open(rm, "rb") as src_fh, open(bak, "xb") as dst_fh: dst_fh.write(src_fh.read())
    for pre, f, h, flds in SPEC:
        dst = tree / "results_src" / "p5a3" / f; dst.parent.mkdir(parents=True, exist_ok=True)
        with open(inp / f, "rb") as src_fh, open(dst, "xb") as dst_fh: dst_fh.write(src_fh.read())
        fsha[f"results_src/p5a3/{f}"] = h
    for k, (v, rel, path) in new.items(): vals[k] = v; srcs[k] = dict(file=rel, path=path)
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp)); cache = {}
    for k, (v, rel, path) in new.items():
        if rel not in cache: cache[rel] = json.load(open(tree / rel))
        if get(cache[rel], path) != chk["values"][k]: print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new): print("GATE (exit 3): key count mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    print(json.dumps(dict(added=len(new), per_prefix=per_prefix, n_values=d["n_values"], new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

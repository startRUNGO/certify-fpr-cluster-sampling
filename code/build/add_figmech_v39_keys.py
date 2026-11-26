# P5a v3.9 ADD-ONLY results_master builder: Figure 2 (fig:mechanism) keys figmech_v39/*.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Same mechanics as add_p5a_keys4.py v1.3: write targets checked absent first; indent=1; every new key re-read from its source path; every old key/source unchanged; key count checked.
# Inputs (--inputs = results_src/c1card/figmech_v39_out_20261124T131527Z, pulled via the scratchpad, SHA256SUMS 69cf0a48 verified twice; export v1.1 89e909a5, rc 0):
#   figmech_v39_numbers.json  a4b23923...  -> keys figmech_v39/<p1>|...|<pn>, every int / float / bool leaf (strings such as field ids, candidate names and shas, and None, are not registered)
#   figmech_v39_arrays.npz    3ee12fc9...  -> no keys; copied into the tree and listed in files_sha256 so that the figure script can pin it
# Both are copied to <tree>/results_src/figmech_v39/<file> and added to files_sha256.
# The 42 old fig_mechanism_numbers/* keys are NOT touched (add-only).  The old figure-side copy <tree>/figures/fig_mechanism_numbers.json (0b7eb8b2, the pre-THRFIX panel (a) numbers)
#   is renamed to <tree>/figures/NOT_USED_pre_thrfix_fig_mechanism_numbers.json (no rm) as the LAST step, after the new master is in place; the master sources of the old keys point to
#   results_src/f0/fig_mechanism_numbers.json, not to this copy, so no master source changes.
# Gates (exit 3, nothing written): master sha b5faa828; the two input shas; old json sha 0b7eb8b2; no existing key starts with figmech_v39/; every new key new; no non-finite value;
#   all write targets (backup, tmp, two copies, renamed old json) absent.
# Backup: results_master_before_figmech_v39_b5faa828.json ('x').
# Usage: python3 -B add_figmech_v39_keys.py --tree <tree_bj_v1_...> --inputs <figmech_v39_out_...>
import argparse, hashlib, json, math, os, sys
from pathlib import Path
MASTER_SHA = "b5faa82817383bc977e03617e6e85e53ad15b9347a079d06dc55d19f585ac04d"
NUM = ("figmech_v39_numbers.json", "a4b23923a65d4284b9580529b70a1845565aa0e5460177e55eb8b85c5442945c")
ARR = ("figmech_v39_arrays.npz", "3ee12fc9a32e285d6c181e2d03a19286ade9ca789b2c2efc2eaa55ea1f27c9d5")
OLD = ("figures/fig_mechanism_numbers.json", "0b7eb8b23db7d341b9c306d0182837c3987fe9f70a5e01fefb35427b1e73976f", "figures/NOT_USED_pre_thrfix_fig_mechanism_numbers.json")
PRE = "figmech_v39/"; DST = "results_src/figmech_v39"


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
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not b5faa828")
        for f, h in (NUM, ARR):
            if sha(inp / f) != h: raise ValueError(f"input sha differs: {f}")
        if sha(tree / OLD[0]) != OLD[1]: raise ValueError("old fig_mechanism_numbers.json sha differs")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        if any(k.startswith(PRE) for k in vals): raise ValueError(f"existing keys already use prefix {PRE}")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    obj = json.load(open(inp / NUM[0])); rel = f"{DST}/{NUM[0]}"; new = {}
    for path, v in leaves(obj, []):
        if isinstance(v, bool) or isinstance(v, (int, float)):
            if isinstance(v, float) and not math.isfinite(v): print("GATE (exit 3): non-finite", path, file=sys.stderr); return 3
            key = PRE + "|".join(path)
            if key in vals or key in new: print("GATE (exit 3): duplicate key", key, file=sys.stderr); return 3
            new[key] = (v, rel, path)
    bak = tree / "results_master_before_figmech_v39_b5faa828.json"; tmp = tree / "results_master.json.figmech_v39_tmp"
    targets = [bak, tmp, tree / DST / NUM[0], tree / DST / ARR[0], tree / OLD[2]]
    existing = [str(x) for x in targets if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}; old_fsha = dict(fsha)
    with open(rm, "rb") as src_fh, open(bak, "xb") as dst_fh: dst_fh.write(src_fh.read())
    (tree / DST).mkdir(parents=True, exist_ok=True)
    for f, h in (NUM, ARR):
        with open(inp / f, "rb") as src_fh, open(tree / DST / f, "xb") as dst_fh: dst_fh.write(src_fh.read())
        if sha(tree / DST / f) != h: print("GATE (exit 3): copy sha differs", f, file=sys.stderr); return 3
        fsha[f"{DST}/{f}"] = h
    for k, (v, rel_, path) in new.items(): vals[k] = v; srcs[k] = dict(file=rel_, path=path)
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp)); src_obj = json.load(open(tree / rel))
    for k, (v, rel_, path) in new.items():
        if get(src_obj, path) != chk["values"][k]: print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    for k, v in old_fsha.items():
        if chk["files_sha256"].get(k) != v: print("GATE (exit 3): an old files_sha256 entry changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new): print("GATE (exit 3): key count mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    os.rename(tree / OLD[0], tree / OLD[2])
    print(json.dumps(dict(added=len(new), n_values=d["n_values"], new_master_sha256=sha(rm), renamed=[OLD[0], OLD[2]]), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

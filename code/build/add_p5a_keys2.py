# P5a (CertBound) ADD-ONLY results_master builder, round 2 (constants for sections 1-4: Z+ alpha / m_p / G, candidate A design constants, Proposition minm+ numbers).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Same mechanics as add_p5a_keys.py v1.1 (312a3391): all write targets checked absent before anything is written; indent=1; after writing every new key is re-read from its source path and every old key/source must be unchanged.
# Unlike round 1, only the LISTED top-level fields are registered (the run_meta / predictions files also hold large per-cell blocks that are not needed): each numeric/boolean leaf under a listed field becomes "<prefix><field>|<p2>|...".
# Inputs (P5a_CertBound/inputs_20261124/, pulled from the VM via the scratchpad; each sha equals the line in its own SHA256SUMS on the VM):
#   prefix          file                        sha256                                                              fields
#   zplus_F1meta/   zF1meta/run_meta.json       eaefcb02dd93981068a2b7529df073596ca029b45a4fc09bcd5a511d7f6508cf    ndraw, mps, alphas         (Z+ F1 run c1zf1_20261123T103940Z)
#   zplus_F2meta/   zF2meta/run_meta.json       67b043827a493757cd1cc5725514f21bff52a833d3c4255648ba0eb1587b2754    ndraw, mps, alphas         (Z+ F2 run c1zrun_20261122T212400Z)
#   zplus_pred/     zpred/predictions.json      c0dbf0f3d7e149fc037334a2062c4f432911ed8ab12abacc9e0cd41cab00a430    G, F, nsite                (sealed Z+ predictions c1zpred_20261121T200746Z)
#   candAmeta/      candAmeta/run_meta.json     9ade90d604a21eb4651256ccf374a4c475c0867dcda1fef6e9511b371432f88a    ndraw, delta_ref, delta_2s, alpha_T, m_T, units_done   (candidate A run c1A_20261123T144758Z)
#   minmplus/       minmplus/minm_plus.json     MINM_SHA (filled after make_minm_plus.py has run and been reported; empty -> exit 3)   all leaves
# Gates (exit 3, nothing written): master sha 58f900ee (round-1 result); input shas; no existing key starts with a new prefix; every listed field exists; every new key new.
# Backup: results_master_before_p5a2_58f900ee.json ('x').  Inputs copied to <tree>/results_src/p5a2/<file> and added to files_sha256.
# Usage: python3 -B add_p5a_keys2.py --tree <tree_bj_v1_...> --inputs <inputs_20261124>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "58f900ee11eb4cd2e909e35a720c2ca9e2fdce7e04215f8bf0683d5f518cc56d"
MINM_SHA = "86a774aa5d34455161106d9e28305e3a2470f0c4216adf5ed968b1efbfff18be"   # v1.2: output of make_minm_plus.py v1.1 b528d919 (rc 0)
SPEC = [("zplus_F1meta/", "zF1meta/run_meta.json", "eaefcb02dd93981068a2b7529df073596ca029b45a4fc09bcd5a511d7f6508cf", ("ndraw", "mps", "alphas")),
        ("zplus_F2meta/", "zF2meta/run_meta.json", "67b043827a493757cd1cc5725514f21bff52a833d3c4255648ba0eb1587b2754", ("ndraw", "mps", "alphas")),
        ("zplus_pred/", "zpred/predictions.json", "c0dbf0f3d7e149fc037334a2062c4f432911ed8ab12abacc9e0cd41cab00a430", ("G", "F", "nsite")),
        ("candAmeta/", "candAmeta/run_meta.json", "9ade90d604a21eb4651256ccf374a4c475c0867dcda1fef6e9511b371432f88a", ("ndraw", "delta_ref", "delta_2s", "alpha_T", "m_T", "units_done")),
        ("minmplus/", "minmplus_v11/minm_plus.json", MINM_SHA, None)]   # None = all leaves


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
        if MINM_SHA == "": raise ValueError("MINM_SHA is empty (make_minm_plus.py not yet run / reported)")
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not the round-1 master 58f900ee")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        for pre, f, h, flds in SPEC:
            if sha(inp / f) != h: raise ValueError(f"input sha differs: {f}")
            if any(k.startswith(pre) for k in vals): raise ValueError(f"existing keys already use prefix {pre}")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    new = {}; per_prefix = {}
    for pre, f, h, flds in SPEC:
        obj = json.load(open(inp / f)); rel = f"results_src/p5a2/{f}"; n = 0
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
    bak = tree / "results_master_before_p5a2_58f900ee.json"; tmp = tree / "results_master.json.p5a2_tmp"
    targets = [bak, tmp] + [tree / "results_src" / "p5a2" / f for pre, f, h, flds in SPEC]
    existing = [str(x) for x in targets if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}
    with open(rm, "rb") as src_fh, open(bak, "xb") as dst_fh: dst_fh.write(src_fh.read())
    for pre, f, h, flds in SPEC:
        dst = tree / "results_src" / "p5a2" / f; dst.parent.mkdir(parents=True, exist_ok=True)
        with open(inp / f, "rb") as src_fh, open(dst, "xb") as dst_fh: dst_fh.write(src_fh.read())
        fsha[f"results_src/p5a2/{f}"] = h
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

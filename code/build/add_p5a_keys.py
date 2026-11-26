# P5a (CertBound) ADD-ONLY results_master builder.  v1.1 (实验管理 review): all write targets (backup, 9 input copies, tmp) are checked to be absent BEFORE anything is written; output written with indent=1 like the original; after writing, all old values and sources must be unchanged key by key.  v1 e01ddc25 kept.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Adds the new, already-read results as keys of results_master.json in the BJ tree, without changing any existing key.
# Inputs (each sha asserted; pulled from the VM through the scratchpad, SHA256SUMS verified, copied to P5a_CertBound/inputs_20261123/):
#   prefix        file (relative to inputs dir)        sha256 (full)
#   zplus_F2/     zF2/reading.json                     9200413e45768f89cc2703f7561d7a53f9f6b7e780bfdb361bc3615f35f1930c   (Z+ F2/F3 frozen reader)
#   zplus_F1/     zF1/reading.json                     3a3f37821ebc91cd95c9eac4e5d04efdbbfe6aae8b842c21ac06ab28548b4ac8   (Z+ F1 frozen reader)
#   c1card/       c1read/reading.json                  88af06e3410cec433013d91e6bf11182fdd345ba721fccd8723e0144ff06ee7b   (card 1 frozen reader)
#   c1root/       c1root/decomposition.json            f276bcfa04576ab43f5329afc2763f5d7bf0881212f735aeb43c998fe53f6e9f   (card 1 root cause T0-T6)
#   m51/          m51read/reading.json                 a748363951c853c4971f92874c6dae0ae8048150d814915d277c166c728fcbda   (M5-1 frozen reader)
#   m51x/         m51x/explore_a2a1.json               61f01ba1c51f23f4ebffd997ceb9992bee12e43253e4bb44570e274215ea7506   (M5-1 post-read exploration, descriptive)
#   candA/        candA/reading.json                   8487c968b7822f80f26ee098b4ee7ae2c6820b83abada3cc36facbe4fbad84ef   (candidate A frozen reader)
#   r1r3v1/       r1r3v1/r1r3.json                     6ee8d0dae7e3971e4ad83454598a66c0500baefa44f8cfd5705ea5a7e158838d   (R1'/R3', the verdict)
#   r1r3v2/       r1r3v2/r1r3.json                     85616d8a67522fa966ccb567d5a815e7db901c0b9fbecbc517a57d9b218717e1   (R1'/R3' + G*, report only)
# Rule: every numeric (int / float, finite) and boolean leaf of each JSON becomes a key "<prefix><p1>|<p2>|...|<pn>" (list indices as decimal strings) with source {file: results_src/p5a/<file>, path: [p1, ..., pn]};
#   strings (shas, names) and None are not registered; non-finite floats are skipped and counted.  The input files are copied into <tree>/results_src/p5a/ and added to files_sha256.
# Gates (any failure -> exit 3, nothing written): tree results_master sha must be 1f6d7a05... (the r4h master, unchanged); no existing key may start with any new prefix; every new key must be new; after writing, every new key is re-read from its source path and must equal the stored value.
# Backup: the original results_master.json is copied to results_master_before_p5a_1f6d7a05.json before it is replaced ('x' mode for all new files).
# Usage: python3 -B add_p5a_keys.py --tree <tree_bj_v1_...> --inputs <inputs_20261123>
import argparse, hashlib, json, math, shutil, sys
from pathlib import Path
MASTER_SHA = "1f6d7a05c3f7b6f79513f6553aa75387f1c568aa9308053e9f2d1ccb8e0d1165"
SPEC = [("zplus_F2/", "zF2/reading.json", "9200413e45768f89cc2703f7561d7a53f9f6b7e780bfdb361bc3615f35f1930c"),
        ("zplus_F1/", "zF1/reading.json", "3a3f37821ebc91cd95c9eac4e5d04efdbbfe6aae8b842c21ac06ab28548b4ac8"),
        ("c1card/", "c1read/reading.json", "88af06e3410cec433013d91e6bf11182fdd345ba721fccd8723e0144ff06ee7b"),
        ("c1root/", "c1root/decomposition.json", "f276bcfa04576ab43f5329afc2763f5d7bf0881212f735aeb43c998fe53f6e9f"),
        ("m51/", "m51read/reading.json", "a748363951c853c4971f92874c6dae0ae8048150d814915d277c166c728fcbda"),
        ("m51x/", "m51x/explore_a2a1.json", "61f01ba1c51f23f4ebffd997ceb9992bee12e43253e4bb44570e274215ea7506"),
        ("candA/", "candA/reading.json", "8487c968b7822f80f26ee098b4ee7ae2c6820b83abada3cc36facbe4fbad84ef"),
        ("r1r3v1/", "r1r3v1/r1r3.json", "6ee8d0dae7e3971e4ad83454598a66c0500baefa44f8cfd5705ea5a7e158838d"),
        ("r1r3v2/", "r1r3v2/r1r3.json", "85616d8a67522fa966ccb567d5a815e7db901c0b9fbecbc517a57d9b218717e1")]


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
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not the r4h master 1f6d7a05")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        for pre, f, h in SPEC:
            if sha(inp / f) != h: raise ValueError(f"input sha differs: {f}")
            if any(k.startswith(pre) for k in vals): raise ValueError(f"existing keys already use prefix {pre}")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    new = {}; skipped_nonfinite = 0; per_prefix = {}
    for pre, f, h in SPEC:
        obj = json.load(open(inp / f)); rel = f"results_src/p5a/{f}"; n = 0
        for path, v in leaves(obj, []):
            if isinstance(v, bool) or (isinstance(v, (int, float)) and not isinstance(v, bool)):
                if isinstance(v, float) and not math.isfinite(v): skipped_nonfinite += 1; continue
                key = pre + "|".join(path)
                if key in vals or key in new: print("GATE (exit 3): duplicate key", key, file=sys.stderr); return 3
                new[key] = (v, rel, path); n += 1
        per_prefix[pre] = n
    bak = tree / "results_master_before_p5a_1f6d7a05.json"; tmp = tree / "results_master.json.p5a_tmp"
    targets = [bak, tmp] + [tree / "results_src" / "p5a" / f for pre, f, h in SPEC]
    existing = [str(x) for x in targets if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}
    with open(rm, "rb") as src_fh, open(bak, "xb") as dst_fh: dst_fh.write(src_fh.read())
    for pre, f, h in SPEC:
        dst = tree / "results_src" / "p5a" / f; dst.parent.mkdir(parents=True, exist_ok=True)
        with open(inp / f, "rb") as src_fh, open(dst, "xb") as dst_fh: dst_fh.write(src_fh.read())
        fsha[f"results_src/p5a/{f}"] = h
    for k, (v, rel, path) in new.items(): vals[k] = v; srcs[k] = dict(file=rel, path=path)
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp)); cache = {}
    for k, (v, rel, path) in new.items():
        if rel not in cache: cache[rel] = json.load(open(tree / rel))
        if get(cache[rel], path) != chk["values"][k]: print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():   # add-only: every old key keeps its value and its source
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new): print("GATE (exit 3): key count mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    print(json.dumps(dict(added=len(new), per_prefix=per_prefix, skipped_nonfinite=skipped_nonfinite, n_values=d["n_values"], new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

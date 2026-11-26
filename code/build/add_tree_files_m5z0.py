# P5a ADD-ONLY tree packaging builder (实验管理 2026-11-25, packaging track (a)): copies two cited result files into the tree and registers them in files_sha256.
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  No value key is added or changed.
# Files (source under Paper5_StainRef/, pinned):
#   results_src/m5/results/m5_confirm.json  55c39915  -> <tree>/results_src/m5/results/m5_confirm.json, new files_sha256 entry (none exists before);
#   results_src/z0/z0.json                  fa9ec2fb  -> <tree>/results_src/z0/z0.json; its files_sha256 entry already exists and must equal the pin (only the file was missing).
# Gates (exit 3, nothing written): master sha MASTER_SHA; both source shas equal the pins; the m5 entry absent from files_sha256 and the z0 entry equal to its pin;
#   both target files absent; backup and tmp absent.  After the copy: copy shas equal the pins.  Then the master is rewritten (indent=1) with only the new m5 entry
#   added: every value, source and other files_sha256 entry unchanged, n_values unchanged.
# Backup: results_master_before_m5z0_<first 8 of MASTER_SHA>.json ('x').
# Usage: python3 -E -B add_tree_files_m5z0.py --tree <tree_bj_v1_...> --p5 <Paper5_StainRef>
import argparse, hashlib, json, sys
from pathlib import Path
MASTER_SHA = "a374574cf7a286e95ee3846aca57b262a87148b52636326a2410b80c3adc6423"
FILES = {"results_src/m5/results/m5_confirm.json": ("55c3991596a7d9fbeaaa0bf234c239f56059929613ca3d66172ee3b8e7690eed", "new"),
         "results_src/z0/z0.json": ("fa9ec2fbae557c556f754b3dc7725fd5d5ae30a1138a4253fe286c6de6e3f797", "existing")}


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); ap.add_argument("--p5", required=True); a = ap.parse_args()
    tree, p5 = Path(a.tree), Path(a.p5); rm = tree / "results_master.json"
    if sha(rm) != MASTER_SHA: return fail("tree results_master differs from the pin")
    d = json.load(open(rm)); fsha = d["files_sha256"]
    for rel, (h, kind) in FILES.items():
        if sha(p5 / rel) != h: return fail("source sha differs from the pin", rel)
        if (tree / rel).exists(): return fail("target exists", rel)
        if kind == "new" and rel in fsha: return fail("files_sha256 entry exists already", rel)
        if kind == "existing" and fsha.get(rel) != h: return fail("existing files_sha256 entry differs from the pin", rel)
    bak = tree / f"results_master_before_m5z0_{MASTER_SHA[:8]}.json"; tmp = tree / "results_master.json.m5z0_tmp"
    if bak.exists() or tmp.exists(): return fail("backup or tmp exists")
    old = json.loads(rm.read_text())
    with open(rm, "rb") as s_, open(bak, "xb") as d_: d_.write(s_.read())
    for rel, (h, kind) in FILES.items():
        (tree / rel).parent.mkdir(parents=True, exist_ok=True)
        with open(p5 / rel, "rb") as s_, open(tree / rel, "xb") as d_: d_.write(s_.read())
        if sha(tree / rel) != h: return fail("copy sha differs", rel)
    fsha["results_src/m5/results/m5_confirm.json"] = FILES["results_src/m5/results/m5_confirm.json"][0]
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp))
    if chk["values"] != old["values"] or chk["sources"] != old["sources"] or chk.get("n_values") != old.get("n_values"): return fail("a value, source or n_values changed")
    if {k: v for k, v in chk["files_sha256"].items() if k != "results_src/m5/results/m5_confirm.json"} != old["files_sha256"]: return fail("an old files_sha256 entry changed")
    tmp.replace(rm)
    print(json.dumps(dict(copied=sorted(FILES), new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

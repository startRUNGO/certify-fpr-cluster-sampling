# P5a bib patch for apalike (实验管理 2026-11-24, approved list): writes ONLY the new bib (mode 'x'); the input is not touched.
# Five field edits, each OLD must occur exactly once; gates: input sha pinned, @key set unchanged, only these five title / pages strings differ.
# Usage: python3 -I -B patch_bib_v63.py --src <references_bj_v39.bib b828897e> --out <new .bib path>
import argparse, hashlib, re, sys
from pathlib import Path
SRC_SHA = "b828897eb1dcff098688e3601ecec976ad2267ad66d831c5f5edf41e6b9f4a20"
EDITS = [
    ("sengar2022vit", "analysis of {\\emph{Plasmodium vivax}} life cycle", "analysis of \\emph{{Plasmodium vivax}} life cycle"),
    ("lacunathin2024", "title        = {Lacuna Malaria Datasets},", "title        = {{Lacuna Malaria Datasets}},"),
    ("bates2021", "volume = {68}, number = {6}, pages = {43},", "volume = {68}, number = {6}, pages = {43:1--43:34},"),
    ("platanios2016agreement", "Estimating Accuracy from Unlabeled Data: A Bayesian Approach", "Estimating Accuracy from Unlabeled Data: A {Bayesian} Approach"),
    ("garg2022atc", "{ATC: An Accurate", "{{ATC}: An Accurate"),
]


def entries(t): return {m.group(2): m.group(0) for m in re.finditer(r"@(\w+)\{([^,]+),.*?\n\}", t, flags=re.S)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--src", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    src, out = Path(a.src), Path(a.out)
    if hashlib.sha256(src.read_bytes()).hexdigest() != SRC_SHA: print("GATE (exit 3): input sha differs", file=sys.stderr); return 3
    if out.exists(): print("GATE (exit 3): output exists", out, file=sys.stderr); return 3
    t = src.read_text(encoding="utf-8"); old = entries(t); new_t = t
    for key, o, n in EDITS:
        e = entries(new_t)[key]
        if e.count(o) != 1 or new_t.count(o) != 1: print(f"GATE (exit 3): {key}: OLD found {e.count(o)} times in the entry, {new_t.count(o)} in the file", file=sys.stderr); return 3
        new_t = new_t.replace(o, n)
    new = entries(new_t)
    if set(new) != set(old): print("GATE (exit 3): @key set changed", sorted(set(new) ^ set(old)), file=sys.stderr); return 3
    changed = sorted(k for k in old if old[k] != new[k])
    if changed != sorted(k for k, _, _ in EDITS): print("GATE (exit 3): changed entries differ from the edit list", changed, file=sys.stderr); return 3
    with open(out, "x", encoding="utf-8") as fh: fh.write(new_t)
    print("written", out, "| entries", len(new), "| changed", changed, "| sha", hashlib.sha256(out.read_bytes()).hexdigest()); return 0


if __name__ == "__main__": sys.exit(main())

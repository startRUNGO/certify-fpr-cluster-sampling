# P5a BJ v4.0 (copy of check_thousands_v39 cf9f5ee3 with three names v39 -> v40: WL_FILE and the two top files; the whitelisted line is unchanged, WL_SHA 44e082fe).
# P5a BJ v3.9 (copy of check_thousands_v38 4154d556; names v39; whitelist line unchanged by v3.9, the SI layout pairs do not touch it, so the same sha 44e082fe): READ-ONLY post-build check that no thousands comma is left in the filled main text and SI.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Files: the recursive \input / \include chain of main_bj_v37.tex and si_main_bj_v37.tex in the work dir (uncommented lines; .tex added when no extension), i.e. the FILLED files that the build compiled.
# Rule: on every uncommented line (first non-blank character not '%'), the regex \d(,|\{,\})\d{3}(?!\d) must not match.
# Whitelist (exact whole-line match only, 实验管理): the one line of sections/si_bj_v37.tex whose sha256 is WL_SHA (the "clip(...-1,0,255)" line of the front-end description; it was line 41 in v36).
# Exit 0 and print the counts when clean; exit 3 and list file:line:match for every other hit.  Writes nothing.
# Usage: python3 -B check_thousands_v37.py --work <bj_v37_...>
import argparse, hashlib, re, sys
from pathlib import Path
WL_FILE = "sections/si_bj_v40.tex"
# v1.1 (实验管理 2026-11-24): the clip formula is now a display line of its own, "\\[ I'_{ic}=...-1,0,255)\\rfloor. \\]"; still exactly one whitelisted line
WL_SHA = "2927bc10917707619d2704add282de588467e81f2bbc535683e15d9bfa51377a"
PAT = re.compile(r"\d(,|\{,\})\d{3}(?!\d)")


def chain(W, rel, seen):
    if rel in seen or not (W / rel).is_file(): return
    seen.add(rel)
    for line in (W / rel).read_text(encoding="utf-8", errors="replace").split("\n"):
        line = re.sub(r"((?:^|[^\\])(?:\\\\)*)%.*$", r"\1", line)
        for m in re.finditer(r"\\(?:input|include)\{([^}]+)\}", line):
            r = m.group(1).strip(); chain(W, r if Path(r).suffix == ".tex" else r + ".tex", seen)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--work", required=True); W = Path(ap.parse_args().work)
    files = set()
    for top in ("main_bj_v40.tex", "si_main_bj_v40.tex"):
        if not (W / top).is_file(): print("GATE (exit 3): missing", top, file=sys.stderr); return 3
        chain(W, top, files)
    hits, wl = [], 0
    for rel in sorted(files):
        for i, line in enumerate((W / rel).read_text(encoding="utf-8", errors="replace").split("\n"), 1):
            if line.lstrip().startswith("%") or not PAT.search(line): continue
            if rel == WL_FILE and hashlib.sha256(line.encode("utf-8")).hexdigest() == WL_SHA: wl += 1; continue
            hits += [f"{rel}:{i}:{m.group(0)}" for m in PAT.finditer(line)]
    if wl != 1: print("GATE (exit 3): whitelisted line found", wl, "times (expected 1)", file=sys.stderr); return 3
    if hits: print("GATE (exit 3): thousands commas left:", *hits, sep="\n  ", file=sys.stderr); return 3
    print("clean:", len(files), "files scanned; whitelisted line matched once"); return 0


if __name__ == "__main__": sys.exit(main())

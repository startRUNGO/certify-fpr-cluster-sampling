# P5a v5.2: READ-ONLY post-build check that the old rule names are gone from the PDF text (实验管理 2026-11-24).  WRITTEN, NOT RUN, until 实验管理 reviews it.  Writes nothing.
# pdftotext of main_bj_v40.pdf and si_main_bj_v40.pdf in --work; whitespace is collapsed and SI pointers "Section(s) / Table(s) / Figure(s) / Fig. S<n>[.<k>]" are removed before the scan.
# Old names searched: S3 as a word, SPLIT6, "rule Z" as words, Z+ (Z$^+$ as rendered: "Z+" or "Z +"), Z-chance / -lto / -c2 / -real / -all.
# Main text: any hit fails (exit 3).  SI: hits are only reported (the three pre-rendered SI figures fig_s_zpi, fig_c1law and fig_c1heat still carry the old names; the SI becomes
#   failing once they are redone, by setting SI_FAILS = True through a reviewed diff).
# v1.2 (实验管理 2026-11-24): SI section headings ('S<n> [line numbers] <title>', numbers and titles read from the .aux) are removed before the scan.
# v1.3 (实验管理 2026-11-24, Wiley USG): the letter-spaced caption labels "T A B L E S<n>.<k>" / "F I G U R E S<n>.<k>" are removed, and an SI section number that
#   stands between two line numbers ("160 S3 161", the USG heading layout) is removed, for every section number read from the .aux.
# Usage: python3 -B check_oldnames_v52.py --work <build dir>
import argparse, re, subprocess, sys
from pathlib import Path

PATS = [r"\bS3\b", r"SPLIT6", r"\brule Z\b", r"\bZ ?\+", r"\bZ-(?:chance|lto|c2|real|all)\b"]
# v1.1 (实验管理 2026-11-24): the SI tables and the three SI figures are renamed, so the SI now fails too
SI_FAILS = True


def headings(aux):
    """(number, first two words of the title) of every SI section / subsection heading, from the .aux table-of-contents lines
    \\contentsline {section}{\\numberline {S<n>[.<k>]}<title>}{<page>}{...}; macros, braces and $ are stripped from the title (Gemini review)."""
    if not aux.is_file(): return []
    out = []
    for num, title in re.findall(r"\\numberline \{(S\d+(?:\.\d+)?)\}(.*?)\}\{\d+\}\{", aux.read_text(encoding="utf-8", errors="replace")):
        w = re.findall(r"[A-Za-z]+", re.sub(r"\\[A-Za-z]+\*?|[{}$]", " ", title))
        if len(w) >= 2: out.append((num, w[0], w[1]))
    return out


def text(pdf, heads=()):
    t = subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, text=True, check=True).stdout
    t = re.sub(r"(?:Sections?|Tables?|Figures?|Fig\.) S\d+(?:\.\d+)?", "", re.sub(r"\s+", " ", t))
    # v1.3 (实验管理): USG prints caption labels letter-spaced
    t = re.sub(r"(?:T A B L E|F I G U R E) S\d+(?:\.\d+)?", "", t)
    # v1.2 (实验管理): SI section headings print as "S<n> [line numbers] <title>" with no "Section" before them
    for num, w1, w2 in heads: t = re.sub(r"\b" + re.escape(num) + r"(?: \d+)* (?=" + re.escape(w1) + r" " + re.escape(w2) + r"\b)", "", t)
    # v1.3 (实验管理): USG heading layout, the section number between two line numbers
    for num in {h[0] for h in heads}: t = re.sub(r"(?<=\d) " + re.escape(num) + r"(?= \d+\b)", " ", t)
    return t


def hits(t):
    return [(p, t[max(0, m.start() - 40):m.end() + 30]) for p in PATS for m in re.finditer(p, t)]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--work", required=True); W = Path(ap.parse_args().work); rc = 0
    for name, fails in (("main_bj_v40.pdf", True), ("si_main_bj_v40.pdf", SI_FAILS)):
        if not (W / name).is_file(): print("GATE (exit 3): missing", name, file=sys.stderr); return 3
        h = hits(text(W / name, headings(W / name.replace(".pdf", ".aux"))))
        print(f"{name}: {len(h)} old-name hit(s)" + ("" if fails else " (report only)"))
        for p, ctx in h: print(f"  {p}: ...{ctx}...")
        if h and fails: rc = 3
    if rc: print("GATE (exit 3): old rule names in a failing document (see the hit list above)", file=sys.stderr)
    return rc


if __name__ == "__main__": sys.exit(main())

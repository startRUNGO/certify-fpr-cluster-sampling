# P5a BJ v3.9 / v4.0: READ-ONLY post-build check that every SI number hard-coded in the filled main text matches the SI aux (实验管理, 2026-11-24).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Reason: the main text cannot \ref into the SI (no xr in the main), so it writes "Section S8", "Table S4.2" literally; the SI numbers come from LaTeX and could shift.
# SI side: <work>/si_main_bj_v39.aux, lines \newlabel{<label>}{{<number>}...}; labels sec:s1 ... sec:s14 and tab:s4_zplus_d must all be present.
# Main side: the recursive \input / \include chain of <work>/main_bj_v39.tex (uncommented text, .tex added when no extension), i.e. the FILLED files the build compiled.
#   Hard-coded references: regex (Sections?|Tables?|Figures?) followed by S<n>[.<k>], including lists "S12 and S13" / "S2, S6".
#   Section S<n> must satisfy aux[sec:s<n>] == "S<n>"; Table S<x> must be "S4.2" and equal aux[tab:s4_zplus_d] (the only SI table the main text cites by number);
#   any other hard-coded table or any SI figure number is a failure (unmapped).
# v1.1: every token matched by ALL must be inside a REF match or listed in EXCL (excluded tokens are printed); otherwise unclassified -> exit 3.
# Exit 0 and list every reference with its aux value when all match; exit 3 and list every mismatch / missing label / unmapped reference.  Writes nothing.
# --tag v38 runs the same check on the v3.8 work dir (validation of the checker on a known-good build).
# v1.2 (assemble v1.24, v4.0): --tag v40; with v40, the SI aux must not define tab:dev or tab:esparse (both tables are in the main text from v4.0 on; the SI cites them through xr).
# v1.3 (实验管理 2026-11-24, assemble v6 / Wiley USG): SI sections S15 (proofs) and S16 (further development-site analyses) are valid targets (sec:s1 ... sec:s16);
#   with v40, the SI aux must not define tab:dev but must define tab:esparse (Table S16.1).
# v1.4 (实验管理 2026-11-25, SI trim v1): after the trim the SI sections are renumbered (S10 -> S9, S12 -> S10, S15 -> S11, S16 -> S12) while the labels keep their names, so a
#   pointer is no longer matched by label name: "Section S<n>" must equal the aux number of some sec: label, "Table S<x>" of some tab: label, "Figure S<x>" of some fig: label
#   (every match is listed with its label).  Required in the SI aux: the labels of the kept sections (sec:s1-s8, s10, s12, s15, s16) and tab:s4_zplus_d; the labels of the deleted
#   sections (sec:s9, s11, s13, s14) must be absent.
# v1.5 (实验管理 2026-11-25, SI trim v2): required sec:s1-s4, s7, s8, s10, s15, s16 and tab:s4_zplus_d; absent sec:s5, s6, s9, s11, s12, s13, s14.
# Usage: python3 -B check_si_numbers_v39.py --work <bj_v39_...> [--tag v38]
import argparse, re, sys
from pathlib import Path
# v1.4: labels of the kept and of the deleted SI sections after the trim
# v1.5 (实验管理 2026-11-25, SI trim v2): S5 (sec:s5), S6 (sec:s6) and S10 (sec:s12) are deleted as well
SECS = [f"sec:s{i}" for i in (1, 2, 3, 4, 7, 8, 10, 15, 16)]; SECS_GONE = [f"sec:s{i}" for i in (5, 6, 9, 11, 12, 13, 14)]; TAB = "tab:s4_zplus_d"; TAB_TEXT = "S4.2"
REF = re.compile(r"\b(Sections?|Tables?|Figures?)(?:~|\s)+(S\d+(?:\.\d+)?(?:(?:\s*,\s*and\s+|\s*,\s*|\s+and\s+)S\d+(?:\.\d+)?)*)")
TOK = re.compile(r"S\d+(?:\.\d+)?")
# v1.1 (实验管理): catch-all for every SI-like token; each must be classified by REF above, or be an EXCLUDED rule name (listed in the output), else "unclassified" -> exit 3
ALL = re.compile(r"(?<![A-Za-z0-9_:{/\\-])S\d+(?:\.\d+)?\b")
EXCL = {"S3": "rule name S3 (not an SI number)"}


def strip(line): return re.sub(r"((?:^|[^\\])(?:\\\\)*)%.*$", r"\1", line)


def chain(W, rel, seen):
    if rel in seen or not (W / rel).is_file(): return
    seen.append(rel)
    for line in (W / rel).read_text(encoding="utf-8", errors="replace").split("\n"):
        for m in re.finditer(r"\\(?:input|include)\{([^}]+)\}", strip(line)):
            r = m.group(1).strip(); chain(W, r if Path(r).suffix == ".tex" else r + ".tex", seen)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--work", required=True); ap.add_argument("--tag", default="v40", choices=("v38", "v39", "v40")); a = ap.parse_args(); W = Path(a.work)
    aux_p, top = W / f"si_main_bj_{a.tag}.aux", f"main_bj_{a.tag}.tex"
    if not aux_p.is_file() or not (W / top).is_file(): print("GATE (exit 3): missing", aux_p.name, "or", top, file=sys.stderr); return 3
    aux = {}
    for line in aux_p.read_text(encoding="utf-8", errors="replace").split("\n"):
        m = re.match(r"\\newlabel\{([^}]+)\}\{\{([^}]*)\}", line)
        if m: aux[m.group(1)] = m.group(2)
    bad = [f"missing label in SI aux: {l}" for l in SECS + [TAB] if l not in aux] + [f"label of a deleted section still in the SI aux: {l}" for l in SECS_GONE if l in aux]
    # v1.4: number -> labels of each kind, from the aux
    bynum = {k: {} for k in ("Section", "Table", "Figure")}
    for l_, v_ in aux.items():
        k_ = {"sec:": "Section", "tab:": "Table", "fig:": "Figure"}.get(l_[:4])
        if k_: bynum[k_].setdefault(v_, []).append(l_)
    # v1.3: tab:dev stays in the main text, tab:esparse is back in the SI (S16)
    if a.tag == "v40": bad += [f"SI aux still defines {l} (table in the main text)" for l in ("tab:dev",) if l in aux] + [f"SI aux does not define {l} (Table S9.1)" for l in ("tab:esparse",) if l not in aux]
    files = []; chain(W, top, files); hits = []; excl = []
    for rel in files:
        text = "\n".join(strip(l) for l in (W / rel).read_text(encoding="utf-8", errors="replace").split("\n"))
        seen = set()
        for m in REF.finditer(text):
            ln = text.count("\n", 0, m.start()) + 1; kind = m.group(1).rstrip("s")
            seen |= {m.start(2) + q.start() for q in TOK.finditer(m.group(2))}
            for t in TOK.findall(m.group(2)):
                # v1.4: matched through the aux numbers of the labels of that kind
                labs = bynum[kind].get(t, [])
                if not labs: bad.append(f"{rel}:{ln}: {kind} {t} is not the number of any SI {kind.lower()} label")
                else: hits.append(f"{rel}:{ln}: {kind} {t} = aux {', '.join(sorted(labs))}")
        for q in ALL.finditer(text):
            if q.start() in seen: continue
            ln = text.count("\n", 0, q.start()) + 1
            if q.group(0) in EXCL: excl.append(f"{rel}:{ln}: {q.group(0)} excluded ({EXCL[q.group(0)]})")
            else: bad.append(f"{rel}:{ln}: unclassified SI-like token {q.group(0)}")
    if bad: print("GATE (exit 3):", *bad, sep="\n  ", file=sys.stderr); return 3
    print(f"clean: {len(files)} files, {len(hits)} hard-coded SI references, all equal to the SI aux; unclassified 0; excluded {len(excl)}", *hits, *excl, sep="\n  "); return 0


if __name__ == "__main__": sys.exit(main())

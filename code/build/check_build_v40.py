# P5a BJ v4.0 (assemble v1.24): READ-ONLY post-build check of the main text and SI build (实验管理, 2026-11-24).  WRITTEN, NOT RUN, until 实验管理 reviews it.  Writes nothing.
# For main_bj_v40 and si_main_bj_v40 (.log and .aux in --work):
# v1.1 (实验管理): multiply defined labels whitelisted: exactly LastPage (w-art + xr, also in v3.9) and RF1 (rotating-package label of a sidewaystable, defined in both documents
#   since Table tab:dev moved to the main text); nobody \ref's them; any other duplicate fails.
#   no LaTeX error line ('! '); no undefined reference or citation (per-item warnings and TeX's summary lines, which do not wrap); no multiply defined label; no 'Float too large' (Figure fig:dev is 166 mm tall);
#   main aux: fig:overview ... fig:confirm numbered 1 ... 10 in the order of FIGS (first-reference order); tab:guarantees, tab:sites, tab:dev, tab:esparse numbered 1 ... 4;
#   SI aux: tab:dev and tab:esparse not defined (the tables are in the main text).
# v1.3 (实验管理 2026-11-24): every 'Overfull \hbox' wider than 5 pt in the main or SI log is listed and fails.
# v1.4 (实验管理 2026-11-24, assemble v6 draft 2, Wiley USG): main aux figures fig:sim ... fig:confirm numbered 1 ... 6 and tables tab:guarantees, tab:sites, tab:dev numbered 1 ... 3;
#   SI aux: tab:dev not defined, tab:esparse defined; fig:minfields, fig:census, fig:esparse, fig:overview = S16.1 ... S16.4; sec:s15 = S15, sec:s16 = S16;
#   every 'Overfull \vbox' fails except the one the USG template itself produces (exactly 50.72pt, reported while \output is active; at most once per document);
#   every plain-text SI pointer in the main-text sources ('Section / Figure / Table S<n>[.<k>]') must be the number of an SI aux label of that kind (sec: / fig: / tab:).
# v1.5 (实验管理 2026-11-24, assemble v6 draft 4): images/Wiley_logo-eps-converted-to.pdf must still have its pinned sha after the build (no Ghostscript conversion overwrote it).
# v1.6 (实验管理 2026-11-25): SI_NUMS after the SI trim (S11 proofs, S12 further analyses: figures S12.1-S12.4, Table S12.1).
# v1.8 (实验管理 2026-11-25): main figures fig:sim ... fig:confirm numbered 1 ... 7 with F-a (fig:funnel) as Figure 4.
# v1.7 (实验管理 2026-11-25): SI_NUMS after the second SI trim (S8 proofs, S9 further analyses: figures S9.1-S9.4, Table S9.1; S5, S6, S7 and their tables).
# v1.9 (实验管理 2026-11-25): main figures fig:sim ... fig:confirm numbered 1 ... 8 with F-c (fig:simbox) as Figure 2; SI figures fig:simheat-s, fig:simbox-s = S3.1, S3.2.
# Exit 0 and print the numbers when all hold; exit 3 and list every failure.
# Usage: python3 -B check_build_v40.py --work <bj_v40_... build dir>
import argparse, hashlib, re, sys
from pathlib import Path
# v1.2 (实验管理, main text v5.0): figure order 1 overview, 2 sim, 3 drift, 4 mechanism, 5 pooled, 6 dev, 7 minfields, 8 census, 9 esparse, 10 confirm
# v1.4: main figures and tables after the S16 move; the S16 figures and S15 / S16 numbers in the SI
# v1.8 (实验管理 2026-11-25, assemble draft 7): F-a (fig:funnel) is Figure 4; F-d and F-e keep the labels fig:pooled and fig:confirm
# v1.9 (实验管理 2026-11-25, assemble draft 9): F-c (fig:simbox) is Figure 2 after F-b (fig:sim)
FIGS = ["fig:sim", "fig:simbox", "fig:drift", "fig:mechanism", "fig:funnel", "fig:pooled", "fig:dev", "fig:confirm"]
# v1.10 (实验管理 2026-11-25, draft 14): the notation table tab:notation is Table 1, so the guarantees map, sites and development tables are Tables 2, 3 and 4.
# v1.11 (实验管理 2026-11-25, draft 15): Table 5 tab:why (baseline families) is new.
TABS = ["tab:notation", "tab:guarantees", "tab:sites", "tab:dev", "tab:why"]
# v1.6 (实验管理 2026-11-25, SI trim v1): S15 -> S11, S16 -> S12; tab:esparse is Table S12.1
# v1.7 (实验管理 2026-11-25, SI trim v2): S7 -> S5, S8 -> S6, S9 -> S7, S11 -> S8, S12 -> S9; tab:esparse is Table S9.1
SI_NUMS = {"fig:minfields": "S9.1", "fig:census": "S9.2", "fig:esparse": "S9.3", "fig:overview": "S9.4", "tab:esparse": "S9.1", "sec:s15": "S8", "sec:s16": "S9",
           "sec:s7": "S5", "sec:s8": "S6", "sec:s10": "S7", "tab:baselines": "S5.1", "tab:s8baselines": "S6.1", "tab:confirm": "S7.1", "tab:confirm-iml": "S7.2", "tab:m5bound": "S7.3",
           "fig:simheat-s": "S3.1", "fig:simbox-s": "S3.2"}
VBOX_WL = "50.72"
# v1.5: the template logo PDF copied by assemble v6 draft 4
LOGO = ("images/Wiley_logo-eps-converted-to.pdf", "166cac545f3f4b2ee5bc27bdfc2686f0d8feba6b9c7e7c7b4f300da8a3174976")
MAIN_SRC = ["sections/bj_main_v40.tex", "sections/bj_abstract_v40.tex", "sections/tab_guar_v50.tex", "sections/tab_sites_v39.tex", "sections/tab_dev_v40.tex", "sections/declarations_v53.tex"]
PTR = r"\b(Section|Figure|Table)s?(?:~|\s)+S(\d+(?:\.\d+)?)\b"
KIND = {"Section": "sec:", "Figure": "fig:", "Table": "tab:"}
WL = {"LastPage", "RF1"}
OVERFULL_PT = 5.0


def aux(p):
    d = {}
    for line in p.read_text(encoding="utf-8", errors="replace").split("\n"):
        m = re.match(r"\\newlabel\{([^}]+)\}\{\{([^}]*)\}", line)
        if m: d[m.group(1)] = m.group(2)
    return d


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--work", required=True); W = Path(ap.parse_args().work); bad = []; info = []
    for d in ("main_bj_v40", "si_main_bj_v40"):
        lp, ap_ = W / f"{d}.log", W / f"{d}.aux"
        if not lp.is_file() or not ap_.is_file(): bad.append(f"missing {d}.log / .aux"); continue
        log = lp.read_text(encoding="utf-8", errors="replace")
        for pat, what in ((r"^! ", "LaTeX error"), (r"(Reference|Citation) `[^']*' on page \d+ undefined", "undefined reference / citation"),
                          (r"Float too large", "Float too large"), (r"There were undefined (references|citations)", "undefined summary line")):
            n = len(re.findall(pat, log, flags=re.M))
            (bad if n else info).append(f"{d}: {what}: {n}")
        # multiply defined labels: join the 79-column wrapped log lines, then every duplicated label must be in WL
        j = "".join(l if len(l) == 79 else l + "\n" for l in log.split("\n"))
        dup = sorted(set(re.findall(r"Label `([^']*)' multiply defined", j)))
        # v1.3 (实验管理): every Overfull \hbox wider than OVERFULL_PT is listed and fails (main and SI)
        # Gemini review: also the line-less form "... has occurred while \output is active"
        ov = [(float(x), rest.strip()) for x, rest in re.findall(r"Overfull \\hbox \(([\d.]+)pt too wide\)([^\n]*)", j)]
        big = [f"{w:.2f}pt {rest[:60]}" for w, rest in ov if w > OVERFULL_PT]
        (bad if big else info).append(f"{d}: Overfull \\hbox > {OVERFULL_PT}pt: {big if big else 0}")
        # v1.4 (实验管理): every Overfull \vbox fails, except one template vbox of exactly 50.72pt reported while \output is active
        vb = re.findall(r"Overfull \\vbox \(([\d.]+)pt too high\)([^\n]*)", j)
        wl = [x for x in vb if f"{float(x[0]):.2f}" == VBOX_WL and "has occurred while \\output is active" in x[1]]
        other = [f"{x[0]}pt {x[1].strip()[:60]}" for x in vb if x not in wl]
        if other or len(wl) > 1: bad.append(f"{d}: Overfull \\vbox outside the whitelist: {other}, template vbox count {len(wl)}")
        else: info.append(f"{d}: Overfull \\vbox: template {VBOX_WL}pt x {len(wl)}, other 0")
        if set(dup) - WL: bad.append(f"{d}: multiply defined labels outside the whitelist: {sorted(set(dup) - WL)}")
        else: info.append(f"{d}: multiply defined labels (all whitelisted): {dup}")
    if (W / "main_bj_v40.aux").is_file():
        A = aux(W / "main_bj_v40.aux")
        for i, l in enumerate(FIGS, 1):
            (info if A.get(l) == str(i) else bad).append(f"main aux {l} = {A.get(l)} (expected {i})")
        for i, l in enumerate(TABS, 1):
            (info if A.get(l) == str(i) else bad).append(f"main aux {l} = {A.get(l)} (expected {i})")
    if (W / "si_main_bj_v40.aux").is_file():
        S = aux(W / "si_main_bj_v40.aux")
        (bad if "tab:dev" in S else info).append(f"SI aux defines tab:dev: {'tab:dev' in S}")
        (info if "tab:esparse" in S else bad).append(f"SI aux defines tab:esparse: {'tab:esparse' in S}")
        for l, n in SI_NUMS.items(): (info if S.get(l) == n else bad).append(f"SI aux {l} = {S.get(l)} (expected {n})")
        # v1.4: plain-text SI pointers in the main-text sources must name an SI number of the same kind
        nums = {k: {v for l, v in S.items() if l.startswith(p)} for k, p in KIND.items()}
        for rel in MAIN_SRC:
            if not (W / rel).is_file(): bad.append(f"missing {rel}"); continue
            src = (W / rel).read_text(encoding="utf-8", errors="replace")
            for kind, num in re.findall(PTR, src):
                (info if "S" + num in nums[kind] else bad).append(f"pointer {kind} S{num} in {rel}: {'S' + num in nums[kind]}")
    lp = W / LOGO[0]
    lsha = hashlib.sha256(lp.read_bytes()).hexdigest() if lp.is_file() else None
    (info if lsha == LOGO[1] else bad).append(f"{LOGO[0]} sha after build: {lsha} (expected {LOGO[1][:8]})")
    if bad: print("GATE (exit 3):", *bad, sep="\n  ", file=sys.stderr); return 3
    print("clean:", *info, sep="\n  "); return 0


if __name__ == "__main__": sys.exit(main())

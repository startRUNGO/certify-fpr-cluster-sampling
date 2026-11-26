# P5a main-text Table 4: data sets and confirmation sites (one row per site).  WRITTEN, NOT RUN, until 实验管理 reviews it.  Writes a LaTeX SOURCE with [[key :: fmt]] placeholders (filled later by tools/fill_numbers.py).
# Columns and keys (every key already in the master unless marked NEW):
#   Fields F      dataset/<site>|fields_all (d); M5: m5/n_fields                                  [results_src/*: e8_plasmo_counts, r8 counts.json, m5_confirm.json 55c39915]
#   Cells         dataset/<site>|cells (d); M5: m5/n_crops
#   Negatives     dataset/<site>|negatives (d); M5: none registered -> '--'
#   Prevalence    NEW p5aprev/<site> (f3) = positives / cells (make_site_prev.py; equals rev5 pi_site where that exists); M5 -> '--'
#   Label source  text, NOT a key: [[MISSING: label source of <site> (writing session)]] (drafts in the review message)
#   Confirmation sites only (development rows '--'), copied from the reviewed tab_confirm4.tex (b30dad38, the P5 confirmation table):
#     PlasmoCount  criterion 'Cov. >= facts/e8_success_cov_min; DeltaSens^a >= facts/e8_success_dsens_min'  Cov. e8_plasmo_results/a2|plasmo|22|S3|0.1|cov (f3)  DeltaSens e8_plasmo_results/a2|plasmo|22|S3|0.1|ds_c0ref05 (sf3)  verdict 'met'
#     IML Malaria  criterion 'Cov. >= ...; DeltaSens_eps >= ... or >= 1/2 DeltaSens_eps(R1-2)'               Cov. r8_rev1_results/agg|a2|iml|22|S3|0.1|cov (f3)         DeltaSens g3/iml|a2|S3|0.1|dsens|point (sf3)                verdict 'met'
#     M5 (SPLIT6)  criterion 'Cov. >= m5/criterion|coverage; DeltaSens_own >= S3 + m5/criterion|gain; no flag'  Cov. m5/judgement|coverage (f3)                         DeltaSens_own - S3: m5/judgement|ds_own_diff_vs_S3 (sf3)   verdict 'confirmed'
#   The verdict words are the reviewed text of tab_confirm4 (they are not keys).  All at target FPR <= protocol/eps_main.
# Gates (exit 3, nothing written): master sha (the one after add_p5a_keys4; filled in when that master exists, empty -> exit 3); tab_confirm4.tex sha b30dad38 and the three criterion strings found verbatim in it;
#   every placeholder key in the master.  Output (new dir, 'x'): tab_sites.tex (sidewaystable from v1.4, label tab:sites, caption [[MISSING: caption of Table 4 (writing session)]]).
# v1.4 (实验管理 v3.1 layout): sidewaystable; Labels and criterion columns p{3.8cm} and p{6.5cm} (v1.3 broke words in 0.14/0.22 textwidth columns).
# Usage: python3 -B make_tab_sites.py --tree <tree_bj_v1_...> --out <new dir>
import argparse, hashlib, json, re, sys
from pathlib import Path
MASTER_SHA = "b5faa82817383bc977e03617e6e85e53ad15b9347a079d06dc55d19f585ac04d"   # v1.2: master after add_p5a_keys4 v1.3
CONFIRM4_SHA = "b30dad38a14e044445986f4154c6f51394696d4e5450dd5a7c31a4c680b247be"   # v1.1 (Gemini): full sha
P = lambda k, f: f"[[{k} :: {f}]]"
CRIT = {"plasmo": r"Cov.\ $\ge[[facts/e8_success_cov_min :: f2]]$; $\Delta$Sens$^{a} \ge[[facts/e8_success_dsens_min :: f2]]$",
        "iml": r"Cov.\ $\ge[[facts/e8_success_cov_min :: f2]]$; $\Delta$Sens$_\epsilon \ge[[facts/e8_success_dsens_min :: f2]]$ or $\ge\frac12\Delta$Sens$_\epsilon$(R1-2)",
        "m5": r"Cov.\ $\ge[[m5/criterion|coverage :: f2]]$; $\Delta$Sens$_{\mathrm{own}}\ge$ S3 $+[[m5/criterion|gain :: f2]]$; no flag"}
ROWS = [("BBBC-ALL", "bbbc041_all", None), ("MP-IDB", "mpidb_full", None), ("Lacuna-B", "lacuna_B", None),
        ("PlasmoCount", "plasmo", (P("e8_plasmo_results/a2|plasmo|22|S3|0.1|cov", "f3"), P("e8_plasmo_results/a2|plasmo|22|S3|0.1|ds_c0ref05", "sf3"), "met")),
        ("IML Malaria", "iml", (P("r8_rev1_results/agg|a2|iml|22|S3|0.1|cov", "f3"), P("g3/iml|a2|S3|0.1|dsens|point", "sf3"), "met")),
        ("M5 (SPLIT6)", "m5", (P("m5/judgement|coverage", "f3"), P("m5/judgement|ds_own_diff_vs_S3", "sf3"), "confirmed"))]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(); T, out = Path(a.tree), Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    if MASTER_SHA == "" or hashlib.sha256((T / "results_master.json").read_bytes()).hexdigest() != MASTER_SHA: print("GATE (exit 3): master sha empty or differs", file=sys.stderr); return 3
    c4 = (T / "sections_src/tab_confirm4.tex").read_bytes()
    if hashlib.sha256(c4).hexdigest() != CONFIRM4_SHA or any(c not in c4.decode() for c in CRIT.values()): print("GATE (exit 3): tab_confirm4 pin or criterion text", file=sys.stderr); return 3
    v = json.loads((T / "results_master.json").read_text())["values"]; rows = []
    for name, s, conf in ROWS:
        if s == "m5": f_, c_, n_, pr = P("m5/n_fields", "d"), P("m5/n_crops", "d"), "--", "--"
        else: f_, c_, n_, pr = P(f"dataset/{s}|fields_all", "d"), P(f"dataset/{s}|cells", "d"), P(f"dataset/{s}|negatives", "d"), P(f"p5aprev/{s}", "f3")
        lab = f"[[MISSING: label source of {name} (writing session)]]"
        tail = f"{CRIT[s]} & {conf[0]} & {conf[1]} & {conf[2]}" if conf else "-- & -- & -- & --"
        rows.append(f"{name} & {f_} & {c_} & {n_} & {pr} & {lab} & {tail} \\\\")
        if name == "Lacuna-B": rows.append(r"\midrule")
    T4 = [r"\begin{sidewaystable}", r"\centering", r"\footnotesize", r"\caption{[[MISSING: caption of Table 4 (writing session)]]}", r"\label{tab:sites}", r"\setlength{\tabcolsep}{3pt}",
          r"\begin{tabular}{@{}lrrrrp{3.8cm}p{6.5cm}rrl@{}}", r"\toprule",
          r" & & & & & & \multicolumn{4}{c}{Confirmation, target FPR $\le " + P("protocol/eps_main", "f2") + r"$} \\", r"\cmidrule(l){7-10}",
          r"Site & Fields & Cells & Negatives & Prevalence & Labels & Sealed criterion & Cov. & $\Delta$Sens & Verdict \\", r"\midrule"] + rows + [r"\bottomrule", r"\end{tabular}", r"\end{sidewaystable}", ""]
    bad = sorted({k for k in re.findall(r"\[\[\s*([^\[\]]+?)\s*::\s*\w+\s*\]\]", "\n".join(T4)) if k not in v})   # v1.3: [^[]] so a [[MISSING: ...]] marker cannot start a match
    if bad: print("GATE (exit 3): keys not in master:", bad, file=sys.stderr); return 3
    out.mkdir(parents=True)
    with open(out / "tab_sites.tex", "x") as fh: fh.write("\n".join(T4))
    print(hashlib.sha256((out / "tab_sites.tex").read_bytes()).hexdigest()); return 0


if __name__ == "__main__": sys.exit(main())

# P5a SI tables S4, S5, S6, S11 as LaTeX SOURCES with [[key :: fmt]] placeholders, laid out from the spec neon_SI_new_20261123T210103Z/SI_tables_spec_20261123T210510Z.md (00adfd97).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Only layout is decided here; every number is a placeholder (filled later by tools/fill_numbers.py); captions are [[MISSING: caption ... (writing session)]] because the spec gives none; column / row wording is
#   taken from the spec's row and column names (shortened only where the spec gives a long phrase; the full list is in the review message).  Booleans (holds / passes) use fmt raw (prints True / False).
# Gates (exit 3, nothing written): master sha f03855aa; spec sha 00adfd97; every placeholder key exists in the master, except the two p5aconst/candA|strata_cut keys (registered by add_p5a_keys4 after review);
#   exactly 8 M4 cells; the one key known to be absent, m51/result|M1|cells|lacuna_B|0.1|b1u_coverage (design-effect fallback share 0 there), is printed as '--'.
# Output (new dir, 'x'): tab_s4_zplus.tex, tab_s5_variant.tex, tab_s6_unequal.tex, tab_s11_mstar.tex.
# v1.2 (实验管理): captions from the writing session's file (e059f480, --captions); header names as named in those captions (S4 A/B, S5 A/B, S6 A two columns, S11); caption placeholders without underscores; gate regex [^[]]; f3 placeholders whose value has 0 < |x| < 0.0005 become f4 (the -0.000 in S6 B2).
# v1.3 (实验管理): S6 caption from 40ab5685 (--s6caption); S6 headers and panel titles A0/A1/A2 -> B/N/U (keys untouched), group head 'Design-effect bound selected'; long S4 A / S5 A / S6 A, B2 header cells broken over two lines.
# v1.5: S4D caption from the pinned delivery file 0b1357e3 (section 2), --s4dcaption_sha removed.
# v1.4 (实验管理): S4 panel D (30 Z+ development-site combinations) as a fifth output tab_s4_zplus_d.tex (label tab:s4_zplus_d), caption from --s4dcaption (+ --s4dcaption_sha).
# Usage: python3 -B make_si_tables.py --master <tree>/results_master.json --spec <SI_tables_spec...md> --captions <交P5_正文sha与全部表注...md> --s6caption <S6表注与tab_baselines引用...md> --s4dcaption <交付_S4D表注_AI声明_SI记号...md> --out <new dir>
import argparse, hashlib, json, re, sys
from pathlib import Path
MASTER_SHA = "b5faa82817383bc977e03617e6e85e53ad15b9347a079d06dc55d19f585ac04d"; SPEC_SHA = "00adfd97"
CAP_SHA = "e059f48050ea23581d8481ddbc5e85a0dd60c944c61cbfee2f2e258b71ee6454"   # v1.2: captions file of the writing session (交P5_正文sha与全部表注, sections TABS4 / TABS5 / TABS6 / TABS11)
S6CAP = "S6表注与tab_baselines引用_写作指导_20261123T215058Z.md"   # v1.3: the S6 caption is taken from this file (40ab5685), not from TABS6 of e059f480
S6CAP_SHA = "40ab5685d7ac1bf946ad272e874324dca4d172e17da39a7ef64a097c08b3fd0f"
S4DCAP_SHA = "0b1357e3ed1a3bf47212b02349c5a63ae2ca8a3422a5a5b4ff4291bf996f8c6e"   # v1.5: delivery file of the writing session (S4D caption in section 2)
CAPSEC = {"tab_s4_zplus.tex": ("S4", "TABS4"), "tab_s5_variant.tex": ("S5", "TABS5"), "tab_s6_unequal.tex": ("S6", "TABS6"), "tab_s11_mstar.tex": ("S11", "TABS11")}
PENDING = {"p5aconst/candA|strata_cut|0", "p5aconst/candA|strata_cut|1"}; ABSENT = {"m51/result|M1|cells|lacuna_B|0.1|b1u_coverage"}
SITES = (("bbbc041_all", "BBBC-ALL"), ("mpidb_full", "MP-IDB"), ("lacuna_B", "Lacuna-B")); EPS = (("0.1", "protocol/eps_main"), ("0.05", "protocol/eps_secondary"))
P = lambda k, f: f"[[{k} :: {f}]]"
# v1.3: a header cell broken over lines (nested tabular, no extra package)
H = lambda *parts: r"\begin{tabular}[b]{@{}c@{}}" + r"\\".join(parts) + r"\end{tabular}"


def head(label, ncol_spec):   # v1.2: caption placeholder names the table (S4, S5, ...) without the label, so it has no underscore
    return [r"\begin{table}[!htbp]", r"\centering", r"\footnotesize", r"\caption{[[MISSING: caption of SI table " + label.split(":")[1].split("_")[0].upper() + r" (writing session)]]}", r"\label{" + label + "}", r"\setlength{\tabcolsep}{3pt}"]


def s4():
    F1 = "zplus_F1/result|"; F2 = "zplus_F2/result|"; mp = ("22", "40"); rows = []
    strata = [(f"Pooled, $\\alpha={P('zplus_F1meta/alphas|0','f2')}$", F1 + "F1_main_alpha005_pooled|", None), (f"Pooled, $\\alpha={P('zplus_F1meta/alphas|1','f2')}$", F1 + "alpha010_pooled|", None)]
    strata += [(f"$m_p={P(f'zplus_F1meta/mps|{i}','d')}$", F1 + f"per_mp_alpha005|{m}|", F1 + f"mean_clicks|{m}") for i, m in enumerate(mp)]
    strata += [(f"{n}, $\\epsilon={P(e[1],'f2')}$", F1 + f"per_site_eps_alpha005|{s}|{e[0]}|", None) for s, n in SITES for e in EPS]
    for lab, pre, mc in strata:
        rows.append(f"{lab} & {P(pre+'n_draws','d')} & {P(pre+'n_accepted','d')} & {P(pre+'acceptance_rate','f3')} & {P(pre+'coverage_among_accepted','f3')} & {P(mc,'f2') if mc else '--'} \\\\")
    A = [r"\begin{tabular}{@{}lrrrrr@{}}", r"\toprule", r"\multicolumn{6}{@{}l}{\emph{A. Coverage among certified draws (development sites)}} \\", r"Stratum & Draws & " + H("Certified", "draws") + " & " + H("Certification", "rate") + " & " + H("Coverage among", "certified draws") + " & " + H("Mean annotated", "positives") + r" \\", r"\midrule"] + rows + [r"\midrule",
         f"\\multicolumn{{6}}{{@{{}}l}}{{Gate on the pooled row at $\\alpha={P('zplus_F1meta/alphas|0','f2')}$: coverage $\\ge {P(F1+'F1_main_alpha005_pooled|gate_ge','f2')}$}} \\\\", r"\bottomrule", r"\end{tabular}"]
    b = [("False certification, target cells, pooled", f"{P(F2+'F2a|n_target_cells','d')} / {P(F2+'F2a|n_draws','d')}", P(F2 + "F2a|false_cert_rate_alpha005_pooled", "f3"), f"{P(F2+'F2a|upper95_cluster_bootstrap_over_cells','f3')}; gate $\\le {P(F2+'F2a|gate_le','f2')}$")]
    b += [(f"\\quad $m_p={P(f'zplus_F1meta/mps|{i}','d')}$", "--", P(F2 + f"F2a|per_mp|{m}", "f3"), "--") for i, m in enumerate(mp)]
    b += [(f"\\quad $\\alpha={P('zplus_F1meta/alphas|1','f2')}$", "--", P(F2 + "F2a|alpha010_pooled", "f3"), "--"),
          ("Coverage among certified target draws", P(F2 + "coverage_among_accepted_target_alpha005|n_accepted", "d"), P(F2 + "coverage_among_accepted_target_alpha005|coverage", "f3"), "--"),
          ("Certified draws with FPR above $\\epsilon$ (report only)", "--", P(F2 + "F2a|accepted_and_fpr_gt_eps_report_only", "f3"), "--"),
          ("Rejection, clear cells, pooled", f"{P(F2+'F2b|n_clear_cells','d')} / {P(F2+'F2b|n_draws','d')}", P(F2 + "F2b|rejection_rate_alpha005_pooled", "f3"), f"gate $\\ge {P(F2+'F2b|gate_ge','f2')}$")]
    b += [(f"\\quad $m_p={P(f'zplus_F1meta/mps|{i}','d')}$", "--", P(F2 + f"F2b|per_mp|{m}", "f3"), "--") for i, m in enumerate(mp)]
    b += [(f"Rejection, $m_p={P(f'zplus_F1meta/mps|{i}','d')}$, $\\alpha={P(f'zplus_F1meta/alphas|{j}','f2')}$", "--", P(F2 + f"F2c|{m}|{al}", "f3"), "--") for i, m in enumerate(mp) for j, al in enumerate(("0.05", "0.1"))]
    B = [r"\begin{tabular}{@{}lrrr@{}}", r"\toprule", r"\multicolumn{4}{@{}l}{\emph{B. False certification and rejection (simulation)}} \\", r"Criterion or check & Cells or draws & Estimate & Upper limit or gate \\", r"\midrule"] + [" & ".join(r) + r" \\" for r in b] + [r"\bottomrule", r"\end{tabular}"]
    c = [f"$m_p={P(f'zplus_F1meta/mps|{i}','d')}$ & {P(F2+f'F3|{m}|n_cells','d')} & {P(F2+f'F3|{m}|mean_clicks','f2')} & {P(F2+f'F3|{m}|ratio_same_fields','f3')} & {P(F2+f'F3|{m}|ratio_to_22_fields','f3')} & "
         f"{P(F2+f'F3|{m}|by_pi_same_fields|0.03','f3')} & {P(F2+f'F3|{m}|by_pi_same_fields|0.08','f3')} & {P(F2+f'F3|{m}|by_pi_to_22|0.03','f3')} & {P(F2+f'F3|{m}|by_pi_to_22|0.08','f3')} \\\\" for i, m in enumerate(mp)]
    C = [r"\begin{tabular}{@{}lrrrrrrrr@{}}", r"\toprule", r"\multicolumn{9}{@{}l}{\emph{C. Annotation burden (report only)}} \\", r" & & & \multicolumn{2}{c}{Ratio} & \multicolumn{2}{c}{Ratio, same fields} & \multicolumn{2}{c}{Ratio, to 22 fields} \\",
         r"\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}", r"Field count & Cells & Mean positives & Same fields & To 22 fields & $\pi=0.03$ & $\pi=0.08$ & $\pi=0.03$ & $\pi=0.08$ \\", r"\midrule"] + c + [r"\bottomrule", r"\end{tabular}"]
    return head("tab:s4_zplus", None) + A + [r"\par\medskip"] + B + [r"\par\medskip"] + C + [r"\end{table}", ""]


def s5():
    a = "candA/applicable_summary|"; C = "candA/"
    crit = [("Marginal coverage", C + "C1|a|value", C + "C1|a|gate", C + "C1|a|holds", P(a + "marginal_mcse", "f3")), ("Joint coverage", C + "C1|b|value", C + "C1|b|gate", C + "C1|b|holds", P(a + "joint_mcse", "f3")),
            ("Share of cells meeting the cell-level gate", C + "C1|c|share_cells_ok", C + "C1|c|gate", C + "C1|c|holds", P(C + "C1|c|n_cells", "d") + " cells"),
            ("Share of draws without a limit", C + "C1|cap|no_limit_share", C + "C1|cap|gate", C + "C1|cap|holds", "--"), ("Mean gap", C + "C2|mean_gap", C + "C2|gate", C + "C2|holds", "--"),
            ("Certification drop", C + "C3|cert_drop", C + "C3|gate", C + "C3|holds", "--"), ("Sensitivity drop", C + "C3|sens_drop", C + "C3|gate", C + "C3|holds", "--")]
    rows = [f"{q} & {P(v,'f3')} & {mc} & {P(g,'f2')} & {P(h,'raw')} \\\\" for q, v, g, h, mc in crit]
    rep = [("Marginal coverage, certified draws", "marginal_certified_only"), ("Joint coverage, certified draws", "joint_certified_only"), ("Joint coverage given valid FPR", "joint_given_fpr_ok"),
           ("Coverage on $L_c$", "coverage_same_Lc"), ("Coverage on $L_s$", "coverage_Ls")]
    rows += [r"\midrule", r"\multicolumn{5}{@{}l}{\emph{Report only}} \\"] + [f"{q} & {P(a+k,'f3')} & -- & -- & -- \\\\" for q, k in rep]
    rows += [f"Clopper--Pearson coverage where applicable (no guarantee) & {P(a+'coverage_conservative_where_applicable','f3')} & {P(a+'n_conservative_applicable','d')} draws & -- & -- \\\\",
             f"Share of draws with the conservative limit at SE $=0$ & {P(a+'se0_conservative_share','f3')} & -- & -- & -- \\\\",
             f"Certification rate, reference / two-sided & {P(a+'cert_ref','f3')} / {P(a+'cert_2s','f3')} & -- & -- & -- \\\\",
             f"Sensitivity, reference / two-sided & {P(a+'sens_ref','f3')} / {P(a+'sens_2s','f3')} & -- & -- & -- \\\\", f"Draws & {P(a+'n','d')} & -- & -- & -- \\\\"]
    A = [r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule", r"\multicolumn{5}{@{}l}{\emph{A. Preregistered criteria and coverage summaries}} \\", r"Quantity & Value & " + H("Monte Carlo SE", "or count") + r" & Gate & Outcome \\", r"\midrule"] + rows + [r"\bottomrule", r"\end{tabular}"]
    st = [(f"At or below {P('p5aconst/candA|strata_cut|0','f1')}", "tpr<=0.5"), (f"Between {P('p5aconst/candA|strata_cut|0','f1')} and {P('p5aconst/candA|strata_cut|1','f1')}", "0.5<tpr<=0.8"), (f"Above {P('p5aconst/candA|strata_cut|1','f1')}", "tpr>0.8")]
    B = [r"\begin{tabular}{@{}lrr@{}}", r"\toprule", r"\multicolumn{3}{@{}l}{\emph{B. Mean gap by true TPR stratum}} \\", r"True TPR stratum & Draws & Mean gap \\", r"\midrule"] + [f"{l} & {P(C+f'C2|by_true_tpr|{k}|n','d')} & {P(C+f'C2|by_true_tpr|{k}|mean_gap','f3')} \\\\" for l, k in st] + [r"\bottomrule", r"\end{tabular}"]
    return head("tab:s5_variant", None) + A + [r"\par\medskip"] + B + [r"\end{table}", ""]


def s6(m4cells):
    M = "m51/result|"; r = []
    for s, n in SITES:
        for e, ek in EPS:
            pre = M + f"M1|cells|{s}|{e}|"; b1 = pre + "b1u_coverage"
            r.append(f"{n}, $\\epsilon={P(ek,'f2')}$ & {P(pre+'n_K1','d')} & {P(pre+'A0_K1_coverage','f3')} & {P(pre+'A1_K1_coverage','f3')} & {P(pre+'A2_K1_coverage','f3')} & {P(pre+'b1u_share','f3')} & {'--' if b1 in ABSENT else P(b1,'f3')} \\\\")
    A = [r"\begin{tabular}{@{}lrrrrrr@{}}", r"\toprule", r"\multicolumn{7}{@{}l}{\emph{A. Per-cell validity (development sites); gate on U: coverage $\ge " + P(M + "M1|gate_ge", "f2") + r"$}} \\",
         r" & & \multicolumn{3}{c}{Coverage, difference-bound draws} & \multicolumn{2}{c}{Design-effect bound selected} \\", r"\cmidrule(lr){3-5}\cmidrule(lr){6-7}",
         r"Site and target & Draws & B & N & U & " + H("Share selecting", "design-effect bound") + " & " + H("Coverage on", "those draws") + r" \\", r"\midrule"] + r + [r"\bottomrule", r"\end{tabular}"]
    b1 = [f"{n} & {P(M+f'M2|sites|{s}|n_pairs','d')} & {P(M+f'M2|sites|{s}|mean_hw_A0','f3')} & {P(M+f'M2|sites|{s}|mean_hw_A2','f3')} & {P(M+f'M2|sites|{s}|ratio_A2_over_A0','f3')} & {P(M+f'M2|sites|{s}|ratio_A1_over_A0_secondary','f3')} & {P(M+f'M2|sites|{s}|passes','raw')} \\\\" for s, n in SITES]
    b2 = [f"{n} & {P(M+f'M3|sites|{s}|mean_diff','f3')} & {P(M+f'M3|sites|{s}|lower95_seed_bootstrap','f3')} & {P(M+f'M3|sites|{s}|K1_branch_only_mean_diff','f3')} & {P(M+f'M3|sites|{s}|A1_minus_A0_mean_diff','f3')} & {P(M+f'M3|sites|{s}|passes','raw')} \\\\" for s, n in SITES]
    B = [r"\begin{tabular}{@{}lrrrrrr@{}}", r"\toprule", r"\multicolumn{7}{@{}l}{\emph{B1. Half-width at the reference threshold; gate: U/B $\le " + P("m51gate/M2|hw_ratio_max", "f2") + r"$}} \\",
         r"Site & Pairs & B & U & U/B & N/B & Met \\", r"\midrule"] + b1 + [r"\midrule", r"\multicolumn{7}{@{}l}{\emph{B2. Sensitivity gain; gate: U$-$B $\ge " + P("m51gate/M3|sens_diff_min", "f2") + r"$ and lower limit $>0$}} \\",
         r"Site & U$-$B & Lower limit & " + H("U$-$B, difference-", "bound draws") + r" & N$-$B & Met & \\", r"\midrule"] + [x.replace(r" \\", r" & \\") for x in b2] + [r"\bottomrule", r"\end{tabular}"]
    c = []
    for cell in m4cells:
        m = re.match(r"^main_pi([0-9.]+)_xi(lo|mid|hi)_cv([0-9.]+)_eps([0-9.]+)_m22$", cell); pre = M + f"M4_report_only|cells|{cell}|"
        c.append(f"{m.group(1)} & {m.group(3)} & {m.group(4)} & {P(pre+'eps','f2')} & {P(pre+'n_K1','d')} & {P(pre+'A2_K1_coverage','f3')} & {P(pre+'hw_ratio_A2_over_A0','f3')} & {P(pre+'passes','raw')} \\\\")
    C = [r"\begin{tabular}{@{}rrrrrrrr@{}}", r"\toprule", r"\multicolumn{8}{@{}l}{\emph{C. Permutation check across synthetic cells (report only)}} \\",
         r"$\pi$ & CV & Cell target & $\epsilon$ evaluated & Draws & U coverage & U/B & Met \\", r"\midrule"] + c + [r"\bottomrule", r"\end{tabular}"]
    return head("tab:s6_unequal", None) + A + [r"\par\medskip"] + B + [r"\par\medskip"] + C + [r"\end{table}", ""]


def s4d(v):   # v1.4 (实验管理): S4 panel D as its own float: the 30 development-site combinations of Z+ (zplus_F1/result|combinations|<i>), site identified by its field count F
    F2s = {v[f"dataset/{s}|fields_all"]: n for s, n in SITES}
    if len(F2s) != 3: raise ValueError("site field counts not distinct")
    rows = []
    for i in range(30):
        pre = f"zplus_F1/result|combinations|{i}|"
        if v.get(pre + "F") not in F2s: raise ValueError(f"combination {i}: F does not identify a site")
        rows.append(f"{F2s[v[pre + 'F']]} & {P(pre+'seed','d')} & {P(pre+'eps','f2')} & {P(pre+'sens','f3')} & {P(pre+'fpr','f3')} & {P(pre+'accepted_22','d')} & {P(pre+'accepted_40','d')} \\\\")
    if f"zplus_F1/result|combinations|30|F" in v: raise ValueError("more than 30 combinations")
    return [r"\begin{table}[!htbp]", r"\centering", r"\footnotesize", r"\caption{[[MISSING: caption of SI table S4D (writing session)]]}", r"\label{tab:s4_zplus_d}", r"\setlength{\tabcolsep}{3pt}",
            r"\begin{tabular}{@{}lrrrrrr@{}}", r"\toprule", r"\multicolumn{7}{@{}l}{\emph{D. Development-site combinations at the census threshold}} \\",
            r" & & & & & \multicolumn{2}{c}{Certified draws} \\", r"\cmidrule(lr){6-7}",
            r"Site & Seed & $\epsilon$ & " + H("Sensitivity", r"at $\tau^+$") + " & " + H("FPR", r"at $\tau^+$") + f" & $m_p={P('zplus_F1meta/mps|0','d')}$ & $m_p={P('zplus_F1meta/mps|1','d')}$ \\\\", r"\midrule"] + rows + [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]


def s11():
    rk = {"bbbc041_all": "prop_minm/r_mean|bbbc041_all", "mpidb_full": "mpidb_hc_b2__mh_s7_results/analytic_fields|a2|mh|ni_over_nmax_mean", "lacuna_B": "lacuna_hc_paired_B__lh_s7_results/analytic_fields|a2|lh|ni_over_nmax_mean"}
    r = [f"{n} & {P(rk[s],'f3')} & {P(f'prop_minm/m_star|delta|{s}','d')} & {P(f'prop_minm/m_star|delta_over_K|{s}','d')} \\\\" for s, n in SITES]
    T = [r"\begin{tabular}{@{}lrrr@{}}", r"\toprule", r"Site & Mean size ratio $\bar r$ & Minimum $m^\ast$ at level $\delta$ & Minimum $m^\ast$ at level $\delta/K$ \\", r"\midrule"] + r + [r"\midrule",
         f"Size-free $n_{{\\min}}$ at level $\\delta$ / $\\delta/K$ & -- & {P('s5_results/analytic|0.1|n_min','d')} & {P('prop_minm/n_min|delta_over_K','d')} \\\\", f"Annotated fields $m$ & -- & {P('protocol/m','d')} & {P('protocol/m','d')} \\\\", r"\bottomrule", r"\end{tabular}"]
    return head("tab:s11_mstar", None) + T + [r"\end{table}", ""]


def main():
    ap = argparse.ArgumentParser(); [ap.add_argument("--" + k, required=True) for k in ("master", "spec", "captions", "s6caption", "s4dcaption", "out")]; a = ap.parse_args(); out = Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    if hashlib.sha256(open(a.master, "rb").read()).hexdigest() != MASTER_SHA or not hashlib.sha256(open(a.spec, "rb").read()).hexdigest().startswith(SPEC_SHA): print("GATE (exit 3): pin differs", file=sys.stderr); return 3
    v = json.load(open(a.master))["values"]
    cap = open(a.captions, encoding="utf-8").read()
    if hashlib.sha256(cap.encode("utf-8")).hexdigest() != CAP_SHA: print("GATE (exit 3): captions file sha differs", file=sys.stderr); return 3
    caps = {}
    for n, (sn, sec) in CAPSEC.items():
        parts = cap.split(f"%%% {sec}\n")
        if len(parts) != 2: print("GATE (exit 3): caption section", sec, file=sys.stderr); return 3
        c = parts[1].split("\n%%% ")[0].strip()
        if not c or "\n\n" in c: print("GATE (exit 3): caption not one paragraph", sec, file=sys.stderr); return 3
        caps[n] = c
    # v1.6: the caption file text is s6txt (v1.5 named it s6 and shadowed the function s6(); the run crashed with TypeError before writing)
    s6txt = open(a.s6caption, encoding="utf-8").read()
    if hashlib.sha256(s6txt.encode("utf-8")).hexdigest() != S6CAP_SHA: print("GATE (exit 3): S6 caption file sha differs", file=sys.stderr); return 3
    s6p = [l for l in s6txt.split("\n## 2 ")[0].split("\n") if l.startswith("Unequal-probability field sampling on the development sites")]
    if len(s6p) != 1: print("GATE (exit 3): S6 caption paragraph", len(s6p), file=sys.stderr); return 3
    caps["tab_s6_unequal.tex"] = s6p[0].strip()
    # v1.5: panel D caption = the one paragraph of section '## 2' of the pinned delivery file 0b1357e3 that starts with 'Positives-only certificate at the census threshold'
    s4df = open(a.s4dcaption, encoding="utf-8").read()
    if hashlib.sha256(s4df.encode("utf-8")).hexdigest() != S4DCAP_SHA or s4df.count("\n## 2 ") != 1: print("GATE (exit 3): S4 panel D caption file", file=sys.stderr); return 3
    paras = [x.strip() for x in s4df.split("\n## 2 ")[1].split("\n", 1)[1].split("\n## ")[0].split("\n\n") if x.strip().startswith("Positives-only certificate at the census threshold")]
    if len(paras) != 1 or "\n" in paras[0]: print("GATE (exit 3): S4 panel D caption paragraph", len(paras), file=sys.stderr); return 3
    s4dc = paras[0]
    m4 = sorted({k.split("|")[3] for k in v if k.startswith("m51/result|M4_report_only|cells|")})
    if len(m4) != 8 or any(not re.match(r"^main_pi[0-9.]+_xi(lo|mid|hi)_cv[0-9.]+_eps[0-9.]+_m22$", c) for c in m4): print("GATE (exit 3): M4 cells", m4, file=sys.stderr); return 3
    if any(k in v for k in ABSENT): print("GATE (exit 3): a key expected absent is present", file=sys.stderr); return 3
    try: t4d = s4d(v)
    except ValueError as e: print("GATE (exit 3): S4 panel D", e, file=sys.stderr); return 3
    tabs = {"tab_s4_zplus.tex": s4(), "tab_s5_variant.tex": s5(), "tab_s6_unequal.tex": s6(m4), "tab_s11_mstar.tex": s11(), "tab_s4_zplus_d.tex": t4d}
    bad = sorted({k for t in tabs.values() for k in re.findall(r"\[\[\s*([^\[\]]+?)\s*::\s*\w+\s*\]\]", "\n".join(t)) if k not in v and k not in PENDING})   # v1.2: [^[]] as in make_tab_sites v1.3
    small = []   # v1.2: a value with 0 < |x| < 0.0005 printed with f3 would show as 0.000 / -0.000; such placeholders get f4 (the rule of make_tab_m51split v1.1)
    def f4fix(m):
        k = m.group(1)
        x = v.get(k)
        if isinstance(x, (int, float)) and not isinstance(x, bool) and 0 < abs(x) < 0.0005: small.append(k); return f"[[{k} :: f4]]"
        return m.group(0)
    ph = r"\caption{[[MISSING: caption of SI table S4D (writing session)]]}"   # v1.4
    if sum(l == ph for l in tabs["tab_s4_zplus_d.tex"]) != 1: print("GATE (exit 3): S4D caption placeholder", file=sys.stderr); return 3
    tabs["tab_s4_zplus_d.tex"] = [(r"\caption{" + s4dc.strip() + "}") if l == ph else l for l in tabs["tab_s4_zplus_d.tex"]]
    for n, (sn, sec) in CAPSEC.items():   # v1.2: the placeholder caption line is replaced by the writing session's caption (exactly one placeholder per table)
        ph = r"\caption{[[MISSING: caption of SI table " + sn + r" (writing session)]]}"
        if sum(l == ph for l in tabs[n]) != 1: print("GATE (exit 3): caption placeholder", n, file=sys.stderr); return 3
        tabs[n] = [(r"\caption{" + caps[n] + "}") if l == ph else l for l in tabs[n]]
    tabs = {n: [re.sub(r"\[\[([^\[\]]+?) :: f3\]\]", f4fix, l) for l in tl] for n, tl in tabs.items()}   # Gemini minor 1: after the captions are in
    # v1.2 (Gemini B2): PENDING keys allowed here too, as in the first gate
    bad2 = sorted({k for t in tabs.values() for k in re.findall(r"\[\[\s*([^\[\]]+?)\s*::\s*\w+\s*\]\]", "\n".join(t)) if k not in v and k not in PENDING})
    if bad2 or any("MISSING" in l for tt in tabs.values() for l in tt): print("GATE (exit 3): caption keys not in master or MISSING left", bad2, file=sys.stderr); return 3
    if bad: print("GATE (exit 3): placeholder keys not in master:", bad, file=sys.stderr); return 3
    out.mkdir(parents=True)
    for n, t in tabs.items():
        with open(out / n, "x", encoding="utf-8") as fh: fh.write("\n".join(t))
    print("f3 -> f4 (|x| < 0.0005):", small); print({n: hashlib.sha256(open(out / n, "rb").read()).hexdigest() for n in tabs}, {n: len(re.findall(r"\[\[", "\n".join(t))) for n, t in tabs.items()}); return 0


if __name__ == "__main__": sys.exit(main())

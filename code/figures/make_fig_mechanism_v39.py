# P5a v3.9 Figure 2 (fig:mechanism): how rule S3 works on one real draw.  Base: make_fig_mechanism_v2.py (89528999, fig0_mechanism); design r2 (ca3bf65e) Figure 2 section; width 139 mm, text >= 7 pt.
# Data, read only, every number shown comes from these two sha-pinned export files (figmech_v39_export.py v1.1 89e909a5, rc 0) or from results_master protocol keys:
#   figmech_v39_numbers.json a4b23923 and figmech_v39_arrays.npz 3ee12fc9 (in <tree>/results_src/figmech_v39/ after add_figmech_v39_keys.py, or the pulled export dir);
#   results_master protocol/eps_main, delta, candidates_K, m, mechanism_draw_seed (gated equal to the export json).  No number is typed into the figure.
# Changes from the base (r2 list): numbers from keys; panel (a) uses the THRFIX threshold; (b) adds the per-field binomial 5-95 % range under the pooled p_hat (cells exchangeable) and DE = n / n_eff;
#   (e) adds the Clopper-Pearson bound that treats the labelled negatives as independent cells (chosen candidate only, vermillion dashed); no shading because the export flag
#   cp_ok_true_over_any is false (the script draws shading only if it is true); candidates in neutral colours (c0 grey, c1s black, c2s sky blue, each with its own line style and marker),
#   the chosen one thicker; S3 orange and SPLIT6 blue are not used here; sentence-case panel titles; Okabe-Ito colours on marks only, all text pure black.
# Output (new directory, 'x' mode): fig2_mechanism_v39.pdf, fig2_mechanism_v39.png (300 dpi), fig2_mechanism_v39_numbers.json (every number printed in the figure, with its source key / path), SHA256SUMS.
# Usage: python3 -B make_fig_mechanism_v39.py --data <dir with the two export files> --master <results_master.json> --out <new dir>
# public-repo patch (2026-11-28): the master is read through public_inputs.load_master, which accepts the public redacted results/results_master.json.gz
#   (uncompressed sha 4b296884) or the pinned internal master; non-public inputs stop with a plain message.  Drawing code unchanged.
import argparse, hashlib, io, json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from scipy.stats import gaussian_kde
# v2.0 (user: IEEE style; module p5a_style_ieee v2.0, imported as P because this script uses `st` as a variable): Times, IEEE panel labels (a)-(e), constrained layout.
#   HERO = the chosen candidate (S3 orange: line in (e), bar in (d)); other candidates and baselines MUTED; FPR above / at or below eps as light bands in (e);
#   BBBC-ALL mass above tau_src as a light vermillion area in (a).  At most one number per panel: (a) BBBC-ALL FPR at tau_src, (b) DE, (d) chosen sensitivity, (c) and (e) none.
#   Inputs, sha pins and gates unchanged; shown keeps only the printed numbers (plus c_shift_ticks, chosen, m).
import p5a_style_ieee as P
import public_inputs as PI
NUM = ("figmech_v39_numbers.json", "a4b23923a65d4284b9580529b70a1845565aa0e5460177e55eb8b85c5442945c")
ARR = ("figmech_v39_arrays.npz", "3ee12fc9a32e285d6c181e2d03a19286ade9ca789b2c2efc2eaa55ea1f27c9d5")
NUM_PUBLIC = "8190467ee1d98929eb45dfc032bdf5d27761a060b191753ec1d11104ff2c9f75"
W, H = P.figsize(150)
BLACK, VERM, SKY, S3C = P.BLACK, P.VERM, P.SKY, P.S3C
G1, G2, G3 = P.MUTED[6], P.MUTED[0], P.MUTED[5]
CEN = ("c0", "c1s", "c2s"); LAB = {"c0": "$c^{0}$", "c1s": "$c^{1s}$", "c2s": "$c^{2s}$"}
STY = {"c0": dict(color=G2, ls=(0, (4, 2)), marker="s"), "c1s": dict(color=BLACK, ls="-", marker="o"), "c2s": dict(color=SKY, ls=(0, (5, 1.5, 1, 1.5)), marker="^")}


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); ap.add_argument("--master", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    D, out = Path(a.data), Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    if sha(D / NUM[0]) not in (NUM[1], NUM_PUBLIC): return fail("input sha differs", NUM[0])
    PI.need_file(D / ARR[0], ARR[1], "figmech_v39_arrays.npz (per-cell scores of the drawn BBBC-ALL fields)")
    J = json.load(open(D / NUM[0])); A = np.load(D / ARR[0]); M, MSHA = PI.load_master(a.master); V = M["values"]
    eps = V["protocol/eps_main"]; delta = V["protocol/delta"]; K = V["protocol/candidates_K"]; m = V["protocol/m"]; seed = V["protocol/mechanism_draw_seed"]
    if not (delta == J["delta"] and K == J["K"] and m == J["m"] and seed == J["seed"] and str(eps) in J["by_eps"]): return fail("protocol keys differ from the export json")
    E = J["by_eps"][str(eps)]; ch = E["chosen"]
    if not E["certified"]: return fail("main-eps draw not certified; the figure design assumes a certified draw")
    js = E["candidates"][ch]["last_index"]; ni = A["n_f"]
    shown = {}

    def num(name, val, src, fmt):
        shown[name] = dict(value=val, source=src, shown=format(val, fmt)); return format(val, fmt)
    P.apply()
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.2], width_ratios=[1.15, 1])
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    gsC = gs[1, 0].subgridspec(1, 2); axC1 = fig.add_subplot(gsC[0, 0]); axC2 = fig.add_subplot(gsC[0, 1], sharey=axC1)
    axD = fig.add_subplot(gs[1, 1]); axE = fig.add_subplot(gs[2, :])
    # (a) source threshold at a new site; the BBBC-ALL mass above tau_src (its FPR) as a light vermillion area; one number: that FPR
    thr = J["tau_src_logit"]; xs = np.linspace(-6, 6, 700)
    for dat, lab, c, lw, fc, al in ((A["nih_val_neg_s"], "NIH validation", G2, 0.9, G2, 0.35), (A["bbbc_all_neg_s"], "BBBC-ALL", BLACK, 1.4, VERM, 0.30)):
        sub = dat if len(dat) <= 20000 else dat[np.random.default_rng(0).choice(len(dat), 20000, replace=False)]
        d = gaussian_kde(sub)(xs); axA.plot(xs, d, color=c, lw=lw, label=lab); axA.fill_between(xs, 0, d, where=xs >= thr, color=fc, alpha=al, lw=0)
    P.ref_line(axA, thr, axis="x"); top = axA.get_ylim()[1] * 1.2; axA.set_ylim(0, top)
    axA.text(thr - 0.15, top * 0.97, r"$\tau_{\mathrm{src}}$", ha="right", va="top", fontsize=P.FS_LEGEND)
    f2 = num("fpr_bbbc_all", J["fpr_at_tau_src_bbbc_all"], "figmech_v39/fpr_at_tau_src_bbbc_all", ".3f")
    axA.annotate(f"FPR {f2}", xy=(1.2, top * 0.08), xytext=(3.2, top * 0.42), fontsize=P.FS_LEGEND, ha="center", arrowprops=dict(arrowstyle="-", lw=0.5, color=BLACK))
    axA.set_xlabel("Score $s$ of negative cells"); axA.set_ylabel("Density"); axA.set_yticks([]); axA.set_xlim(-6, 6)
    axA.legend(loc="upper right", handlelength=1.2)
    P.panel_label(axA, "a", "Source threshold at a new site")
    # (b) fields are not independent; one number: DE
    pf = np.array(E["field_frac"]); q05 = np.array(E["binom_q05_frac"]); q95 = np.array(E["binom_q95_frac"]); order = np.argsort(-pf, kind="stable"); xb = np.arange(m)
    for i, o in enumerate(order): axB.add_patch(Rectangle((i - 0.32, q05[o]), 0.64, q95[o] - q05[o], facecolor=P.MUTED[5], edgecolor="none", zorder=1))
    axB.scatter(xb, pf[order], s=ni[order] / ni.max() * 26 + 5, color=BLACK, edgecolor="white", linewidth=0.4, zorder=3)
    P.ref_line(axB, E["phat_pooled"], axis="y"); axB.set_xticks([]); axB.set_xlim(-1, m)
    ytop = max(pf.max(), q95.max()) * 2.4; axB.set_ylim(0, ytop)
    de = num("DE", E["n_labelled_negatives"] / E["n_eff"], f"figmech_v39/by_eps|{eps}|n_labelled_negatives / n_eff", ".2f")
    axB.text(m - 0.5, ytop * 0.96, f"DE = {de}", ha="right", va="top", fontsize=P.FS_LEGEND)
    axB.set_xlabel("Fields, sorted (size $\\propto$ negatives)"); axB.set_ylabel(r"Share above $\tau^{*}$")
    axB.legend([Patch(facecolor=P.MUTED[5], edgecolor="none"), Line2D([], [], color=BLACK, lw=P.LW_REF, ls=P.DASH_REF)], ["5–95 % if cells exchangeable", r"pooled $\hat p$"],
               loc="upper left", handlelength=1.2, borderaxespad=0.1)
    P.panel_label(axB, "b", "Fields are not independent")
    # (c) centering without labels (no number labels in the panel; the shifts are the tick labels)
    idx = A["field_neg_index"]; six = [int(order[i]) for i in np.linspace(0, m - 1, 6).astype(int)]; rj = np.random.default_rng(1)
    u = A["c0_field_neg_centered"]; cc = A[f"{ch}_field_neg_centered"]; allv = np.concatenate([u[idx == f] for f in six] + [cc[idx == f] for f in six]); lo, hi = np.percentile(allv, [0.5, 99.5])
    for q, f in enumerate(six):
        for ax, v, col in ((axC1, u[idx == f], G2), (axC2, cc[idx == f], G1)):
            ax.scatter(v, np.full(len(v), q) + rj.uniform(-.27, .27, len(v)), s=1.4, color=col, linewidths=0, rasterized=True)
            ax.plot([np.median(v)] * 2, [q - .38, q + .38], color=BLACK, lw=1.2)
    for ax in (axC1, axC2): ax.set_xlim(lo - 0.3, hi + 0.3); ax.set_ylim(-0.7, 5.7)
    shv = A[f"{ch}_shift"]; axC1.set_yticks(range(6)); axC1.set_yticklabels([format(shv[f], "+.2f").replace("-", "−") for f in six])
    shown["c_shift_ticks"] = dict(value=[float(shv[f]) for f in six], source=f"figmech_v39_arrays.npz[{ch}_shift]", shown=[format(shv[f], "+.2f") for f in six])
    axC1.set_ylabel(f"Shift of {LAB[ch]}"); plt.setp(axC2.get_yticklabels(), visible=False); axC1.set_xlabel("Before"); axC2.set_xlabel(f"After ({LAB[ch]})")
    P.panel_label(axC1, "c", "Centering without labels")
    # (d) selection: the chosen candidate in S3 orange (HERO), the others MUTED; one number: the chosen sensitivity
    for i, k in enumerate(CEN):
        c = E["candidates"][k]; v = c["sens_labelled"] if c["certified"] else 0.0
        axD.bar(i, v, width=0.6, color=S3C if k == ch else P.MUTED[2], edgecolor="none", zorder=3)
        if k == ch: axD.text(i, v + 0.03, num(f"sens_labelled_{k}", v, f"figmech_v39/by_eps|{eps}|candidates|{k}|sens_labelled", ".3f"), ha="center", va="bottom", fontsize=P.FS_LEGEND)
        elif not c["certified"]: axD.text(i, 0.03, "not\ncertified", ha="center", va="bottom", fontsize=P.FS_LEGEND)
    axD.set_xticks(range(3)); axD.set_xticklabels([LAB[k] for k in CEN]); axD.set_ylim(0, 1.12); axD.set_yticks([0, 0.5, 1.0]); axD.set_ylabel("Sensitivity (labeled)")
    P.ygrid(axD); axD.set_xlabel(f"Selected: {LAB[ch]}"); P.panel_label(axD, "d", "Selection")
    # (e) certification scan: the chosen design-effect bound is the HERO line (S3 orange); FPR above eps = light vermillion band, at or below = light blue band; no number
    xmax = 0.0; ymax = 0.19
    P.band(axE, eps, ymax, color=VERM, alpha=0.08); P.band(axE, 0, eps, color=SKY, alpha=0.14)
    for k in CEN:
        c = E["candidates"][k]; st = c["stop_index"]; hi_ = st if st is not None else len(A[f"{k}_grid"]) - 1; fr = A[f"{k}_frac"][:hi_ + 1]; ub = A[f"{k}_ub"][:hi_ + 1]; xmax = max(xmax, fr[-1])
        col = S3C if k == ch else STY[k]["color"]
        if k == ch: axE.plot(fr, ub, color=col, lw=P.HERO_KW["linewidth"], zorder=6, label=f"design-effect bound, {LAB[k]} (chosen)")
        else: axE.plot(fr, ub, color=col, lw=P.MUTED_KW["linewidth"], ls=STY[k]["ls"], zorder=3, label=f"design-effect bound, {LAB[k]}")
        if k == ch: axE.plot(fr, A[f"{k}_true_fpr"][:hi_ + 1], color=BLACK, lw=0.8, ls=(0, (1, 1.5)), zorder=2, label=f"true site FPR, {LAB[k]} (oracle, not used)")
        mkw = dict(ls="", marker=STY[k]["marker"], ms=P.HERO_KW["markersize"] if k == ch else 3.6, zorder=7)
        if c["certified"]: axE.plot([fr[c["last_index"]]], [ub[c["last_index"]]], color=col, mec="white" if k == ch else col, mew=0.5, **mkw)
        if st is not None: axE.plot([fr[st]], [ub[st]], mfc="white", mec=col, mew=0.9, **mkw)
    stc = E["candidates"][ch]["stop_index"]; hic = stc if stc is not None else len(A[f"{ch}_grid"]) - 1
    axE.plot(A[f"{ch}_frac"][:hic + 1], A[f"{ch}_cp"][:hic + 1], color=VERM, lw=1.0, ls=(0, (4, 2)), label=f"cells treated as independent, {LAB[ch]}", zorder=5)
    if E["cp_ok_true_over_any"]:
        msk = (A[f"{ch}_cp"][:hic + 1] <= eps) & (A[f"{ch}_true_fpr"][:hic + 1] > eps)
        axE.fill_between(A[f"{ch}_frac"][:hic + 1], 0, 1, where=msk, color=VERM, alpha=0.12, lw=0, transform=axE.get_xaxis_transform())
    xmax *= 1.04; axE.set_xlim(0, xmax); axE.set_ylim(0, ymax)
    P.ref_line(axE, eps, axis="y"); axE.set_yticks([0, 0.05, 0.10, 0.15]); axE.text(0.995 * xmax, eps + 0.003, r"$\epsilon$", ha="right", va="bottom", fontsize=P.FS_TEXT)
    xs_ = A[f"{ch}_frac"][js]; axE.axvline(xs_, color=G2, lw=0.6, ls=(0, (1, 1.5)), zorder=1); axE.text(xs_, ymax * 0.97, r"$\tau^{*}$ ", ha="right", va="top", fontsize=P.FS_LEGEND)
    axE.set_xlabel(r"Share of labeled negatives at or above $\tau$ (threshold decreasing $\rightarrow$)"); axE.set_ylabel("FPR")
    h, l = axE.get_legend_handles_labels()
    h += [Line2D([], [], ls="", marker="o", ms=4.0, color=G2, mec=G2), Line2D([], [], ls="", marker="o", ms=4.0, mfc="white", mec=G2, mew=0.9)]
    l += [r"filled: certified $\tau_k$", "open: scan stops"]
    axE.legend(h, l, loc="upper left", ncol=2, handlelength=1.4, columnspacing=1.2)
    P.panel_label(axE, "e", r"Certification scan, each candidate at level $\delta/K$")
    shown["chosen"] = dict(value=ch, source=f"figmech_v39_numbers.json by_eps|{eps}|chosen"); shown["m"] = dict(value=m, source="protocol/m")
    buf = P.render(fig)
    plt.close(fig)
    payload = json.dumps(dict(script_sha256=sha(__file__), style_sha256=P.STYLE_SHA256, inputs={NUM[0]: sha(D / NUM[0]), ARR[0]: ARR[1]}, master=str(a.master), master_sha256=MSHA, eps_main=eps, shown=shown), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    for ext in ("pdf", "png"):
        with open(out / f"fig2_mechanism_v39.{ext}", "xb") as fh: fh.write(buf[ext])
    with open(out / "fig2_mechanism_v39_numbers.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh:
        for nm in ("fig2_mechanism_v39.pdf", "fig2_mechanism_v39.png", "fig2_mechanism_v39_numbers.json"): fh.write(f"{sha(out / nm)}  {nm}\n")
    print("written", out); return 0


if __name__ == "__main__": sys.exit(main())

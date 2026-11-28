# P5a v3.9 Figure 5 (fig:dev): development sites as the number of annotated fields m varies.  Bases: make_figures_v2.py (5c066c31) Fig. 4 block (fig4_coverage_sensitivity_vs_m) and
#   handover_20261101/fig_dense_v4_20261115T0040Z/make_dense_v4.py (fig4_main_grid_v4, cannot-certify column and per-seed dots / band).  Design r2 (ca3bf65e) Figure 5 section; width 139 mm.
# Data: only results_master keys ebatch_e2p/cells|<site>|<m>|<eps_main>|{S3, SPLIT6, MC, CP}|{coverage, coverage_per_seed|i, coverage_ci90_seedcluster|0/1, ds_own, fpr_mean, cannot_certify,
#   vacuous_share}; ebatch_e2p/ms|i; ebatch_e2p/minm_boundary|prop_minm/m_star|delta_over_K|<site>; protocol/{delta, eps_main}.  MC = marginal split-conformal (cells as independent),
#   CP = exact rule on pooled negative cells; neither name appears as an abbreviation in the figure.
# Layout (r2): 4 rows x 3 columns (BBBC-ALL, MP-IDB, Lacuna-B); rows (a-c) coverage with a broken axis (upper: S3 and SPLIT6 with per-seed dots and the 90 % seed-cluster band;
#   lower: the two baselines), (d-f) Delta Sens_own, (g-i) mean realized FPR (with the baselines), (j-l) cannot-certify rate; log m axis with ticks at the grid values.
#   Open markers: S3 / SPLIT6 values whose draws include fallback draws (cannot_certify > 0), baselines whose draws include vacuous draws (vacuous_share > 0) (rows 1-3).
#   1 - delta and the target eps as dashed lines (values from keys); m* at level delta/K written as text in the cannot-certify panels (no vertical line; m* lies beyond the grid).
# Gates (exit 3): every m* exceeds the largest m (the caption's last sentence); every cell present for the 7 m values x 4 rules; no cell flagged incomplete.
# Output (new directory, 'x' mode): fig5_dev_v39.pdf, .png (300 dpi), fig5_dev_v39_numbers.json, SHA256SUMS.
# Usage: python3 -B make_fig_dev_v39.py --master <results_master.json> --out <new dir>
# public-repo patch (2026-11-28): the master is read through public_inputs.load_master, which accepts the public redacted results/results_master.json.gz
#   (uncompressed sha 4b296884) or the pinned internal master; non-public inputs stop with a plain message.  Drawing code unchanged.
import argparse, hashlib, io, json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
# v2.0 (user: IEEE style; module p5a_style_ieee v2.0): Times, IEEE labels (a)-(l); S3 and SPLIT6 HERO (lw 1.7, white-edged markers), the two baselines MUTED (lw 0.75);
#   regimes as light bands (coverage >= 1 - delta blue, the lower baseline axis and FPR > eps vermillion); legend without numbers (1 - delta, eps, seed count go to the
#   caption); one number per panel at most (m* in the cannot-certify panels).  The hand-placed grid (broken axes, label coordinates) is kept, so this figure is saved
#   with constrained_layout off.  Data, keys, gates and layout unchanged.
import p5a_style_ieee as P
import public_inputs as PI
# v2.1 (实验管理): height 166 -> 148 mm so that the figure and its caption fit the BJ text height (v2.0 gave 'Float too large' by 41.7 pt); legend band (24.6 mm) and
#   bottom margin (10.6 mm) kept in mm, so only the panel rows shrink; data, text sizes and layout otherwise unchanged.
# v2.2 (user 2026-11-24, via 实验管理): rule names in the figure text renamed: S3 -> DE-Grid, SPLIT6 -> Split-Select, Z+ -> Pos-Census, zero-label rule Z -> Zero-Census; keys and logic unchanged.
W, H = P.figsize(148)
BLACK, S3C, SPC, SKY, VERM = P.BLACK, P.S3C, P.SPC, P.SKY, P.VERM
G1, G2, G3 = P.MUTED[6], P.MUTED[0], P.MUTED[5]
SITES = [("bbbc041_all", "BBBC-ALL"), ("mpidb_full", "MP-IDB"), ("lacuna_B", "Lacuna-B")]
RULES = {"S3": dict(color=S3C, marker="o", ls="-"), "SPLIT6": dict(color=SPC, marker="s", ls="-"),
         "MC": dict(color=G2, marker="^", ls=(0, (3, 1.5))), "CP": dict(color=G1, marker="v", ls=(0, (1, 1.2)))}
LABEL = {"S3": "DE-Grid", "SPLIT6": "Split-Select", "MC": "Marginal split-conformal (cells as independent)", "CP": "Exact rule on pooled negative cells"}


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--master", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(); out = Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    M, MSHA = PI.load_master(a.master); V = M["values"]; delta = V["protocol/delta"]; eps = V["protocol/eps_main"]; E = format(eps, "g")
    ms = [V[f"ebatch_e2p/ms|{i}"] for i in range(len([k for k in V if k.startswith("ebatch_e2p/ms|")]))]
    mstar = {s: V[f"ebatch_e2p/minm_boundary|prop_minm/m_star|delta_over_K|{s}"] for s, _ in SITES}
    if not all(mstar[s] > max(ms) for s in mstar): return fail("an m* does not exceed the largest m", mstar, max(ms))

    def g(s, m, r, f, default=None): return V.get(f"ebatch_e2p/cells|{s}|{m}|{E}|{r}|{f}", default)
    for s, _ in SITES:
        for m in ms:
            if V.get(f"ebatch_e2p/cells|{s}|{m}|{E}|incomplete") not in (False, None): return fail("incomplete cell", s, m)
            for r in RULES:
                if g(s, m, r, "coverage") is None or g(s, m, r, "fpr_mean") is None: return fail("missing cell", s, m, r)
    nseed = len([k for k in V if k.startswith(f"ebatch_e2p/cells|{SITES[0][0]}|{ms[0]}|{E}|S3|coverage_per_seed|")])
    P.apply()
    fig = plt.figure(figsize=(W, H), constrained_layout=False)
    gs = fig.add_gridspec(4, 3, left=0.105, right=0.985, top=0.834, bottom=0.072, hspace=0.42, wspace=0.32, height_ratios=[1.35, 1, 1, 1])
    lo_lim = (min(min(g(s, m, r, "coverage") for m in ms) for s, _ in SITES for r in ("MC", "CP")) - 0.05, max(max(g(s, m, r, "coverage") for m in ms) for s, _ in SITES for r in ("MC", "CP")) + 0.05)
    hi_lim = (min(1 - delta, min(min(g(s, m, r, "coverage_ci90_seedcluster|0") for m in ms) for s, _ in SITES for r in ("S3", "SPLIT6"))) - 0.015, 1.006)
    letters = iter("abcdefghijkl"); shown = {"m_star": mstar, "eps": eps, "delta": delta}
    xt = np.array(ms, float)

    def mk(ax, s, r, field, open_if):
        y = [g(s, m, r, field) for m in ms]; st = RULES[r]; hero = r in ("S3", "SPLIT6")
        ax.plot(xt, y, color=st["color"], ls=st["ls"], lw=P.HERO_KW["linewidth"] if hero else P.MUTED_KW["linewidth"], zorder=5 if hero else 3)
        for x_, y_, m in zip(xt, y, ms):
            op = open_if(m); ms_ = 4.0 if hero else 2.8
            ax.plot([x_], [y_], st["marker"], ms=ms_, color=st["color"], mfc="white" if op else st["color"], mec=st["color"] if op else ("white" if hero else st["color"]),
                    mew=0.9 if op else 0.5, zorder=6 if hero else 4)

    fb = lambda s, r: (lambda m: (g(s, m, r, "cannot_certify", 0) or 0) > 0)
    vac = lambda s, r: (lambda m: (g(s, m, r, "vacuous_share", 0) or 0) > 0)
    axes_rows = [[], [], [], []]
    for j, (s, name) in enumerate(SITES):
        sub = gs[0, j].subgridspec(2, 1, height_ratios=[2.1, 1], hspace=0.08); au = fig.add_subplot(sub[0]); al = fig.add_subplot(sub[1], sharex=au)
        for r in ("S3", "SPLIT6"):
            lo = [g(s, m, r, "coverage_ci90_seedcluster|0") for m in ms]; hi = [g(s, m, r, "coverage_ci90_seedcluster|1") for m in ms]
            au.fill_between(xt, lo, hi, color=RULES[r]["color"], alpha=0.18, lw=0, zorder=1)
            for i in range(nseed):
                au.scatter(xt * (1 + (0.025 if r == "SPLIT6" else -0.025)), [g(s, m, r, f"coverage_per_seed|{i}") for m in ms], s=2.5, color=RULES[r]["color"], alpha=0.6, lw=0, zorder=2)
            mk(au, s, r, "coverage", fb(s, r))
        for r in ("MC", "CP"): mk(al, s, r, "coverage", vac(s, r))
        P.ref_line(au, 1 - delta, axis="y"); P.band(au, 1 - delta, hi_lim[1], color=SKY, alpha=0.16); P.band(au, hi_lim[0], 1 - delta, color=VERM, alpha=0.08)
        P.band(al, lo_lim[0], lo_lim[1], color=VERM, alpha=0.08)
        au.set_ylim(*hi_lim); al.set_ylim(*lo_lim); au.spines["bottom"].set_visible(False); al.spines["top"].set_visible(False)
        au.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        d = 0.012
        for ax_, yy in ((au, 0), (al, 1)): ax_.plot([-d * 2, d * 2], [yy - d * 3, yy + d * 3], transform=ax_.transAxes, color=BLACK, lw=0.6, clip_on=False)
        au.text(0.5, 1.05, name, transform=au.transAxes, ha="center", va="bottom", fontsize=P.FS_TEXT)
        P.panel_label(au, next(letters), x=-0.14, y=1.02)
        if j == 0: au.set_ylabel("Coverage", labelpad=2); au.yaxis.set_label_coords(-0.27, 0.2)
        axes_rows[0].append((au, al))
    for row, (field, ylab) in enumerate((("ds_own", r"$\Delta\mathrm{Sens}_{\mathrm{own}}$"), ("fpr_mean", "Mean realized FPR"), ("cannot_certify", "Cannot-certify rate")), start=1):
        for j, (s, name) in enumerate(SITES):
            ax = fig.add_subplot(gs[row, j])
            if row == 1: P.ref_line(ax, 0, axis="y", color=G2, dashed=False, lw=P.LW_AXIS)
            if row == 2: P.ref_line(ax, eps, axis="y")
            rules = ("S3", "SPLIT6", "MC", "CP") if row == 2 else ("S3", "SPLIT6")
            for r in rules:
                op = (vac(s, r) if r in ("MC", "CP") else fb(s, r)) if row < 3 else (lambda m: False)
                mk(ax, s, r, field, op)
            if row == 2:
                # bands after the curves, to the current limits (a band to 1 would rescale the axis)
                yl = ax.get_ylim(); P.band(ax, eps, yl[1], color=VERM, alpha=0.08); P.band(ax, yl[0], eps, color=SKY, alpha=0.16); ax.set_ylim(yl)
            if row == 3:
                ax.text(0.97, 0.95, f"$m^{{*}}$ = {mstar[s]} $\\rightarrow$", transform=ax.transAxes, ha="right", va="top", fontsize=P.FS_LEGEND)
                ax.set_ylim(-0.03, max(0.12, max(max(g(s, m, r, field) for m in ms) for r in rules) * 1.25))
            P.panel_label(ax, next(letters), x=-0.14, y=1.02)
            if j == 0: ax.set_ylabel(ylab, labelpad=2); ax.yaxis.set_label_coords(-0.27, 0.5)
            axes_rows[row].append(ax)
    for row in axes_rows:
        for it in row:
            for ax in (it if isinstance(it, tuple) else (it,)):
                ax.set_xscale("log"); ax.set_xticks(xt); ax.set_xticklabels([str(m) for m in ms]); ax.minorticks_off(); ax.set_xlim(xt[0] / 1.18, xt[-1] * 1.18)
    for ax in axes_rows[3]: ax.set_xlabel("Annotated fields $m$")
    for row in axes_rows[1:3]:
        for ax in row: ax.tick_params(axis="x", labelbottom=True)
    for _, al in axes_rows[0]: al.tick_params(axis="x", labelbottom=True)
    hs = [Line2D([], [], color=RULES[r]["color"], ls=RULES[r]["ls"], lw=1.3 if r in ("S3", "SPLIT6") else P.MUTED_KW["linewidth"], marker=RULES[r]["marker"],
                 ms=4.0 if r in ("S3", "SPLIT6") else 2.8, mec="white" if r in ("S3", "SPLIT6") else RULES[r]["color"], mew=0.5) for r in RULES]
    hs += [Line2D([], [], ls="", marker="o", ms=4.0, mfc="white", mec=G1, mew=0.9), Line2D([], [], color=BLACK, lw=P.LW_REF, ls=P.DASH_REF), Patch(facecolor=S3C, alpha=0.18, edgecolor="none")]
    ls_ = [LABEL[r] for r in RULES] + ["includes fallback or vacuous draws", r"$1-\delta$ (coverage), $\epsilon$ (FPR)", "90 % band over training seeds"]
    fig.legend(hs, ls_, loc="upper center", ncol=2, bbox_to_anchor=(0.53, 0.997), columnspacing=1.2, handlelength=1.4, handletextpad=0.5)
    with matplotlib.rc_context({"figure.constrained_layout.use": False}): buf = P.render(fig)
    plt.close(fig)
    table = {s: {str(m): {r: {f: g(s, m, r, f) for f in ("coverage", "ds_own", "fpr_mean", "cannot_certify", "vacuous_share")} for r in RULES} for m in ms} for s, _ in SITES}
    payload = json.dumps(dict(script_sha256=sha(__file__), style_sha256=P.STYLE_SHA256, master=str(a.master), master_sha256=MSHA, ms=ms, shown=shown, cells=table), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    for ext in ("pdf", "png"):
        with open(out / f"fig5_dev_v39.{ext}", "xb") as fh: fh.write(buf[ext])
    with open(out / "fig5_dev_v39_numbers.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh:
        for nm in ("fig5_dev_v39.pdf", "fig5_dev_v39.png", "fig5_dev_v39_numbers.json"): fh.write(f"{sha(out / nm)}  {nm}\n")
    print("written", out); return 0


if __name__ == "__main__": sys.exit(main())

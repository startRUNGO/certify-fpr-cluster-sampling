# P5a figure F-c (new, 实验管理 / user 2026-11-24 figure set, spec approved 实验管理 2026-11-25): realized site FPR among certified draws in four simulation settings,
#   pooled over 50 synthetic sites of 200 draws.  WRITTEN, NOT RUN, until 实验管理 reviews it.  Full BJ width 139 mm (p5a_style_ieee v2.0, pinned), pure black text,
#   colour only on marks.
# Settings (fixed by design before any number was read, as in d_c_simbox.py 17223230): e10|0 base; e10|3 between-field shift xi = 1; e10|9 heavy tail t5;
#   e12|19 two-stage slides S = 12, omega = 0.8, s = 2 (the stronger slide stress by design).
# Data: only results_master keys (master dba06e5c; simbox/* from D-c, ingested by add_simbox_keys.py c19380c2):
#   simbox/fc|<d>|<ci>|<rule>|<eps>|{n, n_cert, q05, q25, q50, q75, q95, mean, min, max[, strip_values|i]} for the boxes;
#   sim50/<d>|<ci>|<rule>|<eps>|all50|coverage|mean for the coverage row above the boxes (2 decimals);
#   design labels from ebatch_e10/cells|<i>|factors|* and e12/cells|19|cell|*.
# Layout: 2 x 2 panels, one per setting; in each panel one box per rule in Table 3 order (Split-Select, DE-Grid, Cluster bootstrap, GEE, Beta-binomial), x labels
#   on two lines, rotation 0; the y axis (realized site FPR) is shared by the four panels.  DE-Grid and Split-Select in their role colours (box edge and whiskers), the
#   other rules in MUTED greys; box face white.
# Coverage row above the panel frame, between the frame and the panel title; y ticks only up to the highest whisker / eps.
# Encoding: box q25-q75, black median line, whiskers q05 and q95, white dot with black edge at the mean, no outlier points; dashed black line at eps; the coverage of each
#   rule printed in a row above the boxes (label 'Coverage' at the left of the row in the left-hand panels).
#   A rule with n_cert < STRIP_N in a setting is drawn as a strip of its certified values with 'n = <n_cert>' above it instead of a box (decided by the key, not by eye).
# Caption note: box definition, pooled draws, coverage row; 'Every draw certified in these four settings.' only if n_cert == n for every shown entry, else
#   'Boxes show certified draws only.'
# Gates (exit 3, nothing written): master sha MASTER_SHA; style sha STYLE_SHA; --eps in {0.1, 0.05}; every key present and finite; n == DRAWS_TOTAL for every entry;
#   min <= q05 <= q25 <= q50 <= q75 <= q95 <= max when n_cert > 0; strip_values present iff n_cert < STRIP_N, with length n_cert; design labels match the settings;
#   the eps line and every whisker inside the y range; x labels, coverage numbers and panel titles do not overlap each other; every text inside the figure bbox;
#   PDF height <= MAX_H_MM.
# r2 (实验管理 2026-11-25): (1) shared y axis from YLO = floor(0.8 x min q05 of the shown entries, and min strip value) to the tick step Y_STEP, gated so that every q05,
#   q95 and the eps line lie inside; caption 'The y axis does not start at zero.'; (2) light VERM band from eps to the top of the axis (alpha 0.08, as the F-d / F-e
#   bands), caption 'The shaded region marks realized FPR above eps.'; (3) sharex, rule names only under the bottom row, height 112 -> 100 mm (gated <= MAX_H_MM);
#   (4) 'over <n_realisations> synthetic sites of <n / n_realisations> draws' built from sim50/...|all50|n_realisations and |n, equal for every shown entry.
# Output (new directory, 'x' mode): fig_simbox_v66_eps<eps>.pdf, .png (600 dpi), fig_simbox_v66_eps<eps>_numbers.json, SHA256SUMS.
# Usage: python3 -E -B make_fig_simbox_v66.py --master <tree/results_master.json> --eps 0.1 --out <new dir>
import argparse, hashlib, json, math, sys
from matplotlib.ticker import MultipleLocator
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import p5a_style_ieee as P

MASTER_SHA = "dba06e5cf5c9f73b432eb2717bc542b990bb473ed16370f869207c1e6b0e12d1"
STYLE_SHA = "8424218990b9bec4dcb8b3bf6c54d7b2da5671b2fc9444511ccbefa7edd39ce3"
RULES = [("SPLIT6", "Split-\nSelect"), ("S3", "DE-\nGrid"), ("CB_c1s", "Cluster\nbootstrap"), ("GEE_c1s", "GEE"), ("BB_c1s", "Beta-\nbinomial")]
SETTINGS = [("e10", 0), ("e10", 3), ("e10", 9), ("e12", 19)]
STRIP_N = 100
DRAWS_TOTAL = 10000
FIG_H_MM = 100
Y_STEP = 0.01
BAND_ALPHA = 0.08
MAX_H_MM = 150
Q = ("min", "q05", "q25", "q50", "q75", "q95", "max")
STEM = "fig_simbox_v66"


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def load(V, es, shown):
    """Reads and gates every key; returns (titles, B, cov) or an error string."""
    def g(k):
        if k not in V or V[k] is None: raise KeyError(k)
        if isinstance(V[k], float) and not math.isfinite(V[k]): raise KeyError(f"non-finite {k}")
        shown[k] = V[k]; return V[k]
    f = lambda i, x: g(f"ebatch_e10/cells|{i}|factors|{x}")
    base = {x: f(0, x) for x in ("tau", "r", "pi", "tail")}
    if any(f(3, x) != base[x] for x in ("r", "pi", "tail")) or f(3, "tau") != 1.0: return "e10|3 is not the cell with shift xi = 1 only"
    if any(f(9, x) != base[x] for x in ("tau", "r", "pi")) or f(9, "tail") != "t5": return "e10|9 is not the cell with tail t5 only"
    c = "e12/cells|19|cell|"
    if (g(c + "name"), g(c + "slides|0"), g(c + "slides|1"), g(c + "two_stage")) != ("C2_S12_om0.8_s2", 12, 0.8, 2): return "e12|19 is not the two-stage cell S = 12, omega = 0.8, s = 2"
    titles = ["Base", r"Shift $\xi$ = 1", r"Heavy tail $t_5$", r"Two-stage, $S$ = 12, $\omega$ = 0.8, $s$ = 2"]
    B, cov, design = {}, {}, set()
    for d, ci in SETTINGS:
        for r, _ in RULES:
            b = f"simbox/fc|{d}|{ci}|{r}|{es}|"; n, nc = g(b + "n"), g(b + "n_cert")
            if n != DRAWS_TOTAL: return f"n differs from {DRAWS_TOTAL}: {b}"
            e = dict(n=n, n_cert=nc)
            if nc > 0:
                e.update({q: g(b + q) for q in Q}); e["mean"] = g(b + "mean")
                if not all(e[a] <= e[z] for a, z in zip(Q, Q[1:])): return f"quantiles out of order: {b}"
            has_strip = (b + "strip_values|0") in V
            if has_strip != (nc < STRIP_N and nc > 0): return f"strip_values present iff 0 < n_cert < {STRIP_N} fails: {b}"
            if nc < STRIP_N: e["strip"] = [g(b + f"strip_values|{i}") for i in range(nc)]
            B[(d, ci, r)] = e
            cov[(d, ci, r)] = g(f"sim50/{d}|{ci}|{r}|{es}|all50|coverage|mean")
            # r2: sites and draws per site from the master, not hardcoded
            nr, nn = g(f"sim50/{d}|{ci}|{r}|{es}|all50|n_realisations"), g(f"sim50/{d}|{ci}|{r}|{es}|all50|n")
            if nn != n or nn % nr: return f"all50|n differs from simbox n or is not a multiple of n_realisations: {b}"
            design.add((nr, nn // nr))
    if len(design) != 1: return f"sites x draws differ between entries: {sorted(design)}"
    return titles, B, cov, design.pop()


def draw(eps, titles, B, cov):
    """Draws the figure; returns (fig, None) or (None, error string).  Writes nothing."""
    P.apply()
    fig, axs = plt.subplots(2, 2, figsize=P.figsize(FIG_H_MM), sharey=True, sharex=True)
    axs = axs.ravel(); xs = list(range(1, len(RULES) + 1))
    hi = max([eps] + [e.get("q95", 0) for e in B.values()] + [max(e.get("strip", [0]) or [0]) for e in B.values()])
    ytop = hi * 1.08
    # r2: shared lower limit below the lowest whisker / strip point, floored to the tick step
    lowest = min([e["q05"] for e in B.values() if "q05" in e] + [v for e in B.values() for v in e.get("strip", [])] or [0.0])
    ylo = math.floor(0.8 * lowest / Y_STEP) * Y_STEP
    covtexts, xlabs, ptitles = [], [], []
    for ax, (d, ci), title, letter in zip(axs, SETTINGS, titles, "abcd"):
        stats, pos, cols = [], [], []
        for x, (r, _) in zip(xs, RULES):
            e = B[(d, ci, r)]; col = P.ROLE[r]["color"] if r in P.ROLE else P.MUTED[0]
            if e["n_cert"] >= STRIP_N:
                stats.append(dict(med=e["q50"], q1=e["q25"], q3=e["q75"], whislo=e["q05"], whishi=e["q95"], mean=e["mean"], fliers=[], label=r)); pos.append(x); cols.append(col)
            else:
                ax.plot([x] * len(e["strip"]), e["strip"], ls="none", marker="o", ms=2.2, mfc=col, mec=col, zorder=4)
                ax.annotate(f"n = {e['n_cert']}", xy=(x, max(e["strip"] or [0])), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=P.FS_MIN, color=P.BLACK)
        if stats:
            bp = ax.bxp(stats, positions=pos, widths=0.55, showmeans=True, showfliers=False, patch_artist=True,
                        medianprops=dict(color=P.BLACK, lw=1.0), meanprops=dict(marker="o", ms=3.2, mfc=P.WHITE, mec=P.BLACK, mew=0.6))
            for i, c in enumerate(cols):
                bp["boxes"][i].set(facecolor=P.WHITE, edgecolor=c, lw=1.1)
                for part in ("whiskers", "caps"):
                    for ln in bp[part][2 * i:2 * i + 2]: ln.set(color=c, lw=1.0)
        P.ref_line(ax, eps, axis="y", color=P.BLACK)
        # r2: failure region above eps
        P.band(ax, eps, ytop, axis="y", color=P.VERM, alpha=BAND_ALPHA)
        ax.set_xlim(0.4, len(RULES) + 0.6); ax.set_ylim(ylo, ytop)
        # y ticks only up to hi, so no tick label sits at the top edge next to the coverage row
        ax.yaxis.set_major_locator(MultipleLocator(Y_STEP)); ax.set_yticks([v for v in ax.get_yticks() if ylo - 1e-12 <= v <= hi + 1e-12]); ax.set_ylim(ylo, ytop)
        ax.set_xticks(xs); ax.set_xticklabels([h for _, h in RULES], fontsize=P.FS_MIN, linespacing=0.95); ax.tick_params(axis="x", length=0)
        # r2: rule names only under the bottom row
        if letter in "ab": ax.tick_params(axis="x", labelbottom=False)
        for x, (r, _) in zip(xs, RULES):
            covtexts.append(ax.annotate(f"{cov[(d, ci, r)]:.2f}", xy=(x, 1.0), xycoords=("data", "axes fraction"), xytext=(0, 2), textcoords="offset points",
                                        ha="center", va="bottom", fontsize=P.FS_MIN, color=P.BLACK, annotation_clip=False))
        if letter in "ac":
            ax.set_ylabel("Realized site FPR")
            covtexts.append(ax.annotate("Coverage", xy=(0.4, 1.0), xycoords=("data", "axes fraction"), xytext=(-2, 2), textcoords="offset points",
                                        ha="right", va="bottom", fontsize=P.FS_MIN, color=P.BLACK, annotation_clip=False))
        n0 = len(ax.texts); P.panel_label(ax, letter, title, pad_pt=12)
        ptitles.extend(ax.texts[n0:])
        if letter in "cd": xlabs.append(ax.get_xticklabels())
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    for ax in axs:
        y0, y1 = ax.get_ylim()
        if not (y0 <= eps <= y1): return None, "eps line outside the y range"
    # r2: every whisker end and strip point inside the shared y range
    for e in B.values():
        for v in [e[q] for q in ("q05", "q95") if q in e] + e.get("strip", []):
            if not (ylo <= v <= ytop): return None, f"a whisker or strip point outside the y range: {v}"

    def overlap(ts, what):
        bbs = [(t.get_text(), t.get_window_extent(rend)) for t in ts if t.get_text()]
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                if bbs[i][1].overlaps(bbs[j][1]): return f"{what}: '{bbs[i][0]}' overlaps '{bbs[j][0]}'"
        return None
    for ls in xlabs:
        m = overlap(ls, "x labels")
        if m: return None, m
    m = overlap(covtexts + ptitles + [t for ax in axs for t in ax.get_yticklabels()], "coverage numbers / panel titles / y tick labels")
    if m: return None, m
    # coverage numbers stay above every box and whisker of their panel
    for ax in axs:
        top = max([ln.get_window_extent(rend).y1 for ln in ax.lines if ln.get_visible()] + [p.get_window_extent(rend).y1 for p in ax.patches])
        for t in covtexts:
            if t.axes is ax and t.get_text() != "Coverage" and t.get_window_extent(rend).y0 <= top: return None, f"coverage number '{t.get_text()}' touches a box or whisker"
    fb = fig.bbox
    for t in covtexts + ptitles + [t for ls in xlabs for t in ls] + [ax.yaxis.label for ax in axs] + [t for ax in axs for t in ax.get_yticklabels()]:
        if not t.get_text(): continue
        bb = t.get_window_extent(rend)
        if bb.x0 < fb.x0 or bb.y0 < fb.y0 or bb.x1 > fb.x1 or bb.y1 > fb.y1: return None, f"text outside the figure: '{t.get_text()}'"
    return fig, None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--master", required=True); ap.add_argument("--eps", required=True, choices=["0.1", "0.05"]); ap.add_argument("--out", required=True)
    a = ap.parse_args(); out = Path(a.out); es = a.eps; eps = float(es)
    if out.exists(): return fail("output directory exists", out)
    if sha(a.master) != MASTER_SHA: return fail("master sha differs from the pin")
    if P.STYLE_SHA256 != STYLE_SHA: return fail("style module sha differs from the pin")
    V = json.load(open(a.master))["values"]; shown = {}
    try: res = load(V, es, shown)
    except KeyError as e: return fail("missing or non-finite key", e)
    if isinstance(res, str): return fail(res)
    titles, B, cov, (nsite, ndraw) = res
    fig, err = draw(eps, titles, B, cov)
    if fig is None: return fail(err)
    buf = P.render(fig); plt.close(fig)
    if buf["mediabox_mm"][1] > MAX_H_MM: return fail("figure taller than", MAX_H_MM, "mm:", buf["mediabox_mm"])
    all_cert = all(e["n_cert"] == e["n"] for e in B.values())
    caption_note = (f"Realized site FPR pooled over {nsite} synthetic sites of {ndraw} draws: boxes span the 25th to 75th percentiles, the black line marks the median, whiskers the "
                    "5th and 95th percentiles and the white dot the mean, without outlier points. The dashed line marks $\\epsilon$ and the shaded region marks realized FPR above "
                    "$\\epsilon$. The y axis does not start at zero. Numbers above the boxes are coverage. "
                    + ("Every draw certified in these four settings." if all_cert else "Boxes show certified draws only."))
    stem = f"{STEM}_eps{es}"
    payload = json.dumps(dict(script_sha256=sha(__file__), style_sha256=P.STYLE_SHA256, master=str(a.master), master_sha256=MASTER_SHA, eps=eps,
                              settings=[f"{d}|{ci}" for d, ci in SETTINGS], rules=[r for r, _ in RULES], all_certified=all_cert,
                              strips=[f"{d}|{ci}|{r}" for (d, ci, r), e in B.items() if e["n_cert"] < STRIP_N], caption_note=caption_note,
                              mediabox_mm=buf["mediabox_mm"], shown=shown), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    for ext in ("pdf", "png"):
        with open(out / f"{stem}.{ext}", "xb") as fh: fh.write(buf[ext])
    with open(out / f"{stem}_numbers.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh:
        for nm in (f"{stem}.pdf", f"{stem}.png", f"{stem}_numbers.json"): fh.write(f"{sha(out / nm)}  {nm}\n")
    print("written", out); return 0


if __name__ == "__main__": sys.exit(main())

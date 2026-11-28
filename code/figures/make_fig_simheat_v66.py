# P5a figure F-b (new, 实验管理 / user 2026-11-24 figure set, spec approved 实验管理 2026-11-25): coverage heat map of the simulation study over 50 synthetic sites per cell.
#   Replaces fig:sim (实验管理 2026-11-25).  WRITTEN, NOT RUN, until 实验管理 reviews it.  Full BJ width 139 mm (p5a_style_ieee v2.0, pinned; render() enforces the width),
#   pure black or white text.
# Data: only results_master keys (master 9b04913e, sim50 ingested by add_sim50_keys.py 71af83c7):
#   sim50/<d>|<ci>|<rule>|<eps>|all50|coverage|mean (cell text and colour), |all50|coverage|mcse (largest shown value -> caption), |all50|n_realisations (== 50);
#   protocol/delta; setting labels from e12/cells|<i>|cell|* and ebatch_e10/cells|<i>|factors|*; sim50/platt_fallback|e10|<ci>|n_tasks (dagger on the row label).
# r2 (实验管理 2026-11-25): transposed.  Two panels side by side, columns = the 5 rules in Table 3 order (headers above, rotation 0), rows = settings:
#   (a) single-factor design, 10 rows grouped by the varied factor (base | xi | mu | pi | tail); (b) stress design, 24 rows grouped in blocks
#   (base | rare fields | field sizes | slides | two-stage slides | small sites).  Row labels show the factor value ("$\xi$ = 0.25", "$\sigma$ = 0.5"); block labels at the
#   left of the row labels with a vertical bracket line; blocks separated by GAP in y; both panels share the row height.
# Encoding: diverging map centred at 1 - delta (white), SKY above (full at 1.00), VERM below (full at 1 - 2 delta and clipped below it; colour bar extend arrow 'min',
#   horizontal, under both panels).  Each cell shows its coverage to 2 decimals with the leading zero.
# r2: text colour by the higher WCAG contrast ratio against the cell colour: black if (L + 0.05) / 0.05 >= 1.05 / (L + 0.05), i.e. L >= 0.179, else white.
# Cell text size: 7 pt if every cell text fits its cell in display coordinates (width and height, 0.5 pt padding each side), else 6.5 pt, else exit 3 (nothing written).
# Gates (exit 3, nothing written): master sha MASTER_SHA; style sha STYLE_SHA; protocol/delta == 0.1; --eps in {0.1, 0.05}; for every shown (cell, rule): coverage mean and
#   mcse present and finite, 0 <= coverage <= 1, n_realisations == 50; e10 cell factors: each non-base cell differs from the base cell in exactly the factor of its block;
#   e12 cell names carry the block prefix of their position (base, A_, B_, C_, C2_, D_), two-stage cells have S = 12; r2: column headers do not overlap each other;
#   row labels do not overlap each other or the block labels; block labels, panel titles and column headers do not overlap; the colour bar does not overlap a panel;
#   PDF height <= MAX_H_MM.
# r2 (实验管理 2026-11-25): caption note follows the code (a draw without a certificate is evaluated at grid[0] of its candidate, m3_replay_lib.fixed_sequence), with
#   grid[0] named as in the main text (the descending grid tau_1 > ... > tau_G', so grid[0] = tau_1, the largest grid threshold).  The cannot-certify wording waits for D-c.
# r3 (实验管理 2026-11-25): (1) gates for the gap between the panels: every (b) row label, block label and bracket line is disjoint from panel (a)'s axes window
#   extent, and every label, header, title and colour-bar text of the figure lies inside the figure bbox; (2) caption sentence for the dagger from the code
#   (e10_synth.platt_synth: Newton-Raphson on the same two-parameter Platt objective when the frozen BFGS fit raises its 'platt fit' precision assertion; the
#   calibration enters only Split-Select, m3_split_dev_eval.run_site_seed); (3) draws per realisation read as all50|n / all50|n_realisations, equal for every shown
#   (cell, rule) and gated == DRAWS_EXPECTED, and the caption sentence built from it.
# r4 (实验管理 2026-11-25): (1) colour bar moved into the empty area under panel (a), a horizontal inset aligned with (a)'s cell columns (data x 0 .. 5) at CB_Y rows
#   below (a)'s last row; figure height FIG_H_MM 140 -> 130; gates: the colour bar (with its tick labels and label) does not overlap (a)'s cell area, (a)'s labels or
#   panel (b); (2) the redundant row label 'base' is left empty (the block label 'Base' names it); gate: exactly one empty row label per panel, in the Base block;
#   (3) two-stage block label 'Two-stage ($S$ = 12)' on two lines, both centred (one line widened (b)'s label column and the header gate failed).
# r5 (实验管理 2026-11-25): caption definition stated as in the protocol, 'a draw that cannot certify counting as covered', gated on simbox/uncert_total|<eps> == 0
#   (D-c d_c_simbox.py 17223230, ingested by add_simbox_keys.py c19380c2: no uncertified draw with realized FPR above eps in any of the 340 cells); fails if not 0.
#   Master pin 9b04913e -> dba06e5c (adds only simbox/*); the figure itself does not change.
# Output (new directory, 'x' mode): fig_simheat_v66_eps<eps>.pdf, .png (600 dpi), fig_simheat_v66_eps<eps>_numbers.json, SHA256SUMS.
# Usage: python3 -E -B make_fig_simheat_v66.py --master <tree/results_master.json> --eps 0.1 --out <new dir>
# public-repo patch (2026-11-28): the master is read through public_inputs.load_master, which accepts the public redacted results/results_master.json.gz
#   (uncompressed sha 4b296884) or the pinned internal master; non-public inputs stop with a plain message.  Drawing code unchanged.
import argparse, hashlib, json, math, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, to_rgb
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import matplotlib.transforms as mtrans
import p5a_style_ieee as P
import public_inputs as PI

MASTER_SHA = "dba06e5cf5c9f73b432eb2717bc542b990bb473ed16370f869207c1e6b0e12d1"
STYLE_SHA = "74417a86988854226a6d4967425d34b79b51c279b6cfc8d59c7942ca94e9f171"
# public-repo patch r2: STYLE_SHA follows the public style module, whose header line naming another project was redacted (original 84242189).
# r2: short column headers (实验管理 2026-11-25), two lines where a name is long
RULES = [("SPLIT6", "Split-\nSelect"), ("S3", "DE-\nGrid"), ("CB_c1s", "Cluster\nbootstrap"), ("GEE_c1s", "GEE"), ("BB_c1s", "Beta-\nbinomial")]
# (a) single-factor design: (block label, factor key in ebatch_e10/cells|i|factors, cells)
E10_BLOCKS = [("Base", None, [0]), (r"Shift $\xi$", "tau", [1, 2, 3]), (r"Cluster size $\mu$", "r", [4, 5]), (r"Positives $\pi$", "pi", [6, 7, 8]), ("Tail", "tail", [9])]
E10_SYM = {"tau": r"$\xi$", "r": r"$\mu$", "pi": r"$\pi$"}
# (b) stress design: (block label, name prefix, cells)
E12_BLOCKS = [("Base", "base", [0]), ("Rare fields", "A_", [1, 2, 3, 4, 5, 6]), ("Field sizes", "B_", [7, 8, 9, 10]), ("Slides", "C_", [11, 12, 13, 14, 15, 16]),
              ("Two-stage\n($S$ = 12)", "C2_", [17, 18, 19, 20]), ("Small sites", "D_", [21, 22, 23])]
GAP = 0.45
CELL_FS = (7.0, 6.5)
HEAD_FS = 7.0
PAD_PT = 0.5
FIG_H_MM = 130
CB_Y = 1.6
CB_H = 0.55
MAX_H_MM = 150
WCAG_L = 0.179
DRAWS_EXPECTED = 200
STEM = "fig_simheat_v66"


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def num(v): return f"{v:g}"


def lum(rgb):
    c = [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in rgb]
    return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]


def text_colour(rgb):
    # r2: higher WCAG contrast ratio; black wins iff (L + 0.05) / 0.05 >= 1.05 / (L + 0.05)
    L = lum(rgb)
    return P.BLACK if (L + 0.05) / 0.05 >= 1.05 / (L + 0.05) else P.WHITE


def load(V, es, shown):
    """Reads and gates every key; returns (delta, cols, dag, D, mcse) or an error string."""
    def g(k):
        if k not in V or V[k] is None: raise KeyError(k)
        shown[k] = V[k]; return V[k]
    delta = g("protocol/delta")
    if delta != 0.1: return "protocol/delta differs from 0.1"
    cols = {"e10": [], "e12": []}
    base10 = {f: g(f"ebatch_e10/cells|0|factors|{f}") for f in ("tau", "r", "pi", "tail")}
    for blk, fac, cells in E10_BLOCKS:
        for ci in cells:
            fv = {f: g(f"ebatch_e10/cells|{ci}|factors|{f}") for f in base10}
            diff = sorted(f for f in base10 if fv[f] != base10[f])
            if diff != ([] if fac is None else [fac]): return f"e10 cell {ci} does not vary only its block factor {fac}: {diff}"
            if fac == "tail" and fv[fac] != "t5": return f"e10 tail cell {ci} is not t5"
            lab = "" if fac is None else ("$t_5$" if fac == "tail" else f"{E10_SYM[fac]} = {num(fv[fac])}")
            cols["e10"].append((blk, ci, lab))
    for blk, pre, cells in E12_BLOCKS:
        for ci in cells:
            c = f"e12/cells|{ci}|cell|"; nm = g(c + "name")
            ok = nm == "base" if pre == "base" else (nm.startswith(pre) and not (pre == "C_" and nm.startswith("C2_")))
            if not ok: return f"e12 cell {ci} name {nm} does not match its block {pre}"
            if pre == "base": lab = ""
            elif pre == "A_": lab = rf"$f_\mathrm{{r}}$ = {num(g(c + 'rare|0'))}, $d_\mathrm{{r}}$ = {num(g(c + 'rare|1'))}"
            elif pre == "B_": lab = rf"$\sigma$ = {num(g(c + 'size_sigma'))}"
            elif pre == "C_": lab = rf"$S$ = {num(g(c + 'slides|0'))}, $\omega$ = {num(g(c + 'slides|1'))}"
            elif pre == "C2_":
                if g(c + "slides|0") != 12: return f"two-stage cell {ci} does not have S = 12"
                lab = rf"$\omega$ = {num(g(c + 'slides|1'))}, $s$ = {num(g(c + 'two_stage'))}"
            else: lab = rf"$F$ = {num(g(c + 'F'))}"
            cols["e12"].append((blk, ci, lab))
    dag = {}
    for d in ("e10", "e12"):
        for _, ci, _ in cols[d]:
            k = f"sim50/platt_fallback|{d}|{ci}|n_tasks"
            if k in V: dag[(d, ci)] = g(k) > 0
    D = {}; mcse = []; dpr = set()
    for d in ("e10", "e12"):
        for _, ci, _ in cols[d]:
            for r, _ in RULES:
                b = f"sim50/{d}|{ci}|{r}|{es}|all50|"
                cv, mc, nr, nn = g(b + "coverage|mean"), g(b + "coverage|mcse"), g(b + "n_realisations"), g(b + "n")
                if nr != 50: return f"n_realisations differs from 50: {b}"
                # r3: draws per realisation from the master, not hardcoded
                if nn % nr: return f"all50|n is not a multiple of n_realisations: {b}"
                dpr.add(nn // nr)
                if not (math.isfinite(cv) and math.isfinite(mc) and 0 <= cv <= 1): return f"coverage not finite or outside [0, 1]: {b}"
                D[(d, ci, r)] = cv; mcse.append(mc)
    if dpr != {DRAWS_EXPECTED}: return f"draws per realisation {sorted(dpr)} differ from {DRAWS_EXPECTED}"
    # r5: the caption's cannot-certify clause holds only if no uncertified draw exceeds eps in any cell at this eps
    if g(f"simbox/uncert_total|{es}") != 0: return f"simbox/uncert_total|{es} is not 0: the cannot-certify clause of the caption does not hold"
    return delta, cols, dag, D, mcse, DRAWS_EXPECTED


def draw(delta, cols, dag, D):
    """Draws the figure; returns (fig, used cell font) or (None, error string).  Writes nothing."""
    lo, mid = 1 - 2 * delta, 1 - delta
    cmap = LinearSegmentedColormap.from_list("cov", [to_rgb(P.VERM), (1, 1, 1), to_rgb(P.SKY)]); cmap.set_under(P.VERM)
    norm = TwoSlopeNorm(vmin=lo, vcenter=mid, vmax=1.0)
    P.apply()
    fig, axs = plt.subplots(1, 2, figsize=P.figsize(FIG_H_MM))
    ys_all, ymax = {}, 0.0
    for d, blocks in (("e10", E10_BLOCKS), ("e12", E12_BLOCKS)):
        ys, y, prev = [], 0.0, None
        for blk, ci, lab in cols[d]:
            if prev is not None and blk != prev: y += GAP
            ys.append(y); y += 1; prev = blk
        ys_all[d] = ys; ymax = max(ymax, y)
    texts, rowlabs, heads, blabels, ptitles, brackets = {}, {}, {}, [], [], []
    for ax, d, letter, title, blocks in ((axs[0], "e10", "a", "Single-factor design", E10_BLOCKS), (axs[1], "e12", "b", "Stress design", E12_BLOCKS)):
        ys = ys_all[d]; texts[d] = []
        for i, (blk, ci, lab) in enumerate(cols[d]):
            for j, (r, _) in enumerate(RULES):
                v = D[(d, ci, r)]; rgba = cmap(norm(v))
                ax.add_patch(Rectangle((j, ys[i]), 1, 1, facecolor=rgba, edgecolor=P.WHITE, lw=0.6))
                texts[d].append((ax.text(j + 0.5, ys[i] + 0.5, f"{v:.2f}", ha="center", va="center", fontsize=CELL_FS[0], color=text_colour(rgba[:3])), j, ys[i]))
        ax.set_xlim(0, len(RULES)); ax.set_ylim(ymax, 0)
        ax.set_yticks([yy + 0.5 for yy in ys]); ax.set_yticklabels([lab + ("$^{\\dagger}$" if dag.get((d, ci)) else "") for _, ci, lab in cols[d]])
        ax.xaxis.tick_top(); ax.set_xticks([j + 0.5 for j in range(len(RULES))]); ax.set_xticklabels([h for _, h in RULES], fontsize=HEAD_FS, linespacing=0.95)
        ax.tick_params(length=0, pad=2)
        for s in ax.spines.values(): s.set_visible(False)
        n0 = len(ax.texts); P.panel_label(ax, letter, title, pad_pt=21)
        ptitles.extend(ax.texts[n0:])
        rowlabs[d], heads[d] = ax.get_yticklabels(), ax.get_xticklabels()
    # r4: exactly one empty row label per panel, in the Base block
    for d in ("e10", "e12"):
        empty = [blk for blk, _, lab in cols[d] if lab == ""]
        if empty != ["Base"]: return None, f"empty row labels of {d} are not exactly the Base row: {empty}"
    sm = ScalarMappable(norm=norm, cmap=cmap)
    # r4: horizontal colour bar inside the empty area under panel (a), aligned with its cell columns
    y_cb = ys_all["e10"][-1] + 1 + CB_Y
    if y_cb + CB_H > ymax: return None, "no room for the colour bar under panel (a)"
    cax = axs[0].inset_axes([0, y_cb, len(RULES), CB_H], transform=axs[0].transData)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", extend="min", ticks=[lo, mid, 1.0])
    cb.ax.set_xticklabels([f"{lo:.2f}", f"{mid:.2f}", "1.00"]); cb.set_label("Coverage"); cb.outline.set_linewidth(P.LW_AXIS)
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    # block labels and bracket lines left of the widest row label of the panel; offsets in points from the axes edge, so they stay put under constrained layout
    for ax, d, blocks in ((axs[0], "e10", E10_BLOCKS), (axs[1], "e12", E12_BLOCKS)):
        wpt = max(t.get_window_extent(rend).width for t in rowlabs[d]) * 72 / fig.dpi + 4
        tr = mtrans.offset_copy(mtrans.blended_transform_factory(ax.transAxes, ax.transData), fig=fig, x=-wpt, y=0, units="points")
        ys, i0 = ys_all[d], 0
        for blk, _, cells in blocks:
            y0, y1 = ys[i0], ys[i0 + len(cells) - 1] + 1; i0 += len(cells)
            brackets.append(ax.add_line(Line2D([0, 0], [y0 + 0.1, y1 - 0.1], transform=tr, lw=0.6, color=P.BLACK, clip_on=False)))
            blabels.append(ax.annotate(blk, xy=(0, (y0 + y1) / 2), xycoords=("axes fraction", "data"), xytext=(-wpt - 3, 0), textcoords="offset points",
                                       ha="right", va="center", multialignment="center", fontsize=P.FS_TICK, color=P.BLACK, annotation_clip=False, linespacing=0.95))
    used = None
    pad = PAD_PT * fig.dpi / 72
    for fs in CELL_FS:
        for d in texts:
            for t, _, _ in texts[d]: t.set_fontsize(fs)
        fig.canvas.draw(); worst = []
        for d, ax in (("e10", axs[0]), ("e12", axs[1])):
            for t, x0, y0 in texts[d]:
                bb = t.get_window_extent(rend)
                (cx0, cy0), (cx1, cy1) = ax.transData.transform((x0, y0)), ax.transData.transform((x0 + 1, y0 + 1))
                cw, ch = abs(cx1 - cx0) - 2 * pad, abs(cy1 - cy0) - 2 * pad
                if bb.width > cw or bb.height > ch: worst.append((d, round(bb.width - cw, 2), round(bb.height - ch, 2)))
        if not worst: used = fs; break
    if used is None: return None, f"cell text does not fit its cell even at 6.5 pt: {worst[:5]}"
    fig.canvas.draw()

    def overlap(ts, what):
        bbs = [(t.get_text(), t.get_window_extent(rend)) for t in ts]
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                if bbs[i][1].overlaps(bbs[j][1]): return f"{what}: '{bbs[i][0]}' overlaps '{bbs[j][0]}'"
        return None
    for d in ("e10", "e12"):
        for ts, what in ((heads[d], f"column headers of {d}"), (list(rowlabs[d]) + [b for b in blabels if b.axes is (axs[0] if d == "e10" else axs[1])], f"row / block labels of {d}")):
            m = overlap(ts, what)
            if m: return None, m
    m = overlap(blabels + ptitles + list(heads["e10"]) + list(heads["e12"]), "block labels / panel titles / column headers")
    if m: return None, m
    # r3: (b)'s row labels, block labels and bracket lines stay out of panel (a)'s cell area
    a_box = axs[0].get_window_extent(rend)
    for art in list(rowlabs["e12"]) + [b for b in blabels if b.axes is axs[1]] + [l for l in brackets if l.axes is axs[1]]:
        if art.get_window_extent(rend).overlaps(a_box): return None, f"panel (b) label or bracket overlaps panel (a): {getattr(art, 'get_text', lambda: 'bracket line')()}"
    # r3: every text of the figure lies inside the figure bbox
    fb = fig.bbox
    alltext = list(rowlabs["e10"]) + list(rowlabs["e12"]) + list(heads["e10"]) + list(heads["e12"]) + blabels + ptitles + \
        list(cb.ax.get_xticklabels()) + [cb.ax.xaxis.label]
    for t in alltext:
        bb = t.get_window_extent(rend)
        if bb.x0 < fb.x0 or bb.y0 < fb.y0 or bb.x1 > fb.x1 or bb.y1 > fb.y1: return None, f"text outside the figure: '{t.get_text()}'"
    # r4: the colour bar (bar, tick labels, label) stays clear of (a)'s cell area, (a)'s labels and panel (b)
    cbb = cb.ax.get_tightbbox(rend)
    (ax0, ay0), (ax1, ay1) = axs[0].transData.transform((0, 0)), axs[0].transData.transform((len(RULES), ys_all["e10"][-1] + 1))
    a_cells = mtrans.Bbox.from_extents(min(ax0, ax1), min(ay0, ay1), max(ax0, ax1), max(ay0, ay1))
    if cbb.overlaps(a_cells): return None, "colour bar overlaps panel (a)'s cells"
    for t in list(rowlabs["e10"]) + [b for b in blabels if b.axes is axs[0]]:
        if t.get_text() and cbb.overlaps(t.get_window_extent(rend)): return None, f"colour bar overlaps label '{t.get_text()}'"
    if cbb.overlaps(axs[1].get_tightbbox(rend)): return None, "colour bar overlaps panel (b)"
    return fig, used


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--master", required=True); ap.add_argument("--eps", required=True, choices=["0.1", "0.05"]); ap.add_argument("--out", required=True)
    a = ap.parse_args(); out = Path(a.out); es = a.eps
    if out.exists(): return fail("output directory exists", out)
    if P.STYLE_SHA256 != STYLE_SHA: return fail("style module sha differs from the pin")
    M, MSHA = PI.load_master(a.master, (MASTER_SHA,)); V = M["values"]; shown = {}
    try: res = load(V, es, shown)
    except KeyError as e: return fail("missing key", e)
    if isinstance(res, str): return fail(res)
    delta, cols, dag, D, mcse, ndraw = res
    fig, used = draw(delta, cols, dag, D)
    if fig is None: return fail(used)
    buf = P.render(fig); plt.close(fig)
    if buf["mediabox_mm"][1] > MAX_H_MM: return fail("figure taller than", MAX_H_MM, "mm:", buf["mediabox_mm"])
    stem = f"{STEM}_eps{es}"
    # r5 (实验管理 2026-11-25): protocol wording with the cannot-certify clause, backed by the simbox/uncert_total gate in load()
    caption_note = ("Coverage is the share of draws whose certified threshold has realized site FPR at most $\\epsilon$, a draw that cannot certify counting as covered "
                    f"(protocol definition), over 50 synthetic sites of {ndraw} draws per cell. "
                    f"Colours are clipped at $1-2\\delta={1 - 2 * delta:.2f}$.")
    # r3: dagger sentence from e10_synth.platt_synth (the calibration enters only Split-Select)
    if any(dag.values()):
        caption_note += (" $\\dagger$ In at least one site realization of this setting, the frozen fitting routine for the Platt calibration $\\psi$ failed its precision "
                         "check on the synthetic source set, and $\\psi$ was fitted by Newton-Raphson on the same objective. $\\psi$ enters only Split-Select.")
    payload = json.dumps(dict(script_sha256=sha(__file__), style_sha256=P.STYLE_SHA256, master=str(a.master), master_sha256=MSHA, eps=float(es), delta=delta,
                              cell_fontsize_pt=used, rules=[r for r, _ in RULES], rows={d: [dict(block=b, cell=ci, label=l) for b, ci, l in cols[d]] for d in cols},
                              platt_dagger=[f"{d}|{ci}" for (d, ci), v in dag.items() if v], max_coverage_mcse=max(mcse), caption_note=caption_note,
                              mediabox_mm=buf["mediabox_mm"], shown=shown), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    for ext in ("pdf", "png"):
        with open(out / f"{stem}.{ext}", "xb") as fh: fh.write(buf[ext])
    with open(out / f"{stem}_numbers.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh:
        for nm in (f"{stem}.pdf", f"{stem}.png", f"{stem}_numbers.json"): fh.write(f"{sha(out / nm)}  {nm}\n")
    print("written", out, "cell font", used); return 0


if __name__ == "__main__": sys.exit(main())

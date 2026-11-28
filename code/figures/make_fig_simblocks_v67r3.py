# P5a figure F-g (new, draft 25 D7, 实验管理 2026-11-26 after the v6.21 scoring): coverage against the varied factor, one small panel per simulation block, plus the Aim 1
#   fixed-threshold check.  WRITTEN, NOT RUN, until 实验管理 reviews it.  No new run: only results_master keys.  Full BJ width 139 mm (p5a_style_ieee v2.0, pinned), pure black text,
#   colour only on marks (DE-Grid and Split-Select in their role colours as HERO series, cluster bootstrap, GEE and beta-binomial MUTED).
# Panels (eps = 0.10, the main target; one value per scenario = coverage over 50 synthetic sites, bar = +-1 MCSE):
#   (a) single-factor between-field shift xi: e10 scenarios with only tau varied (tau 0, 0.25, 0.5 = base, 1.0), x = xi;
#   (b) rare-field block: e12 1..6, x = (f_r, d_r) on two-line tick labels;            (c) field-size block: e12 0 (base, sigma 0) and 7..10, x = sigma;
#   (d) slide block: e12 11..16, x = (S, omega);                                       (e) two-stage block: e12 17..20, x = (omega, s), S = 12;
#   (f) small-site block: e12 21..23 and 0 (base, F = 300), x = F;
#   (g) Aim 1, fixed-threshold check: certification rate at the threshold with site FPR = eps in the ten single-factor scenarios (pooled over seeds), point estimate rho_hat (B1)
#       against the upper limit rho_U (B1u at level delta), with the dashed line at delta; x = scenario on two-line labels;
#   legend below the panels (all rules, the point estimate of (g) and the reference lines).
# Layout (r1, before any run): row 1 = (a)-(d); row 2 = (e), (f) and (g) spanning two columns, so that the ten scenario labels of (g) fit.
# Encoding: dashed black line at 1 - delta in (a)-(f); light VERM band below 1 - delta (failure region), as the other coverage figures; y grid only; x offsets of +-0.12 per rule so that
#   MCSE bars do not overlap.
# Data (keys only): sim50/<d>|<ci>|<rule>|0.1|all50|coverage|{mean, mcse} and |n_realisations (rules SPLIT6, S3, CB_c1s, GEE_c1s, BB_c1s); design values from
#   ebatch_e10/cells|<ci>|factors|* and e12/cells|<ci>|cell|*; Aim 1 from ebatch_e10_derived/cells|<ci>|fixed_threshold_pooled|0.1|1.0|certify_rate|{B1, B1u@0.10};
#   protocol/conf_level and protocol/delta for the reference lines.
# Gates (exit 3, nothing written): master sha MASTER_SHA (131b06d0, draft 23); style sha STYLE_SHA; every key present and finite; n_realisations equal for every shown entry;
#   design values match the block definitions (e10 scenarios vary only tau in (a); e12 names A_, B_, C_, C2_, D_ and 'base'); every mean +- mcse inside the y range;
#   tick labels within a panel and panel titles do not overlap; every text inside the figure bbox; PDF height <= MAX_H_MM.
# Output (new directory, 'x' mode): fig_simblocks_v67.pdf, .png (600 dpi), fig_simblocks_v67_numbers.json, SHA256SUMS.
# r1 (2026-11-26): bbox gate restricted to drawn y tick labels, after the first approved render stopped at that gate with nothing written.
# r2 (2026-11-26, 实验管理 review of the r1 PNG): (g) point estimate gets its own muted marker 'x' (was 'v', as GEE); (b) f_r printed in percent; (g) y label shortened, p = eps stated in the caption.
# r3 (2026-11-26, 实验管理): layout as two subfigures (fig.subfigures(2, 1)) with independent constrained layouts, so row 1 has four equal panels and the (g) y label no longer widens
#   the gap between (b) and (c); shared legend unchanged on the parent figure; new gate: the legend must not overlap any panel (axes with tick labels).
# Usage: python3 -E -B make_fig_simblocks_v67r3.py --master <tree/results_master.json> --out <new dir>
# public-repo patch (2026-11-28): the master is read through public_inputs.load_master, which accepts the public redacted results/results_master.json.gz
#   (uncompressed sha 4b296884) or the pinned internal master; non-public inputs stop with a plain message.  Drawing code unchanged.
import argparse, hashlib, json, math, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import p5a_style_ieee as P
import public_inputs as PI

# r2: muted grey of P.muted(0) with marker 'x'; 'x' has no face, so it gets a slightly larger size and edge width to stay visible at FS_MIN scale.
RHO_HAT_KW = {**P.muted(0, "x"), "markersize": 3.4, "markeredgewidth": 0.8}
MASTER_SHA = "131b06d0d0aa8ebeaba6f2d8c3fb102494106ca0182fbdb4d4b40c7913547489"
STYLE_SHA = "74417a86988854226a6d4967425d34b79b51c279b6cfc8d59c7942ca94e9f171"
# public-repo patch r2: STYLE_SHA follows the public style module, whose header line naming another project was redacted (original 84242189).
STEM = "fig_simblocks_v67"; FIG_H_MM = 104; MAX_H_MM = 110; EPS = "0.1"
RULES = [("SPLIT6", "Split-Select"), ("S3", "DE-Grid"), ("CB_c1s", "Cluster bootstrap"), ("GEE_c1s", "GEE"), ("BB_c1s", "Beta-binomial")]
DODGE = {"SPLIT6": -0.12, "S3": 0.12, "CB_c1s": -0.06, "GEE_c1s": 0.0, "BB_c1s": 0.06}
PANELS = [("a", "Between-field shift", "e10", [1, 2, 0, 3]), ("b", "Rare field", "e12", [1, 2, 3, 4, 5, 6]), ("c", "Field size", "e12", [0, 7, 8, 9, 10]),
          ("d", "Slides", "e12", [11, 12, 13, 14, 15, 16]), ("e", "Two-stage", "e12", [17, 18, 19, 20]), ("f", "Small site", "e12", [21, 22, 23, 0])]
PREFIX = {"b": "A_", "d": "C_", "e": "C2_", "f": "D_"}


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def num(V, k, shown):
    v = V[k]
    if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v): raise KeyError(k)
    shown[k] = v; return float(v)


def xlabel_of(V, letter, d, ci, shown):
    """Tick label and x axis title of one scenario (design values from the master)."""
    if d == "e10":
        if any(V[f"ebatch_e10/cells|{ci}|factors|{f}"] != V[f"ebatch_e10/cells|0|factors|{f}"] for f in ("r", "pi", "tail")): raise ValueError(f"e10 {ci} varies more than tau")
        return f"{num(V, f'ebatch_e10/cells|{ci}|factors|tau', shown):g}", r"$\xi$"
    nm = V[f"e12/cells|{ci}|cell|name"]; shown[f"e12/cells|{ci}|cell|name"] = nm
    if letter in PREFIX and ci != 0 and (not nm.startswith(PREFIX[letter]) or (letter == "d" and nm.startswith("C2_"))): raise ValueError(f"e12 {ci} ({nm}) not in block {letter}")
    if ci == 0 and nm != "base": raise ValueError("e12 0 is not the base scenario")
    if letter == "b": return f"{round(100 * num(V, f'e12/cells|{ci}|cell|rare|0', shown), 6):g}\n{num(V, f'e12/cells|{ci}|cell|rare|1', shown):g}", r"$f_{\mathrm{r}}$ (%) / $d_{\mathrm{r}}$"
    if letter == "c": return f"{num(V, f'e12/cells|{ci}|cell|size_sigma', shown):g}", r"$\sigma$"
    if letter == "d": return f"{num(V, f'e12/cells|{ci}|cell|slides|0', shown):g}\n{num(V, f'e12/cells|{ci}|cell|slides|1', shown):g}", r"$S$ / $\omega$"
    if letter == "e": return f"{num(V, f'e12/cells|{ci}|cell|slides|1', shown):g}\n{num(V, f'e12/cells|{ci}|cell|two_stage', shown):g}", r"$\omega$ / $s$"
    return f"{num(V, f'e12/cells|{ci}|cell|F', shown):g}", r"$F$"


def load(V, shown):
    conf, delta = num(V, "protocol/conf_level", shown), num(V, "protocol/delta", shown)
    D, labs, nreal = {}, {}, set()
    for letter, _, d, cis in PANELS:
        for ci in cis:
            labs[(letter, ci)] = xlabel_of(V, letter, d, ci, shown)
            for r, _ in RULES:
                base = f"sim50/{d}|{ci}|{r}|{EPS}|all50|"
                D[(d, ci, r)] = (num(V, base + "coverage|mean", shown), num(V, base + "coverage|mcse", shown)); nreal.add(num(V, base + "n_realisations", shown))
    if len(nreal) != 1: return f"n_realisations differ: {sorted(nreal)}"
    A1 = {}
    for ci in range(10):
        b = f"ebatch_e10_derived/cells|{ci}|fixed_threshold_pooled|0.1|1.0|certify_rate|"
        A1[ci] = (num(V, b + "B1", shown), num(V, b + "B1u@0.10", shown))
    a1lab = {}
    for ci in range(10):
        f = {k: V[f"ebatch_e10/cells|{ci}|factors|{k}"] for k in ("tau", "r", "pi", "tail")}
        b0 = {k: V[f"ebatch_e10/cells|0|factors|{k}"] for k in ("tau", "r", "pi", "tail")}
        diff = [k for k in f if f[k] != b0[k]]
        if len(diff) > 1: return f"e10 {ci} varies more than one factor"
        a1lab[ci] = "base" if not diff else {"tau": r"$\xi$", "r": r"$\mu$", "pi": r"$\pi$", "tail": "tail"}[diff[0]] + "\n" + (f"{f[diff[0]]:g}" if diff[0] != "tail" else str(f["tail"]))
    return conf, delta, D, labs, nreal.pop(), A1, a1lab


def draw(conf, delta, D, labs, A1, a1lab):
    P.apply()
    fig = plt.figure(figsize=P.figsize(FIG_H_MM)); top, bot = fig.subfigures(2, 1)
    M = {**top.subplot_mosaic([["a", "b", "c", "d"]]), **bot.subplot_mosaic([["e", "f", "g", "g"]])}; axs = [M[k] for k in "abcdefg"]
    lo = min(m - s for m, s in D.values()); ylo = math.floor(min(lo, 1 - 2 * delta) * 20) / 20; yhi = 1.005
    texts = []
    for ax, (letter, title, d, cis) in zip(axs[:6], PANELS):
        xs = list(range(len(cis)))
        P.band(ax, ylo, conf, axis="y", color=P.VERM, alpha=P.BAND_ALPHA); P.ref_line(ax, conf, axis="y", color=P.BLACK)
        for i, (r, _) in enumerate(RULES):
            kw = P.hero(P.ROLE[r]["color"], P.ROLE[r]["marker"]) if r in P.ROLE else P.muted(i - 2)
            ys = [D[(d, ci, r)][0] for ci in cis]; es = [D[(d, ci, r)][1] for ci in cis]
            ax.errorbar([x + DODGE[r] for x in xs], ys, yerr=es, capsize=0, elinewidth=0.6, **kw)
        ax.set_ylim(ylo, yhi); ax.set_xlim(-0.5, len(cis) - 0.5); P.ygrid(ax)
        ax.set_xticks(xs); ax.set_xticklabels([labs[(letter, ci)][0] for ci in cis], fontsize=P.FS_MIN, linespacing=0.95); ax.tick_params(axis="x", length=0)
        ax.set_xlabel(labs[(letter, cis[0])][1], fontsize=P.FS_MIN, labelpad=1)
        if letter in "ae": ax.set_ylabel("Coverage")
        else: ax.tick_params(axis="y", labelleft=False)
        n0 = len(ax.texts); P.panel_label(ax, letter, title); texts.extend(ax.texts[n0:])
    ax = axs[6]; xs = list(range(10))
    P.ref_line(ax, delta, axis="y", color=P.BLACK)
    ax.plot(xs, [A1[ci][0] for ci in xs], **RHO_HAT_KW); ax.plot(xs, [A1[ci][1] for ci in xs], **P.hero(P.ROLE["S3"]["color"], P.ROLE["S3"]["marker"]))
    ax.set_xticks(xs); ax.set_xticklabels([a1lab[ci] for ci in xs], fontsize=P.FS_MIN, linespacing=0.9); ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, max(max(v) for v in A1.values()) * 1.15); P.ygrid(ax); ax.set_ylabel("Certification rate")
    n0 = len(ax.texts); P.panel_label(ax, "g", "Fixed threshold"); texts.extend(ax.texts[n0:])
    hs = [Line2D([], [], **(P.hero(P.ROLE[r]["color"], P.ROLE[r]["marker"]) if r in P.ROLE else P.muted(i - 2))) for i, (r, _) in enumerate(RULES)]
    hs += [Line2D([], [], **RHO_HAT_KW), Line2D([], [], color=P.BLACK, lw=P.LW_REF, ls=P.DASH_REF)]
    names = [n for _, n in RULES] + [r"(g) point estimate $\hat\rho$", r"$1-\delta$, (g) $\delta$"]
    leg = fig.legend(hs, names, loc="outside lower center", ncol=4, frameon=False, fontsize=P.FS_MIN, handlelength=1.8, columnspacing=1.2)
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    for (d, ci, r), (m, s) in D.items():
        if not (ylo <= m - s and m + s <= yhi + 1e-9): return None, f"coverage bar outside the y range: {d}|{ci}|{r}"
    def overlap(ts, what):
        bbs = [(t.get_text(), t.get_window_extent(rend)) for t in ts if t.get_text()]
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                if bbs[i][1].overlaps(bbs[j][1]): return f"{what}: '{bbs[i][0]}' overlaps '{bbs[j][0]}'"
        return None
    for ax in axs[:7]:
        m = overlap(ax.get_xticklabels(), "x tick labels")
        if m: return None, m
    m = overlap(texts, "panel titles")
    if m: return None, m
    fb = fig.bbox
# r3: the shared legend sits on the parent figure below both subfigures; it must not overlap any panel including its tick and axis labels.
    lb = leg.get_window_extent(rend)
    for ax in axs[:7]:
        if lb.overlaps(ax.get_tightbbox(rend)): return None, "legend overlaps a panel"
# r1: matplotlib keeps labels for ticks outside the view limits but does not draw them, so the bbox gate checks only y tick labels whose tick lies inside the y limits.
    def drawn_y(ax):
        lo, hi = ax.get_ylim(); return [t for t in ax.get_yticklabels() if lo - 1e-9 <= t.get_position()[1] <= hi + 1e-9]
    for t in texts + [t for ax in axs[:7] for t in ax.get_xticklabels() + drawn_y(ax)] + [ax.xaxis.label for ax in axs[:7]] + [ax.yaxis.label for ax in axs[:7]] + leg.get_texts():
        if not t.get_text(): continue
        bb = t.get_window_extent(rend)
        if bb.x0 < fb.x0 or bb.y0 < fb.y0 or bb.x1 > fb.x1 or bb.y1 > fb.y1: return None, f"text outside the figure: '{t.get_text()}'"
    return fig, None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--master", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(); out = Path(a.out)
    if out.exists(): return fail("output directory exists", out)
    if P.STYLE_SHA256 != STYLE_SHA: return fail("style module sha differs from the pin")
    M, MSHA = PI.load_master(a.master, (MASTER_SHA,)); V = M["values"]; shown = {}
    try: res = load(V, shown)
    except (KeyError, ValueError) as e: return fail("key or design check", e)
    if isinstance(res, str): return fail(res)
    conf, delta, D, labs, nreal, A1, a1lab = res
    fig, err = draw(conf, delta, D, labs, A1, a1lab)
    if fig is None: return fail(err)
    buf = P.render(fig); plt.close(fig)
    if buf["mediabox_mm"][1] > MAX_H_MM: return fail("figure taller than", MAX_H_MM, "mm:", buf["mediabox_mm"])
    payload = json.dumps(dict(script_sha256=sha(__file__), style_sha256=P.STYLE_SHA256, master=str(a.master), master_sha256=MSHA, eps=float(EPS), n_realisations=nreal,
                              panels=[(l, t, d, c) for l, t, d, c in PANELS], rules=[r for r, _ in RULES], mediabox_mm=buf["mediabox_mm"], shown=shown), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    for ext in ("pdf", "png"):
        with open(out / f"{STEM}.{ext}", "xb") as fh: fh.write(buf[ext])
    with open(out / f"{STEM}_numbers.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh:
        for nm in (f"{STEM}.pdf", f"{STEM}.png", f"{STEM}_numbers.json"): fh.write(f"{sha(out / nm)}  {nm}\n")
    print("written", out); return 0


if __name__ == "__main__": sys.exit(main())

# P5a figure F-e (new, 实验管理 / user 2026-11-24 figure set): forest plot of the confirmation sites (PlasmoCount, IML Malaria, M5) at both targets.
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  Width 139 mm, IEEE style module p5a_style_ieee v2.0 (pinned), text >= 7 pt, pure black text.
# Data: results_master keys only; the D-e cells are mapped by (site, rule, eps, stat), never by index, through the readable index of the D-e source file
#   <tree>/results_src/deconfirm/d_e_confirm_intervals.json (sha DE_SRC_SHA, registered in files_sha256 of the master):
#   IML Malaria, DE-Grid (S3) and the exact rule with two labeled negatives per field (R1_2), eps 0.10 and 0.05: d_e/cells|<i>|point with the post hoc two-sided 90%
#     seed-cluster interval d_e/cells|<i>|ci90_seedcluster|0/1;
#   PlasmoCount (decision 2a, 实验管理 2026-11-25: per-draw values were not retained; no interval, points only): rev1_results/agg|a2|plasmo|22|<rule>|<eps>|cov / ds_eps;
#   M5, Split-Select at eps 0.10: coverage m5/cells_p|SPLIT6|0.1|coverage with the stored 90% seed-cluster interval |coverage_ci90_seedcluster|0/1;
#     Delta Sens_eps m5/cells_p|SPLIT6|0.1|ds_eps as a point (no interval, decision 2a).
#   PlasmoCount R1_2 at eps 0.05 is shown in (b) only: it could not certify in any draw (rev1_results/agg|a2|plasmo|22|R1_2|0.05|cannot_certify == 1, gated), so a
#     coverage point would describe the fallback threshold.
# Layout: (a) coverage, nominal 1 - delta dashed, light blue / light vermillion bands; (b) Delta Sens_eps, grey zero line.  Same rows in both panels, grouped by site;
#   filled markers eps = 0.10, open markers eps = 0.05; DE-Grid and Split-Select HERO, the exact rule MUTED.  No number in the figure.
# Caption text required by 实验管理: "intervals in this figure are two-sided 90% seed-cluster intervals (post hoc); Section 5.6 reports the preregistered Wilson and
#   seed-level t intervals"; PlasmoCount and M5 Delta Sens_eps: "per-draw values were not retained; no interval".  CAPTION_NOTE is written to the numbers json.
# Gates (exit 3, nothing written): master sha; style sha; D-e source sha and its files_sha256 entry; protocol eps_main 0.1 / eps_secondary 0.05 / delta 0.1; each D-e
#   cell found exactly once by (site, rule, eps, stat), its d_e/cells|<i>|eps equal to eps, and |point - master[master_key]| <= 1e-12; lo <= point <= hi for every interval.
# Output (new directory, 'x' mode): fig_confforest_v65.pdf, .png (600 dpi), fig_confforest_v65_numbers.json, SHA256SUMS.
# Usage: python3 -E -B make_fig_confforest_v65.py --master <tree/results_master.json> --de <tree/results_src/deconfirm/d_e_confirm_intervals.json> --out <new dir>
# public-repo patch (2026-11-28): the master is read through public_inputs.load_master, which accepts the public redacted results/results_master.json.gz
#   (uncompressed sha 4b296884) or the pinned internal master; non-public inputs stop with a plain message.  Drawing code unchanged.
import argparse, hashlib, json, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import p5a_style_ieee as P
import public_inputs as PI

MASTER_SHA = "09e168643816c024fccf77fc7306e3b7cb6cd88ea3cb118c3c4d16bd3f9ebeb5"
STYLE_SHA = "74417a86988854226a6d4967425d34b79b51c279b6cfc8d59c7942ca94e9f171"
# public-repo patch r2: STYLE_SHA follows the public style module, whose header line naming another project was redacted (original 84242189).
DE_SRC_SHA = "21dd87da35a552a950dedb34ff866ac5e87ea3bb811e6c04b1167e69c1c293a5"
DE_REL = "results_src/deconfirm/d_e_confirm_intervals.json"
STEM = "fig_confforest_v65"
# r2 (实验管理 2026-11-25): one more sentence for the "not certified" cell
CAPTION_NOTE = ("Intervals in this figure are two-sided 90% seed-cluster intervals (post hoc); Section 5.6 reports the preregistered Wilson and seed-level t intervals. "
                "PlasmoCount (all values) and M5 Delta Sens_eps: per-draw values were not retained; no interval. "
                "On PlasmoCount at eps = 0.05 the exact rule with two labeled negatives per field certified in no draw, so it has no coverage point in (a); "
                "its Delta Sens_eps in (b), marked †, is evaluated at the fallback threshold (the largest labeled negative score), as for the † entries of Table S7.1, "
                "which also report the fallback coverage. "
                "Intervals narrower than the marker are hidden by it.")
# r3 (实验管理 2026-11-25): the CAPTION_NOTE above gains the fallback half-sentence and the hidden-interval sentence; row names as in F-d (capitalised rule,
#   site as a bold group title above its rows); "not certified" on a white box
NAME = {"S3": "DE-Grid", "R1_2": "Exact rule, 2 negatives/field", "SPLIT6": "Split-Select"}
# rows: (site label, rule, eps string, source of (a), source of (b)); source = ("de", site) | ("plasmo",) | ("m5",) | None
ROWS = [("PlasmoCount", "S3", "0.1", ("plasmo",), ("plasmo",)), ("PlasmoCount", "S3", "0.05", ("plasmo",), ("plasmo",)),
        ("PlasmoCount", "R1_2", "0.1", ("plasmo",), ("plasmo",)), ("PlasmoCount", "R1_2", "0.05", None, ("plasmo",)),
        ("IML Malaria", "S3", "0.1", ("de", "iml"), ("de", "iml")), ("IML Malaria", "S3", "0.05", ("de", "iml"), ("de", "iml")),
        ("IML Malaria", "R1_2", "0.1", ("de", "iml"), ("de", "iml")), ("IML Malaria", "R1_2", "0.05", ("de", "iml"), ("de", "iml")),
        ("M5", "SPLIT6", "0.1", ("m5",), ("m5",))]
TOL = 1e-12


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def main():
    ap = argparse.ArgumentParser()
    for k in ("--master", "--de", "--out"): ap.add_argument(k, required=True)
    a = ap.parse_args(); out = Path(a.out)
    if out.exists(): return fail("output directory exists", out)
    if sha(a.de) != DE_SRC_SHA: return fail("D-e source sha differs from the pin")
    if P.STYLE_SHA256 != STYLE_SHA: return fail("style module sha differs from the pin")
    M, MSHA = PI.load_master(a.master, (MASTER_SHA,)); V = M["values"]; shown = {}
    if M.get("files_sha256", {}).get(DE_REL) != DE_SRC_SHA: return fail("D-e source not registered in files_sha256 with the pinned sha")
    cells = json.load(open(a.de))["cells"]

    def g(k):
        if k not in V or V[k] is None: raise KeyError(k)
        shown[k] = V[k]; return V[k]

    def de(site, rule, eps, stat):
        hit = [i for i, c in enumerate(cells) if (c["site"], c["rule"], c["eps"], c["stat"]) == (site, rule, float(eps), stat)]
        if len(hit) != 1: raise ValueError(f"D-e cell {site} {rule} {eps} {stat} found {len(hit)} times")
        i = hit[0]; b = f"d_e/cells|{i}|"
        if g(b + "eps") != float(eps): raise ValueError(f"d_e/cells|{i}|eps differs from {eps}")
        pt = g(b + "point")
        if abs(pt - g(cells[i]["master_key"])) > TOL: raise ValueError(f"D-e point differs from {cells[i]['master_key']}")
        return pt, g(b + "ci90_seedcluster|0"), g(b + "ci90_seedcluster|1")
    try:
        if g("protocol/eps_main") != 0.1 or g("protocol/eps_secondary") != 0.05 or g("protocol/delta") != 0.1: return fail("protocol eps / delta differ from 0.1 / 0.05 / 0.1")
        delta = g("protocol/delta")
        if g("rev1_results/agg|a2|plasmo|22|R1_2|0.05|cannot_certify") != 1: return fail("PlasmoCount R1_2 at 0.05 certified in some draw; the omitted coverage point must be reconsidered")
        A, B = [], []
        for site, rule, es, sa, sb in ROWS:
            for src, st, lst in ((sa, "cov", A), (sb, "dsens", B)):
                if src is None: lst.append(None); continue
                if src[0] == "de": v = de(src[1], rule, es, st)
                elif src[0] == "plasmo": v = (g(f"rev1_results/agg|a2|plasmo|22|{rule}|{es}|{'cov' if st == 'cov' else 'ds_eps'}"), None, None)
                elif st == "cov": v = (g("m5/cells_p|SPLIT6|0.1|coverage"), g("m5/cells_p|SPLIT6|0.1|coverage_ci90_seedcluster|0"), g("m5/cells_p|SPLIT6|0.1|coverage_ci90_seedcluster|1"))
                else: v = (g("m5/cells_p|SPLIT6|0.1|ds_eps"), None, None)
                if v[1] is not None and not v[1] <= v[0] <= v[2]: return fail("point outside its interval", site, rule, es, st)
                lst.append(v)
    except (KeyError, ValueError) as e: return fail(e)
    P.apply()
    fig, (axA, axB) = plt.subplots(1, 2, figsize=P.figsize(80), sharey=True)
    # r3: a group-title row (site name, bold) above the rows of each site; ticks = titles + rows
    ys, y, prev, tick_y, tick_lab, title_y = [], 0.0, None, [], [], []
    for site, rule, es, *_ in ROWS:
        if site != prev:
            if prev is not None: y += 0.3
            tick_y.append(y); tick_lab.append(site); title_y.append(y); y += 1
        ys.append(y); prev = site; y += 1
        tick_y.append(ys[-1]); tick_lab.append(f"{NAME[rule]}, $\\epsilon={'0.10' if es == '0.1' else '0.05'}$")
    drawn = []
    for ax, vals, st in ((axA, A, "cov"), (axB, B, "dsens")):
        nums = [x for v in vals if v for x in v if x is not None]
        if st == "cov":
            # r3: up to 1.02 as in F-d
            xl = (min(nums + [1 - delta]) - 0.02, 1.02)
            P.band(ax, 1 - delta, xl[1], axis="x", color=P.SKY, alpha=0.16); P.band(ax, xl[0], 1 - delta, axis="x", color=P.VERM, alpha=0.08)
            P.ref_line(ax, 1 - delta, axis="x", color=P.BLACK)
        else:
            xl = (min(nums + [0]) - 0.03, max(nums + [0]) + 0.03)
            P.ref_line(ax, 0, axis="x", color=P.MUTED[0], dashed=False, lw=P.LW_AXIS)
        for yy, (site, rule, es, *_), v in zip(ys, ROWS, vals):
            if v is None:
                # r2 (实验管理 2026-11-25): mark the empty cell; black 7 pt (user rule: figure text pure black or white)
                # r3: on a white box so it does not read as a result inside the blue (nominal) band
                ax.text(xl[1] - 0.005, yy, "not certified", ha="right", va="center", fontsize=P.FS_MIN, color=P.BLACK, zorder=7,
                        bbox=dict(facecolor=P.WHITE, edgecolor="none", pad=1.0))
                continue
            pt, lo, hi = v
            if rule in P.ROLE: c, mk, ms, lw = P.ROLE[rule]["color"], P.ROLE[rule]["marker"], P.HERO_KW["markersize"], 1.2
            else: c, mk, ms, lw = P.MUTED[6], "D", 3.0, 0.8
            face = c if es == "0.1" else "white"
            kw = dict(fmt=mk, color=c, ecolor=c, elinewidth=lw, capsize=1.5, ms=ms, mfc=face, mec=c, mew=0.8, zorder=6)
            if lo is not None:
                ax.errorbar([pt], [yy], xerr=[[pt - lo], [hi - pt]], **kw); drawn.append((ax, f"{site}|{rule}|{es}|{st}", pt, lo, hi, yy, ms))
            else: ax.errorbar([pt], [yy], **kw)
            # r4 (实验管理 2026-11-25): the (b) point of a row with no certified draw (no coverage point in (a)) carries a black dagger, as in Table S7.1
            if st == "dsens" and A[ys.index(yy)] is None:
                ax.annotate("†", xy=(pt, yy), xytext=(ms, 0), textcoords="offset points", ha="left", va="center", fontsize=P.FS_MIN, color=P.BLACK, zorder=7)
        ax.set_xlim(*xl); ax.set_ylim(ys[-1] + 0.6, -0.6); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    axA.set_yticks(tick_y); axA.set_yticklabels(tick_lab)
    for t_, yv in zip(axA.get_yticklabels(), tick_y):
        if yv in title_y: t_.set_fontweight("bold")
    fig.canvas.draw(); hidden = []
    for ax, key, pt, lo, hi, yy, ms in drawn:
        x0, x1, xp = (ax.transData.transform((x, yy))[0] for x in (lo, hi, pt))
        if max(xp - x0, x1 - xp) <= ms * fig.dpi / 72 / 2: hidden.append(key)
    axA.set_xlabel("Coverage"); axB.set_xlabel(r"$\Delta\mathrm{Sens}_\epsilon$")
    P.panel_label(axA, "a", "Coverage"); P.panel_label(axB, "b", r"$\Delta\mathrm{Sens}_\epsilon$")
    buf = P.render(fig); plt.close(fig)
    payload = json.dumps(dict(script_sha256=sha(__file__), style_sha256=P.STYLE_SHA256, master=str(a.master), master_sha256=MSHA, de_source_sha256=DE_SRC_SHA,
                              caption_note=CAPTION_NOTE, rows=[r[:3] for r in ROWS], hidden_intervals=hidden, n_intervals=len(drawn),
                              mediabox_mm=buf["mediabox_mm"], shown=shown), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    for ext in ("pdf", "png"):
        with open(out / f"{STEM}.{ext}", "xb") as fh: fh.write(buf[ext])
    with open(out / f"{STEM}_numbers.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh:
        for nm in (f"{STEM}.pdf", f"{STEM}.png", f"{STEM}_numbers.json"): fh.write(f"{sha(out / nm)}  {nm}\n")
    print("written", out); return 0


if __name__ == "__main__": sys.exit(main())

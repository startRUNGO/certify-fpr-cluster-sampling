# P5a figure F-d (new, 实验管理 / user 2026-11-24 figure set): forest plot of every rule on the three development sites at the main target.
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  Width 139 mm, IEEE style module p5a_style_ieee v2.0 (pinned), text >= 7 pt, pure black text.
# Data: only results_master keys d_d/<rule>|<site>|<eps>|<field> (d_d_dev_intervals.py fa9902c5, VM run ddintv_20261124T232236Z, ingested by add_ddintv_keys.py 0a5fd9ed):
#   coverage and its 90% seed-cluster interval |coverage_ci90_seedcluster|0/1; Delta Sens_eps and |ds_eps_ci90_seedcluster|0/1; eps = protocol/eps_main.
# Layout: 2 x 3 panels; columns = BBBC-ALL, MP-IDB, Lacuna-B; row 1 coverage with the nominal 1 - delta (protocol/delta) as a dashed line, light blue band above it and
#   light vermillion band below; row 2 Delta Sens_eps with a solid grey zero line.  Rows of each panel = the 12 rules in the order of Table 3 (tab:dev); DE-Grid and
#   Split-Select in their role colours and markers (HERO), the other rules grey (MUTED).  Row 2 shows only the rules for which Table 3 reports Delta Sens_eps
#   (DSENS_RULES); the cluster-aware and cell-level risk-control rows stay empty there, as in Table 3 ("--").
# r2 (实验管理 2026-11-25): row 2 shows all 12 rules; rule labels only in the first column (shared y axis); text >= 7 pt (6.5 pt allowed for the labels if needed).
# Gates (exit 3, nothing written): master sha MASTER_SHA; style sha STYLE_SHA; protocol/eps_main == 0.1 and protocol/delta == 0.1 (the axis text says 0.90); for every
#   shown cell: d_d/...|eps equals the eps of the key path, lo <= point <= hi for both intervals, 5 per-seed values present; no number is printed in the figure.
# Output (new directory, 'x' mode): fig_devforest_v65.pdf, .png (600 dpi), fig_devforest_v65_numbers.json, SHA256SUMS.
# Usage: python3 -E -B make_fig_devforest_v65.py --master <tree/results_master.json> --out <new dir>
import argparse, hashlib, json, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import p5a_style_ieee as P

MASTER_SHA = "09e168643816c024fccf77fc7306e3b7cb6cd88ea3cb118c3c4d16bd3f9ebeb5"
STYLE_SHA = "8424218990b9bec4dcb8b3bf6c54d7b2da5671b2fc9444511ccbefa7edd39ce3"
STEM = "fig_devforest_v65"
SITES = [("bbbc041_all", "BBBC-ALL"), ("mpidb_full", "MP-IDB"), ("lacuna_B", "Lacuna-B")]
# (rule key, label as in Table 3); order of Table 3
RULES = [("SPLIT6", "Split-Select"), ("S3", "DE-Grid"), ("R1_2", "Exact rule, 2 negatives/field"), ("B1", r"Design-effect (point ICC), $c^0$"),
         ("B1u_c0", r"Design-effect bound, $c^0$"), ("B1u_c1s", r"Design-effect bound, $c^{1s}$"), ("CB_c1s", "Cluster bootstrap"), ("GEE_c1s", "GEE cluster-robust"),
         ("BB_c1s", "Beta-binomial"), ("R1_all", "Exact rule, pooled negatives"), ("SSRCPS_WSR", "Risk control per cell, WSR bound"), ("SSRCPS_CLT", "Risk control per cell, CLT bound")]
# r2 (实验管理 2026-11-25): Delta Sens_eps for all 12 rules (d_d has ds_eps and its interval for the cluster-aware and risk-control rules too)
DSENS_RULES = {r for r, _ in RULES}
# r3 (实验管理 2026-11-25): caption sentence for intervals hidden by their marker (all intervals exist; hidden_intervals in the numbers json lists them)
CAPTION_NOTE = "Every point has a two-sided 90% seed-cluster interval; intervals narrower than the marker are hidden by it."


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--master", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(); out = Path(a.out)
    if out.exists(): return fail("output directory exists", out)
    if sha(a.master) != MASTER_SHA: return fail("master sha differs from the pin")
    if P.STYLE_SHA256 != STYLE_SHA: return fail("style module sha differs from the pin")
    V = json.load(open(a.master))["values"]; shown = {}

    def g(k):
        if k not in V or V[k] is None: raise KeyError(k)
        shown[k] = V[k]; return V[k]
    try:
        eps, delta = g("protocol/eps_main"), g("protocol/delta")
        if eps != 0.1 or delta != 0.1: return fail("protocol/eps_main or protocol/delta differ from 0.1")
        es = "0.1"; D = {}
        for r, _ in RULES:
            for s, _ in SITES:
                b = f"d_d/{r}|{s}|{es}|"
                if g(b + "eps") != eps: return fail("eps field differs from the key path", b)
                for st, ci in (("coverage", "coverage_ci90_seedcluster"), ("ds_eps", "ds_eps_ci90_seedcluster")):
                    if st == "ds_eps" and r not in DSENS_RULES: continue
                    pt, lo, hi = g(b + st), g(b + ci + "|0"), g(b + ci + "|1")
                    ps = [g(b + f"per_seed_{'coverage' if st == 'coverage' else 'ds_eps'}|{i}") for i in range(5)]
                    if not lo <= pt <= hi: return fail("point outside its interval", b + st)
                    # r3: the caption says every point has an interval, so a zero-width interval fails
                    if not lo < hi: return fail("zero-width interval", b + st)
                    D[(r, s, st)] = (pt, lo, hi)
    except KeyError as e: return fail("missing key", e)
    P.apply()
    fig, axs = plt.subplots(2, 3, figsize=P.figsize(118), sharey=True, height_ratios=[1, 1])
    ys = list(range(len(RULES)))
    HK = dict(ms=P.HERO_KW["markersize"], mec="white", mew=0.5, zorder=6)
    cov_vals = [v for (r, s, st), t in D.items() if st == "coverage" for v in t]
    ds_vals = [v for (r, s, st), t in D.items() if st == "ds_eps" for v in t]
    # r3 (实验管理 2026-11-25): coverage axis up to 1.02 so the markers at 1.0 do not touch the frame; the blue band follows the new limit
    xc = (min(cov_vals + [1 - delta]) - 0.02, 1.02); xd = (min(ds_vals + [0]) - 0.03, max(ds_vals + [0]) + 0.03)
    drawn = []
    for j, (s, sname) in enumerate(SITES):
        for i_row, (st, xl) in enumerate((("coverage", xc), ("ds_eps", xd))):
            ax = axs[i_row, j]
            if st == "coverage":
                P.band(ax, 1 - delta, xl[1], axis="x", color=P.SKY, alpha=0.16); P.band(ax, xl[0], 1 - delta, axis="x", color=P.VERM, alpha=0.08)
                P.ref_line(ax, 1 - delta, axis="x", color=P.BLACK)
            else: P.ref_line(ax, 0, axis="x", color=P.MUTED[0], dashed=False, lw=P.LW_AXIS)
            for y, (r, _) in zip(ys, RULES):
                if (r, s, st) not in D: continue
                pt, lo, hi = D[(r, s, st)]
                if r in P.ROLE:
                    c, mk = P.ROLE[r]["color"], P.ROLE[r]["marker"]
                    ax.errorbar([pt], [y], xerr=[[pt - lo], [hi - pt]], fmt=mk, color=c, ecolor=c, elinewidth=1.2, capsize=1.5, **HK); ms = HK["ms"]
                else:
                    c = P.MUTED[6]; ms = 2.8
                    ax.errorbar([pt], [y], xerr=[[pt - lo], [hi - pt]], fmt="D", color=c, ecolor=c, elinewidth=0.8, capsize=1.2, ms=ms, zorder=3)
                drawn.append((ax, f"{r}|{s}|{st}", pt, lo, hi, y, ms))
            ax.set_xlim(*xl); ax.set_ylim(len(RULES) - 0.5, -0.6); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
            ax.set_xlabel("Coverage" if st == "coverage" else r"$\Delta\mathrm{Sens}_\epsilon$")
            P.panel_label(ax, "abcdef"[3 * i_row + j], sname)
    axs[0, 0].set_yticks(ys); axs[0, 0].set_yticklabels([lab for _, lab in RULES])
    # r3 (实验管理 2026-11-25): every interval exists (d_d gives lo < hi for every cell); list those whose half-widths both lie inside the marker radius,
    #   so the caption can say that intervals narrower than the marker are hidden by it
    fig.canvas.draw(); hidden = []
    for ax, key, pt, lo, hi, y, ms in drawn:
        x0, x1, xp = (ax.transData.transform((x, y))[0] for x in (lo, hi, pt))
        if max(xp - x0, x1 - xp) <= ms * fig.dpi / 72 / 2: hidden.append(key)
    buf = P.render(fig); plt.close(fig)
    payload = json.dumps(dict(script_sha256=sha(__file__), style_sha256=P.STYLE_SHA256, master=str(a.master), master_sha256=MASTER_SHA, eps=eps, rules=[r for r, _ in RULES],
                              dsens_rules=sorted(DSENS_RULES), caption_note=CAPTION_NOTE, hidden_intervals=hidden, n_intervals=len(drawn),
                              mediabox_mm=buf["mediabox_mm"], shown=shown), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    for ext in ("pdf", "png"):
        with open(out / f"{STEM}.{ext}", "xb") as fh: fh.write(buf[ext])
    with open(out / f"{STEM}_numbers.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh:
        for nm in (f"{STEM}.pdf", f"{STEM}.png", f"{STEM}_numbers.json"): fh.write(f"{sha(out / nm)}  {nm}\n")
    print("written", out); return 0


if __name__ == "__main__": sys.exit(main())

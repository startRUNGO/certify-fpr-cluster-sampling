# P5a figure F-a (new, 实验管理 / user 2026-11-24 figure set): funnel plot of the field-level false-positive rate on BBBC-ALL.
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  Width 139 mm, IEEE style module p5a_style_ieee v2.0 (pinned), text >= 7 pt, pure black text.
# Data: only results_master keys d_a/* (d_a_funnel_bbbc.py 0ff2f656, VM run dafunnel_20261124T230844Z, ingested by add_dafunnel_keys.py f54f8bcd):
#   per field f = 0..n_fields-1: d_a/n_f|f (labeled negatives), d_a/X_f|f (negatives above the threshold tau_eps), d_a/frac_f|f = X_f / n_f;
#   funnel: d_a/funnel_n|j with the binomial 90% limits d_a/funnel_q05_frac|j and d_a/funnel_q95_frac|j around the site FPR d_a/fpr_site (exchangeable cells);
#   d_a/n_out_below, d_a/n_out_above, d_a/share_out, d_a/share_out_expected_binomial, d_a/dispersion_x2_over_df, d_a/icc_anova, d_a/design_effect go to the numbers json
#   (the caption reads them from the keys; no number is printed in the figure).
# Drawing: x = labeled negatives per field (log axis), y = field FPR; grey points inside the funnel, vermillion points outside; the funnel limits as two black lines
#   (step at each n), the site FPR as a dashed black line.
# Gates (exit 3, nothing written): master sha MASTER_SHA; style sha STYLE_SHA; d_a/eps == protocol/eps_main; n_f / X_f / frac_f each have exactly d_a/n_fields
#   entries; sum n_f == d_a/n_neg, sum X_f == d_a/k; |frac_f - X_f / n_f| <= 1e-12; |fpr_site - k / n_neg| <= 1e-12; funnel_n strictly increasing and covering every n_f;
#   q05 <= fpr_site <= q95 at every n; the counts of fields below q05 and above q95 recomputed from the keys equal d_a/n_out_below and d_a/n_out_above;
#   share_out == (below + above) / n_fields within 1e-12.
# Output (new directory, 'x' mode): fig_funnel_v65.pdf, .png (600 dpi), fig_funnel_v65_numbers.json, SHA256SUMS.
# Usage: python3 -E -B make_fig_funnel_v65.py --master <tree/results_master.json> --out <new dir>
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
STEM = "fig_funnel_v65"
TOL = 1e-12


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--master", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(); out = Path(a.out)
    if out.exists(): return fail("output directory exists", out)
    if P.STYLE_SHA256 != STYLE_SHA: return fail("style module sha differs from the pin")
    M, MSHA = PI.load_master(a.master, (MASTER_SHA,)); V = M["values"]; shown = {}

    def g(k):
        if k not in V or V[k] is None: raise KeyError(k)
        shown[k] = V[k]; return V[k]
    try:
        if g("d_a/eps") != g("protocol/eps_main"): return fail("d_a/eps differs from protocol/eps_main")
        F, nneg, k, fpr = int(g("d_a/n_fields")), g("d_a/n_neg"), g("d_a/k"), g("d_a/fpr_site")
        n = [g(f"d_a/n_f|{i}") for i in range(F)]; X = [g(f"d_a/X_f|{i}") for i in range(F)]; fr = [g(f"d_a/frac_f|{i}") for i in range(F)]
        J = sum(1 for kk in V if kk.startswith("d_a/funnel_n|"))
        fn = [g(f"d_a/funnel_n|{j}") for j in range(J)]; q05 = [g(f"d_a/funnel_q05_frac|{j}") for j in range(J)]; q95 = [g(f"d_a/funnel_q95_frac|{j}") for j in range(J)]
        nb, na, so = g("d_a/n_out_below"), g("d_a/n_out_above"), g("d_a/share_out")
        for kk in ("d_a/share_out_expected_binomial", "d_a/dispersion_x2_over_df", "d_a/icc_anova", "d_a/design_effect", "d_a/n_bar_w", "d_a/seed"): g(kk)
    except KeyError as e: return fail("missing key", e)
    if any(kk.startswith("d_a/n_f|") and int(kk.split("|")[1]) >= F for kk in V): return fail("more n_f entries than d_a/n_fields")
    if sum(n) != nneg or sum(X) != k: return fail("sums of n_f / X_f differ from d_a/n_neg / d_a/k")
    if any(abs(f_ - x_ / n_) > TOL for f_, x_, n_ in zip(fr, X, n)): return fail("frac_f differs from X_f / n_f")
    if abs(fpr - k / nneg) > TOL: return fail("fpr_site differs from k / n_neg")
    if any(b <= a_ for a_, b in zip(fn, fn[1:])) or not set(n) <= set(fn): return fail("funnel_n not increasing or not covering every n_f")
    if any(not (lo - TOL <= fpr <= hi + TOL) for lo, hi in zip(q05, q95)): return fail("site FPR outside the funnel at some n")
    lim = {m: (lo, hi) for m, lo, hi in zip(fn, q05, q95)}
    below = [i for i in range(F) if fr[i] < lim[n[i]][0]]; above = [i for i in range(F) if fr[i] > lim[n[i]][1]]
    if len(below) != nb or len(above) != na: return fail("fields outside the funnel differ from d_a/n_out_below / d_a/n_out_above", len(below), len(above))
    if abs(so - (nb + na) / F) > TOL: return fail("d_a/share_out differs from (below + above) / n_fields")
    P.apply()
    fig, ax = plt.subplots(figsize=P.figsize(70))
    out_set = set(below) | set(above)
    ins = [i for i in range(F) if i not in out_set]
    ax.plot([n[i] for i in ins], [fr[i] for i in ins], "o", ms=2.0, color=P.MUTED[0], alpha=0.7, mew=0, zorder=2, label="field inside the 90% limits")
    ax.plot([n[i] for i in sorted(out_set)], [fr[i] for i in sorted(out_set)], "o", ms=2.0, color=P.VERM, alpha=0.8, mew=0, zorder=3, label="field outside the 90% limits")
    ax.step(fn, q05, where="mid", color=P.BLACK, lw=P.LW_REF, zorder=4, label="binomial 90% limits")
    ax.step(fn, q95, where="mid", color=P.BLACK, lw=P.LW_REF, zorder=4)
    P.ref_line(ax, fpr, axis="y", color=P.BLACK)
    ax.plot([], [], color=P.BLACK, lw=P.LW_REF, ls=P.DASH_REF, label="site FPR")
    # r2 (render check 2026-11-25): top limit 1.03 so a field at FPR 1 is not cut by the frame
    ax.set_xscale("log"); ax.set_xlim(min(n) * 0.9, max(n) * 1.1); ax.set_ylim(-0.02, 1.03)
    ax.set_xlabel("Labeled negative cells in the field"); ax.set_ylabel("Field false-positive rate")
    ax.legend(loc="upper right", ncol=2)
    buf = P.render(fig); plt.close(fig)
    payload = json.dumps(dict(script_sha256=sha(__file__), style_sha256=P.STYLE_SHA256, master=str(a.master), master_sha256=MSHA,
                              recomputed=dict(n_out_below=len(below), n_out_above=len(above)), mediabox_mm=buf["mediabox_mm"],
                              shown={kk: v for kk, v in shown.items() if not kk.startswith(("d_a/n_f|", "d_a/X_f|", "d_a/frac_f|", "d_a/funnel_"))},
                              per_field_and_funnel_keys=dict(n_f=F, funnel=J)), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    for ext in ("pdf", "png"):
        with open(out / f"{STEM}.{ext}", "xb") as fh: fh.write(buf[ext])
    with open(out / f"{STEM}_numbers.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh:
        for nm in (f"{STEM}.pdf", f"{STEM}.png", f"{STEM}_numbers.json"): fh.write(f"{sha(out / nm)}  {nm}\n")
    print("written", out); return 0


if __name__ == "__main__": sys.exit(main())

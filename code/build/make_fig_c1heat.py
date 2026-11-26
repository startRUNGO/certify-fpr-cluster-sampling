# v1.3 (user 2026-11-24 via 实验管理): rule names in the figure text renamed (S3 -> DE-Grid, SPLIT6 -> Split-Select, Z -> Zero-Census); keys, data and layout unchanged.
# v1.2 (实验管理): cell values with |x| < 0.005 shown with 3 decimals (sign kept), others 2; minus sign in cells is U+2212.  v1.1 bd5032ee kept.
# P5a (CertBound): c1 heat map of the SPLIT6 measured lead over Z on the main grid, pi x CV panels by xi (columns) and eps (rows).  Default place: SI.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Reads ONLY results_master.json of the tree (sha pinned 58f900ee): keys c1card/result|cell_leads|<cell>|{measured, predicted_sign, small}|SPLIT6_single for the 90 main-grid cells, and c1card/result|P1|SPLIT6_single|{n_correct, n_in_denominator, n_small_excluded}.
# Cell names are parsed (main_pi<pi>_xi<xi>_cv<cv>_eps<eps>_m22) into the four factors; missing (pi, CV) combinations are DERIVED from the data (not hard-coded) and drawn light grey with hatching ("not in grid").
# Encoding: diverging colour map (RdBu) centred at 0, one symmetric range +-max|measured| for all panels, one colorbar; each cell shows the value with 2 decimals (font 6), text black or white by background luminance;
#   wrong-sign cells (P1 definition, small excluded) get a thick black frame; small (excluded) cells get a middle dot after the value.  All other text pure black.  Width 6.6 in.  No jitter, no 3-D.
# Gate (before any file is written): over all cells, right = n_correct, wrong = n_in_denominator - n_correct, small = n_small_excluded; any mismatch -> exit 3.  Existing outputs -> refuse.
# Output (figures/, 'x' semantics): fig_c1heat.pdf, fig_c1heat.png (600 dpi), fig_c1heat_numbers.json (every cell's value and key).
# Usage: python3 -B make_fig_c1heat.py --tree <tree_bj_v1_...>
import argparse, hashlib, json, re, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import cm, colors
MASTER_SHA = "58f900ee11eb4cd2e909e35a720c2ca9e2fdce7e04215f8bf0683d5f518cc56d"; RULE = "SPLIT6_single"
PAT = re.compile(r"^main_pi([0-9.]+)_xi(lo|mid|hi)_cv([0-9.]+)_eps([0-9.]+)_m22$")
XIS = ("lo", "mid", "hi"); EPSS = (0.05, 0.10)


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); a = ap.parse_args(); tree = Path(a.tree); rm = tree / "results_master.json"
    if sha(rm) != MASTER_SHA: print("GATE (exit 3): results_master sha differs", file=sys.stderr); return 3
    v = json.load(open(rm))["values"]; pre = "c1card/result|cell_leads|"
    cells = sorted({k[len(pre):].split("|")[0] for k in v if k.startswith(pre) and k[len(pre):].startswith("main_")})
    if len(cells) != 90: print("GATE (exit 3): main grid is not 90 cells", len(cells), file=sys.stderr); return 3
    outs = [tree / "figures" / n for n in ("fig_c1heat.pdf", "fig_c1heat.png", "fig_c1heat_numbers.json")]
    if any(o.exists() for o in outs): print("GATE (exit 3): output exists", file=sys.stderr); return 3
    D = {}; cnt = {"right": 0, "wrong": 0, "small": 0}; numbers = {}
    for c in cells:
        m = PAT.match(c)
        if not m: print("GATE (exit 3): cell name not parsed", c, file=sys.stderr); return 3
        pi, xi, cv, eps = float(m.group(1)), m.group(2), float(m.group(3)), float(m.group(4))
        meas = v[f"{pre}{c}|measured|{RULE}"]; sgn = bool(v[f"{pre}{c}|predicted_sign|{RULE}"]); sm = bool(v[f"{pre}{c}|small|{RULE}"])
        kind = "small" if sm else ("right" if (meas > 0) == sgn else "wrong"); cnt[kind] += 1
        D[(eps, xi, pi, cv)] = (meas, kind); numbers[c] = dict(pi=pi, xi=xi, cv=cv, eps=eps, measured=meas, kind=kind, key=f"{pre}{c}|measured|{RULE}")
    P1 = {k: v[f"c1card/result|P1|{RULE}|{k}"] for k in ("n_correct", "n_in_denominator", "n_small_excluded")}
    if not (cnt["right"] == P1["n_correct"] and cnt["wrong"] == P1["n_in_denominator"] - P1["n_correct"] and cnt["small"] == P1["n_small_excluded"]):
        print("GATE (exit 3): counts differ from P1", cnt, P1, file=sys.stderr); return 3
    PIS = sorted({k[2] for k in D}); CVS = sorted({k[3] for k in D})
    if set(k[0] for k in D) != set(EPSS) or set(k[1] for k in D) != set(XIS): print("GATE (exit 3): eps / xi levels differ", file=sys.stderr); return 3
    vmax = max(abs(x[0]) for x in D.values()); norm = colors.Normalize(-vmax, vmax); cmap = matplotlib.colormaps["RdBu"]
    plt.rcParams.update({"font.family": "serif", "font.size": 7, "text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black", "axes.linewidth": 0.5, "hatch.linewidth": 0.4})
    fig, axes = plt.subplots(2, 3, figsize=(6.6, 4.6), constrained_layout=True, sharex=True, sharey=True); missing_drawn = False
    for r, eps in enumerate(EPSS):
        for col, xi in enumerate(XIS):
            ax = axes[r, col]
            for i, pi in enumerate(PIS):
                for j, cv in enumerate(CVS):
                    key = (eps, xi, pi, cv)
                    if key not in D:
                        ax.add_patch(Rectangle((i, j), 1, 1, facecolor="#e6e6e6", edgecolor="#a0a0a0", hatch="///", linewidth=0.5)); missing_drawn = True; continue
                    val, kind = D[key]; rgba = cmap(norm(val)); lum = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
                    ax.add_patch(Rectangle((i, j), 1, 1, facecolor=rgba, edgecolor="white", linewidth=0.5))
                    if kind == "wrong": ax.add_patch(Rectangle((i + 0.04, j + 0.04), 0.92, 0.92, facecolor="none", edgecolor="black", linewidth=1.4))
                    ax.text(i + 0.5, j + 0.5, (f"{val:.3f}" if abs(val) < 0.005 else f"{val:.2f}").replace("-", "\u2212") + ("·" if kind == "small" else ""), ha="center", va="center", fontsize=6, color=("white" if lum < 0.45 else "black"))
            ax.set_xlim(0, len(PIS)); ax.set_ylim(0, len(CVS)); ax.set_xticks([i + 0.5 for i in range(len(PIS))]); ax.set_xticklabels([f"{p:g}" for p in PIS])
            ax.set_yticks([j + 0.5 for j in range(len(CVS))]); ax.set_yticklabels([f"{c:g}" for c in CVS]); ax.set_aspect("equal"); ax.tick_params(length=0)
            ax.set_title(f"$\\xi$ = {xi}, $\\varepsilon$ = {eps:g}", fontsize=7, color="black")
            if r == 1: ax.set_xlabel("Prevalence $\\pi$")
            if col == 0: ax.set_ylabel("Field CV")
    sm_ = cm.ScalarMappable(norm=norm, cmap=cmap); cb = fig.colorbar(sm_, ax=axes, shrink=0.8, pad=0.02); cb.set_label("Measured lead of Split-Select over Zero-Census", color="black"); cb.ax.tick_params(labelsize=6, colors="black")
    handles = [Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black", linewidth=1.4, label="sign predicted wrongly"), Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="white", label="· |predicted| < 0.01 (excluded)")]
    if missing_drawn: handles.append(Rectangle((0, 0), 1, 1, facecolor="#e6e6e6", hatch="///", edgecolor="#a0a0a0", label="not in grid"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.45, -0.05), fontsize=6)
    fig.savefig(outs[0], bbox_inches="tight"); fig.savefig(outs[1], dpi=600, bbox_inches="tight")
    with open(outs[2], "x") as fh: json.dump(dict(master_sha256=MASTER_SHA, script_sha256=sha(__file__), vmax=vmax, counts=cnt, P1=P1, pi_levels=PIS, cv_levels=CVS, cells=numbers), fh, indent=1)
    print({o.name: sha(o) for o in outs}, cnt); return 0


if __name__ == "__main__": sys.exit(main())

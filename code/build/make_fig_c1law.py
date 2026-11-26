# v1.3 (user 2026-11-24 via 实验管理): rule names in the figure text renamed (S3 -> DE-Grid, SPLIT6 -> Split-Select, Z -> Zero-Census); keys, data and layout unchanged.
# P5a (CertBound): c1-law figure (predicted vs measured lead over the zero-label rule Z), main grid, S3_single and SPLIT6_single.  v1.1 (实验管理): count gate against P1 (right = n_correct, right + wrong = n_in_denominator, small = n_small_excluded) before anything is drawn or written; legend = union of both panels' handles, de-duplicated by label.  v1 4747fb6f kept.  v1.2 (实验管理): one common axis range for both panels (union of both panels' data, 5% pad) with sharex / sharey; every marker gets a 0.3 pt white edge, alpha 0.85, triangles s = 12 (circles stay 14); no jitter (plotted values unchanged).  v1.1 0d5687a5 kept, its outputs renamed NOT_USED_v1.1_.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Reads ONLY results_master.json of the tree (keys c1card/result|cell_leads|<cell>|{predicted, measured, predicted_sign, small}|<rule> and c1card/result|P1|<rule>|{n_correct, n_in_denominator}); master sha pinned.
# Marks: correct sign (filled circle), wrong sign (filled triangle), |predicted| < 0.01 excluded from P1 (open circle, grey).  Text pure black; colour only on marks; reference lines y = x (dashed) and the zero lines (thin).
# Output (into <tree>/figures, 'x' mode): fig_c1law.pdf, fig_c1law.png (600 dpi), fig_c1law_numbers.json (the plotted numbers with their keys).
# Usage: python3 -B make_fig_c1law.py --tree <tree_bj_v1_...>
import argparse, hashlib, json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
MASTER_SHA = "58f900ee11eb4cd2e909e35a720c2ca9e2fdce7e04215f8bf0683d5f518cc56d"
RULES = (("S3_single", "DE-Grid"), ("SPLIT6_single", "Split-Select"))
COL = {"right": "#1f5fa8", "wrong": "#d9711c", "small": "#9a9a9a"}


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); a = ap.parse_args(); tree = Path(a.tree); rm = tree / "results_master.json"
    if sha(rm) != MASTER_SHA: print("GATE (exit 3): results_master sha differs", file=sys.stderr); return 3
    v = json.load(open(rm))["values"]; pre = "c1card/result|cell_leads|"
    cells = sorted({k[len(pre):].split("|")[0] for k in v if k.startswith(pre) and k[len(pre):].startswith("main_")})
    if len(cells) != 90: print("GATE (exit 3): main grid is not 90 cells", len(cells), file=sys.stderr); return 3
    outs = [tree / "figures" / n for n in ("fig_c1law.pdf", "fig_c1law.png", "fig_c1law_numbers.json")]
    if any(o.exists() for o in outs): print("GATE (exit 3): output exists", file=sys.stderr); return 3
    plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.linewidth": 0.6, "text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black"})
    PTS = {}; numbers = {}
    for rule, label in RULES:   # gate BEFORE any drawing or writing
        pts = []
        for c in cells:
            g = lambda f: v[f"{pre}{c}|{f}|{rule}"]
            p, m, s, sm = g("predicted"), g("measured"), g("predicted_sign"), g("small")
            kind = "small" if sm else ("right" if (m > 0) == bool(s) else "wrong"); pts.append((p, m, kind)); numbers[f"{rule}|{c}"] = dict(predicted=p, measured=m, kind=kind)
        n = {kk: sum(1 for _, _, x in pts if x == kk) for kk in ("right", "wrong", "small")}
        P1 = {kk: v[f"c1card/result|P1|{rule}|{kk}"] for kk in ("n_correct", "n_in_denominator", "n_small_excluded")}
        if not (n["right"] == P1["n_correct"] and n["right"] + n["wrong"] == P1["n_in_denominator"] and n["small"] == P1["n_small_excluded"]):
            print("GATE (exit 3): counts differ from P1", rule, n, P1, file=sys.stderr); return 3
        PTS[rule] = pts
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.1), constrained_layout=True, sharex=True, sharey=True)
    allv = [z for r in PTS.values() for x, y, _ in r for z in (x, y)]; LO, HI = min(allv), max(allv); PAD = 0.05 * (HI - LO); LO -= PAD; HI += PAD
    for ax, (rule, label), tag in zip(axes, RULES, ("(a)", "(b)")):
        pts = PTS[rule]
        lo, hi = LO, HI
        ax.plot([lo, hi], [lo, hi], ls="--", lw=0.6, color="black", zorder=1); ax.axhline(0, lw=0.4, color="black", zorder=1); ax.axvline(0, lw=0.4, color="black", zorder=1)
        for kind, mk, fc in (("right", "o", COL["right"]), ("wrong", "^", COL["wrong"]), ("small", "o", "none")):
            xs = [x for x, _, k in pts if k == kind]; ys = [y for _, y, k in pts if k == kind]
            ax.scatter(xs, ys, s=(12 if mk == "^" else 14), marker=mk, facecolors=fc, edgecolors=(COL["small"] if kind == "small" else "white"), linewidths=(0.6 if kind == "small" else 0.3), alpha=0.85, zorder=3,
                       label={"right": "sign predicted correctly", "wrong": "sign predicted wrongly", "small": "|predicted| < 0.01 (excluded)"}[kind])
        nc = v[f"c1card/result|P1|{rule}|n_correct"]; nd = v[f"c1card/result|P1|{rule}|n_in_denominator"]
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
        ax.set_xlabel(f"Predicted lead of {label} over Zero-Census"); ax.set_ylabel(f"Measured lead of {label} over Zero-Census")
        ax.set_title(f"{tag} {label}: sign correct in {nc} of {nd} cells", fontsize=8, loc="left", color="black")
        numbers[f"P1|{rule}"] = dict(n_correct=nc, n_in_denominator=nd)
    seen = {}
    for ax in axes:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in seen: seen[ll] = hh
    fig.legend(list(seen.values()), list(seen.keys()), loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.08), fontsize=7)
    fig.savefig(outs[0], bbox_inches="tight"); fig.savefig(outs[1], dpi=600, bbox_inches="tight")
    with open(outs[2], "x") as fh: json.dump(dict(master_sha256=MASTER_SHA, script_sha256=sha(__file__), numbers=numbers), fh, indent=1)
    print({o.name: sha(o) for o in outs}); return 0


if __name__ == "__main__": sys.exit(main())

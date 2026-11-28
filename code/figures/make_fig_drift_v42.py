# P5a v3.9 Figure 4 (fig:drift): FPR at the source threshold against AUROC for the 35 surveyed checkpoints on five sites.  Base: make_figures_v2.py (5c066c31) Fig. 2 block (fig2_fpr_vs_auroc).
# Design r2 (ca3bf65e) Figure 4 section; width 139 mm, text >= 7 pt, pure black text.
# Data (site versions as in the main text, checked against the base, which already used the rebuilt Lacuna-B and MP-IDB rows):
#   (a) NIH held-out fold and (b) BBBC-ALL: rows of results_src/s6/s6a_survey.json (sha 9d194a56, pinned; these rows are not master keys), site nih_locked_test / bbbc041_all;
#   (c) MP-IDB: master family mpidb_hc_b2__mh_s6a_survey/<i>|{arm, auroc, fpr} (rebuilt MP-IDB; NOT the old mpidb_full rows of s6a_survey.json, rho -0.405);
#   (d) Lacuna-B: e7_a_survey/rows|<i>|{arm, auroc, fpr} (source lh_survey.json, the rebuilt Lacuna-B); (e) PlasmoCount: e8_plasmo_survey/rows|<i>|{arm, auroc, fpr}.
#   rho_s shown from keys (not recomputed for display): s6a_spearman/nih_locked_test|rho, s6a_spearman/bbbc041_all|rho, mpidb_hc_b2__mh_s6a_spearman/mh|rho,
#   lacuna_hc__lh_s6a_spearman/lh|rho, e8_plasmo_survey/summary|spearman_rho.  protocol/eps_main, protocol/seeds, facts/ckpt_survey_n.
# Gates (exit 3): file sha; 35 = facts/ckpt_survey_n rows per site; 7 arms x protocol/seeds rows each; rho recomputed from the plotted rows equals the shown key to 1e-12;
#   for NIH and BBBC the a2 / a0-ctrl FPR and AUROC means equal survey_arm/<site>|<arm>|{fpr_mean, auroc_mean} to 1e-12 (ties the pinned file to the master).
# Arm names (SI Table tab:arms; codes checked against _opt_20261021/FACTS.md lines 19-25): a2 normalized checkpoint, a0-ctrl raw-front-end control, b0 parent, c0 identity pseudo-site control,
#   c1 pseudo-site, no front end, a3 normalized, frozen gating, a4 normalized, frozen gating, GRL.
# Changes from the base (r2): arm names instead of codes; legend title and target value from keys; normalized checkpoint black filled diamond, raw-front-end control dark grey filled circle,
#   the other five arms small light-grey open markers (vermillion and blue are reserved for limits and SPLIT6); panel (e) "PlasmoCount (confirmation site)"; sentence-case titles.
# Output (new directory, 'x' mode): fig4_drift_v39.pdf, .png (300 dpi), fig4_drift_v39_numbers.json, SHA256SUMS.
# Usage: python3 -B make_fig_drift_v39.py --master <results_master.json> --s6a <results_src/s6/s6a_survey.json> --out <new dir>
# public-repo patch (2026-11-28): the master is read through public_inputs.load_master, which accepts the public redacted results/results_master.json.gz
#   (uncompressed sha 4b296884) or the pinned internal master; non-public inputs stop with a plain message.  Drawing code unchanged.
import argparse, hashlib, io, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# v4.2 (draft 23 O3a, 2026-11-26): plain legend labels (deployed classifier, raw-image control, five training variants) and legend title; the caption defines them, Table S1.1 keeps the full key.  Data, keys and gates unchanged.
# v4.1 (draft 21 m8, 2026-11-26): the Spearman symbol in the x labels is rho_S, as in the text.  Data, keys and gates unchanged.
# v2.0 (user: IEEE style; module p5a_style_ieee v2.0): Times, (a)-(e) labels, constrained layout; the normalized checkpoint (the deployed classifier) HERO
#   (black diamonds, white edge, top layer), the raw-front-end control and the five other arms MUTED; FPR above the target eps as a light vermillion band, at or below it a light
#   blue band; one number per panel: rho_s in the x label; "(5 seeds each)" and "eps = 0.10" leave the legend (caption).  Data, keys and gates unchanged.
import p5a_style_ieee as P
import public_inputs as PI
S6A_SHA = "9d194a56738dbdf7626f87f6f3b026b67ab3912b8cf4aaa4aca6de65d9a3a1c2"
W, H = P.figsize(100)
BLACK, G1, G2 = P.BLACK, P.MUTED[6], P.MUTED[0]
ARMS = ["a2", "a0-ctrl", "b0", "c0", "c1", "a3", "a4"]
NAME = {"a2": "Deployed classifier", "a0-ctrl": "Raw-image control", "b0": "Original training", "c0": "Empty augmentation", "c1": "Augmentation only",
        "a3": "Random frozen gating", "a4": "Adversarial frozen\ngating"}
STY = {"a2": dict(marker="D", s=26, facecolor=BLACK, edgecolor="white", linewidth=0.5, zorder=6), "a0-ctrl": dict(marker="o", s=14, facecolor=G1, edgecolor=G1, linewidth=0.3, zorder=4)}
for a, mk in zip(["b0", "c0", "c1", "a3", "a4"], ["^", "v", "s", "<", ">"]): STY[a] = dict(marker=mk, s=10, facecolor="none", edgecolor=G2, linewidth=0.6, zorder=3)
SITES = [("nih_locked_test", "NIH held-out fold (source)"), ("bbbc041_all", "BBBC-ALL"), ("mpidb_full", "MP-IDB"), ("lacuna_B", "Lacuna-B"), ("plasmo", "PlasmoCount (confirmation site)")]
RHO = {"nih_locked_test": "s6a_spearman/nih_locked_test|rho", "bbbc041_all": "s6a_spearman/bbbc041_all|rho", "mpidb_full": "mpidb_hc_b2__mh_s6a_spearman/mh|rho",
       "lacuna_B": "lacuna_hc__lh_s6a_spearman/lh|rho", "plasmo": "e8_plasmo_survey/summary|spearman_rho"}


def sha(p):
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def fam(V, pre):
    idx = sorted({int(k[len(pre):].split("|")[0]) for k in V if k.startswith(pre)})
    return [(V[f"{pre}{i}|arm"], V[f"{pre}{i}|auroc"], V[f"{pre}{i}|fpr"]) for i in idx]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--master", required=True); ap.add_argument("--s6a", default=None); ap.add_argument("--out", required=True); a = ap.parse_args(); out = Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    PI.need_file(a.s6a, S6A_SHA, "s6a_survey.json (per-checkpoint survey rows of panels a and b)")
    M, MSHA = PI.load_master(a.master); V = M["values"]; eps = V["protocol/eps_main"]; seeds = V["protocol/seeds"]; nck = V["facts/ckpt_survey_n"]
    rows = json.load(open(a.s6a))
    data = {s: [(r["arm"], r["auroc"], r["fpr"]) for r in rows if r["site"] == s] for s in ("nih_locked_test", "bbbc041_all")}
    data["mpidb_full"] = fam(V, "mpidb_hc_b2__mh_s6a_survey/"); data["lacuna_B"] = fam(V, "e7_a_survey/rows|"); data["plasmo"] = fam(V, "e8_plasmo_survey/rows|")
    shown = {}
    for s, _ in SITES:
        pts = data[s]
        if len(pts) != nck: return fail("rows", s, len(pts))
        if sorted(set(p[0] for p in pts)) != sorted(ARMS) or any(sum(p[0] == arm for p in pts) != seeds for arm in ARMS): return fail("arms x seeds", s)
        rho = V[RHO[s]]; r2 = spearmanr([p[1] for p in pts], [p[2] for p in pts]).statistic
        if not abs(r2 - rho) <= 1e-12: return fail("rho recomputed differs from key", s, r2, rho)
        shown[s] = dict(rho=rho, rho_key=RHO[s], rho_shown=format(rho, "+.2f"))
    for s in ("nih_locked_test", "bbbc041_all"):
        for arm in ("a2", "a0-ctrl"):
            for j, nm in ((2, "fpr_mean"), (1, "auroc_mean")):
                k = f"survey_arm/{s}|{arm}|{nm}"
                if k in V and not abs(np.mean([p[j] for p in data[s] if p[0] == arm]) - V[k]) <= 1e-12: return fail("survey_arm mean differs", k)
                if k not in V: return fail("missing key", k)
    P.apply()
    fig, axs = plt.subplots(2, 3, figsize=(W, H), sharey=True)
    for k, (ax, (s, title)) in enumerate(zip(axs.flat, SITES)):
        P.band(ax, eps, 1.0, color=P.VERM, alpha=0.08); P.band(ax, -0.02, eps, color=P.SKY, alpha=0.16)
        for arm in ARMS[::-1]:
            q = [(au, f) for (ar, au, f) in data[s] if ar == arm]
            ax.scatter([p[0] for p in q], [p[1] for p in q], **STY[arm])
        P.ref_line(ax, eps, axis="y")
        P.panel_label(ax, "abcde"[k], title)
        ax.set_xlabel("AUROC (" + f"$\\rho_{{\\mathrm{{S}}}}$ = {shown[s]['rho_shown']})".replace("-", "−")); ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    for r in range(2): axs[r, 0].set_ylabel("FPR at source threshold")
    axs[0, 0].set_ylim(-0.02, 0.96); axs[0, 0].set_yticks([0, 0.1, 0.2, 0.4, 0.6, 0.8])
    leg = axs[1, 2]; leg.axis("off")
    hs = [Line2D([], [], ls="", marker=STY[x]["marker"], markerfacecolor=STY[x]["facecolor"], markeredgecolor=(STY[x]["edgecolor"] if STY[x]["facecolor"] == "none" else STY[x]["facecolor"]),
                 ms=np.sqrt(STY[x]["s"]) * 1.1, markeredgewidth=STY[x]["linewidth"]) for x in ARMS] + [Line2D([], [], color=BLACK, lw=P.LW_REF, ls=P.DASH_REF)]
    leg.legend(hs, [NAME[x] for x in ARMS] + [r"target $\epsilon$"], loc="center left", title="Classifier", title_fontsize=P.FS_LEGEND, alignment="left",
               borderaxespad=0.0, handletextpad=0.4, handlelength=1.0)
    buf = P.render(fig)
    plt.close(fig)
    payload = json.dumps(dict(script_sha256=sha(__file__), style_sha256=P.STYLE_SHA256, master=str(a.master), master_sha256=MSHA, s6a_sha256=S6A_SHA, eps=eps, seeds=seeds, ckpt_survey_n=nck, rho=shown,
                              rows={s: [list(p) for p in data[s]] for s, _ in SITES}), indent=1, allow_nan=False)
    out.mkdir(parents=True)
    for ext in ("pdf", "png"):
        with open(out / f"fig4_drift_v39.{ext}", "xb") as fh: fh.write(buf[ext])
    with open(out / "fig4_drift_v39_numbers.json", "x") as fh: fh.write(payload)
    with open(out / "SHA256SUMS", "x") as fh:
        for nm in ("fig4_drift_v39.pdf", "fig4_drift_v39.png", "fig4_drift_v39_numbers.json"): fh.write(f"{sha(out / nm)}  {nm}\n")
    print("written", out); return 0


if __name__ == "__main__": sys.exit(main())

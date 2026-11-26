# P5a BJ v4.0 (assemble v1.24): re-render the ten IEEE figures into a new directory and check them against the reviewed renders (实验管理, 2026-11-24).
# WRITTEN, NOT RUN, until 实验管理 reviews it.
# Gates (exit 3): every figure script, the style module p5a_style_ieee.py, the master and s6a_survey.json have their pinned sha before anything runs; each script exits 0;
#   each new PNG is byte-identical to the reviewed render (sha pinned below).  PDFs are not compared byte for byte (matplotlib writes the creation date into the PDF);
#   their shas are recorded in <out>/FIGS_SUMS, which assemble_bj_v4.py checks.
# Writes only under --out (new, not under CloudStorage); each figure script writes its own new subdirectory in 'x' mode; FIGS_SUMS in 'x' mode.
# v1.2 (实验管理): Fig 7 = make_fig_dev_v39 v2.1 (148 mm), script e2a1070a, reviewed PNG 5009a607.
# v1.3 (实验管理 2026-11-24, rule rename): seven scripts re-pinned with the renamed figure text (overview v2.1, sim v2.1, pooled v2.3, dev v2.2, census v2.2, esparse v2.2, confirm v2.1); reviewed PNGs
#   from the approved trial renders (rename_trial_20261124T202515Z, census_v22_202813); minfields, mechanism and drift unchanged.
# Usage: python3 -B render_figs_v40.py --p5a <P5a_CertBound> --stainref <Paper5_StainRef> --out <new dir>
import argparse, hashlib, subprocess, sys
from pathlib import Path
STYLE = ("p5a_style_ieee.py", "8424218990b9bec4dcb8b3bf6c54d7b2da5671b2fc9444511ccbefa7edd39ce3")
MASTER = ("tree_bj_v1_20261123T191030Z/results_master.json", "1f68e729ba8c2f499021363ea2bc59f6d6776bc64483559d0a0a5f37dc6e13bc")
S6A = ("results_src/s6/s6a_survey.json", "9d194a56738dbdf7626f87f6f3b026b67ab3912b8cf4aaa4aca6de65d9a3a1c2")
FIGDATA = "tree_bj_v1_20261123T191030Z/results_src/figmech_v39"
# subdir, script, script sha (= script_sha256 of the reviewed render), extra args kind, output PNG, reviewed PNG sha (render dir in the scratchpad)
FIGS = [("fig01_overview", "make_fig_overview_v39.py", "5725e62337312003d46aee79207e4a19e7aa65bdb7be89c534b1acfa89c65c37", "none", "fig1_overview_v39",
         "177c4aa3874e94eba8a8dacf6fcfd1ff2f9eaa3a907e2df5914700fbb1d8b0db"),
        ("fig02_minfields", "make_fig_minfields_v39.py", "5f42bac4fd9a175f4c6ad7f75527664ff6624cedc8d8ed3bb9719182aaa9559c", "master", "fig3_minfields_v39",
         "61ce5d92384fed1865fa42891c0eed86172553d910bcddc986a2c21d241018e0"),
        ("fig03_mechanism", "make_fig_mechanism_v39.py", "228d0789f618aca9e04d6836949e1b4b62a756fb314f48cc575c0800a52c0f8d", "mech", "fig2_mechanism_v39",
         "ec09152dc6e007be8842eb30d45c4a293db9d5b38f1d879383e32af41f448ca7"),
        ("fig04_sim", "make_fig_sim_v39.py", "fe4a37c42defbee918019a1efdc05de676759d45037813c0a546212cb0e06343", "master", "fig3_sim_v39",
         "2ddb48a01eafed7182385852eafe5539c1053031560a95f9ae67f5e79afdbaf5"),
        ("fig05_drift", "make_fig_drift_v39.py", "22a69c56599e1c4740097f8625752b59b8bd0b22d41a06bb255e63cf42a5afc1", "drift", "fig4_drift_v39",
         "5382080c22b7ad7eea7d0b72bd088461451897968349e0cef4cf75cea76f85e7"),
        ("fig06_pooled", "make_fig_pooled_v39.py", "faf119617ae6a7aa7e552e5e3f953fa4fe17546a7d88cc4ce00d8c3e09174957", "master", "fig7_pooled_v39",
         "2e42f591e178f683a221d8a528e447cdc8d0a12ad906e0671d86dd0b139a737f"),
        ("fig07_dev", "make_fig_dev_v39.py", "060c1ef3227283636d9e653d545dff3a2d75c1f29a89f01c33003311de1543c8", "master", "fig5_dev_v39",
         "4d3645ee4b292943123d803d50cccdfceb6519791bf9d730948c499f18126970"),
        ("fig08_census", "make_fig_census_v39.py", "9aada9a1b554be902c3745b1bfc4db467e7e9b80fa1b560bcaab040eb1984f21", "master", "fig8_census_v39",
         "3f96269369f84eaab9f9dada0e77f4e3fc2b1616f977914b7fd67343ee69a1c6"),
        ("fig09_esparse", "make_fig_esparse_v39.py", "e62a84970d34ffffc9252e23c9a5f78cf5bfae72e2c1ad2246a193ebdf919cea", "master", "fig9_esparse_v39",
         "030ccc165333101b706f0ff72ccbb59e372749288a05a2507b666d5592fad30e"),
        ("fig10_confirm", "make_fig_confirm_v39.py", "1af2cd2e36bc0523858817b2739941c94e90c91f109f61421820626e9e186097", "master", "fig10_confirm_v39",
         "2efa68e9d5936b8d6fb10d292924a351df9a742c2eb3abd5c620beabda8aa1fa")]


def sha(p):
    # missing file -> None, so every pin / output check below fails through fail() instead of a traceback (Gemini review)
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--p5a", required=True); ap.add_argument("--stainref", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    P, R, O = Path(a.p5a).resolve(), Path(a.stainref).resolve(), Path(a.out).resolve(); F = P / "figs_v39"
    if O.exists(): return fail("output directory exists", O)
    # v1.1 (实验管理): real path test instead of a substring test (the session scratchpad name contains the word "CloudStorage")
    if (Path.home() / "Library" / "CloudStorage").resolve() in O.resolve().parents: return fail("output directory under CloudStorage", O)
    for base, (rel, h) in ((F, STYLE), (P, MASTER), (R, S6A)):
        if sha(base / rel) != h: return fail("pin", rel)
    for sub, scr, h, *_ in FIGS:
        if sha(F / scr) != h: return fail("script pin", scr)
    M = str(P / MASTER[0]); O.mkdir(parents=True); sums = []
    for sub, scr, h, kind, stem, png_sha in FIGS:
        args = {"none": [], "master": ["--master", M], "mech": ["--data", str(P / FIGDATA), "--master", M],
                "drift": ["--master", M, "--s6a", str(R / S6A[0])]}[kind]
        r = subprocess.run([sys.executable, "-B", str(F / scr)] + args + ["--out", str(O / sub)], cwd=str(F), capture_output=True, text=True)
        if r.returncode != 0: return fail(scr, "exit", r.returncode, r.stderr[-400:])
        got = sha(O / sub / f"{stem}.png")
        if got != png_sha: return fail(scr, "PNG differs from the reviewed render", got)
        pdf = sha(O / sub / f"{stem}.pdf")
        if pdf is None: return fail(scr, "no PDF written")
        sums += [f"{pdf}  {sub}/{stem}.pdf", f"{got}  {sub}/{stem}.png"]
        print("ok", sub, "PNG identical to the reviewed render")
    with open(O / "FIGS_SUMS", "x") as fh: fh.write("\n".join(sums) + "\n")
    print("written", O / "FIGS_SUMS"); return 0


if __name__ == "__main__": sys.exit(main())

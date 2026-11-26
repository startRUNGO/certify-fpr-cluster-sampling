# P5a v5.2: re-render the three SI figures with the renamed rule names (fig_s_zpi, fig_c1law, fig_c1heat) into a new directory and check them against the reviewed renders.
# WRITTEN, NOT RUN, until 实验管理 reviews it.  Writes only under --out (new, not under CloudStorage).
# Gates (exit 3): the three scripts and their two inputs have their pinned sha before anything runs; each script exits 0; each PNG is byte-identical to the reviewed trial render
#   (PDFs carry a creation date and are not compared byte for byte; their shas go to <out>/SIFIGS_SUMS, which assemble_bj_v4.py checks).
# The c1 scripts read <tree>/results_master.json pinned 58f900ee and write <tree>/figures/ ('x' semantics), so a scratch tree is built under --out with a copy of the pinned master.
# Usage: python3 -B render_sifigs_v52.py --stainref <Paper5_StainRef> --out <new dir>
import argparse, hashlib, shutil, subprocess, sys
from pathlib import Path

SCRIPTS = {"zpi": ("handover_20261101/fig_zpi_20261116T1241Z/make_zpi.py", "cb2b1dd7882ca711bab834b8da2e882b20b0dbbd03824f03438a4b7eb22d6a61"),
           "c1law": ("P5a_CertBound/build/make_fig_c1law.py", "cbca41e272cadc3c1cf748231c1ad1b9eec9378a59df886bc2a1bda44ea96185"),
           "c1heat": ("P5a_CertBound/build/make_fig_c1heat.py", "57a1d7351391555acb0f22d7a7704cafdee898448fc84449222fed3fac5b69e0")}
ZLIMIT = ("results_src/z0/out_zlimit_20261116T1232Z/zlimit.json", "1cb1910506eb3a267ca88eed598716e46eb46dc5031aa8a76439dc035cfcf0ab")
MASTER = ("P5a_CertBound/tree_bj_v1_20261123T191030Z/results_master_before_p5a2_58f900ee.json", "58f900ee11eb4cd2e909e35a720c2ca9e2fdce7e04215f8bf0683d5f518cc56d")
# reviewed trial renders (scratchpad sifig_trial_205258: zpi6/ (v2.1, legend in panel a), tree/figures/)
PNG = {"fig_s_zpi": "7f58957b1e364d78a3fb5aaa285f87684c6f2848f22e59a3a3465180e187879f", "fig_c1law": "a4e83800e2140a734bd82ecd27f2af1d4b7185474b330e0f12f78e3ddba43008",
       "fig_c1heat": "c279201f94600cbe536259ca21400447028b0fc265ef16337d51a6b40fdd718e"}


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def fail(*m): print("GATE (exit 3):", *m, file=sys.stderr); return 3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--stainref", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    R, O = Path(a.stainref).resolve(), Path(a.out).resolve()
    if O.exists(): return fail("output directory exists", O)
    if (Path.home() / "Library" / "CloudStorage").resolve() in O.parents: return fail("output directory under CloudStorage", O)
    for k, (rel, h) in list(SCRIPTS.items()) + [("zlimit", ZLIMIT), ("master", MASTER)]:
        if sha(R / rel) != h: return fail("pin", k, rel)
    O.mkdir(parents=True); tree = O / "tree"; (tree / "figures").mkdir(parents=True); shutil.copyfile(R / MASTER[0], tree / "results_master.json")
    if sha(tree / "results_master.json") != MASTER[1]: return fail("master copy differs")
    runs = [("zpi", [str(R / ZLIMIT[0]), str(O / "zpi")], O / "zpi"), ("c1law", ["--tree", str(tree)], tree / "figures"), ("c1heat", ["--tree", str(tree)], tree / "figures")]
    sums = []
    for k, args, d in runs:
        r = subprocess.run([sys.executable, "-B", str(R / SCRIPTS[k][0])] + args, capture_output=True, text=True)
        if r.returncode != 0: return fail(k, "exit", r.returncode, r.stderr[-400:])
    for stem, d in (("fig_s_zpi", O / "zpi"), ("fig_c1law", tree / "figures"), ("fig_c1heat", tree / "figures")):
        if sha(d / f"{stem}.png") != PNG[stem]: return fail(stem, "PNG differs from the reviewed render", sha(d / f"{stem}.png"))
        if sha(d / f"{stem}.pdf") is None: return fail(stem, "no PDF")
        sums += [f"{sha(d / f'{stem}.pdf')}  {d.relative_to(O)}/{stem}.pdf", f"{PNG[stem]}  {d.relative_to(O)}/{stem}.png"]
        print("ok", stem, "PNG identical to the reviewed render")
    with open(O / "SIFIGS_SUMS", "x") as fh: fh.write("\n".join(sums) + "\n")
    print("written", O / "SIFIGS_SUMS"); return 0


if __name__ == "__main__": sys.exit(main())

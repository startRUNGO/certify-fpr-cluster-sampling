"""Generates the C5b evaluation scripts (PREREG_P5_C5B_FRAMING_AUG_20261026.md + errata; main-session spec of 2026-10-26) from the REVIEWED scripts (sha256 asserted below):
  c5b_bbbc_rev1.py  <- rev1.py (original, BBBC-ALL only: SITES = ["bbbc041_all"], so key8 = 8 and sidx = 0 exactly as the original stream for that site)
  c5b_mh_rev1.py    <- mh_rev1.py            (site index 31, unchanged)
  c5b_lh_rev1.py    <- lh_paired_rev1.py     (paired stream 8 / 2, unchanged; the OLD/NEW negative-score site check is removed: the models differ, so scores differ by design)
  c5b_lbsame_rev1.py<- lb_same_paired_rev1.py (paired stream 8 / 2, unchanged)
  c5b_bbbc_s11.py   <- _opt_20261021/s11.py (BBBC-ALL only; seed [..., 8, ai, sd, m, si] with si = 0 as in the original)
  c5b_mh_s11.py     <- mh_s11.py
  c5b_lh_s11.py     <- lh_paired_s11.py
Edits (each must hit exactly once): (1) ARMS = ["a2"] so that the arm index ai stays 0, i.e. every random stream is the frozen A2's (pairing); the name "a2" is only the key used inside the results, the files carry the arm label LAB (env C5B_ARM: c5b-f / c5b-fs; "a2" is allowed for the smoke only);
(2) the task grid loops over one arm; (3) the logits are read from C5B_LOGDIR/<LAB>_s<seed>_<site>_logits.npz (key logits) instead of the A2 files; the ids stay the ones the original script read; (4) the result files carry LAB in their name; the per-draw file is <script stem>_perdraw_<LAB>_m22.npz;
(5) BBBC scripts and mh_rev1 / mh_s11 gain the per-draw output exactly as build_paired.py added it to the lh / lb_dedup scripts; (6) a header line. Run from this directory (results_src/c5b). Existing files are never overwritten."""
import difflib, hashlib
from pathlib import Path
RS = Path("..")
SRC = {"bbbc_rev1": (RS / "rev1/rev1.py", "7bcb1635464757458968be79c2a33cd4bd2b1c77e7a4d31bd9fa870d49fc3ae0"), "mh_rev1": (RS / "mh/mh_rev1.py", "47f79f1df632e5efd1be5ea9485826d082eeaac38f89ee6f60a6aed2292116c5"),
       "lh_rev1": (RS / "paired/lh_paired_rev1.py", "9768296ca546891d24db3631f12a339b4f3fbe8fff477455198185eea0eda2d6"), "lbsame_rev1": (RS / "samepop/lb_same_paired_rev1.py", "9641fa7e6adce94cca2a89ec5b66a8fa2d567007c77412e5e4d17fd4a8929b58"),
       "bbbc_s11": (RS / ".." / "_opt_20261021/s11.py", "c7d3ad6228254855dde16d08087b47d911558cd33869d12e769617daa36fca6f"), "mh_s11": (RS / "mh/mh_s11.py", "8d294b6a09e97eba545df0f6a346bad05022b71e35aa802492c2a4ba98b334a1"),
       "lh_s11": (RS / "paired/lh_paired_s11.py", "669f112b57b10400974ed19a2245ca2796560f6ac90e79b9f6defbaaea1e873c")}
PRELUDE = 'import os as _os; LAB = _os.environ["C5B_ARM"]; LOGDIR = Path(_os.environ["C5B_LOGDIR"]); assert LAB in ("a2", "c5b-f", "c5b-fs"), LAB   # C5b: arm label of the files; ARMS keeps ONE entry so that ai = 0 (the A2 streams)\n'
ARMS_OLD = 'ARMS = ["a2", "a0-ctrl"];'; ARMS_NEW = PRELUDE + 'ARMS = ["a2"];'
REC_OLD = 'e["cert"].append(cert); e["clicks"].append(clicks); e["nlab"].append(nlab)'; REC_NEW = REC_OLD + '; e["kk"].append(k)'
REC_INIT_OLD = '"cert": [], "clicks": [], "nlab": []}'; REC_INIT_NEW = '"cert": [], "clicks": [], "nlab": [], "kk": []}'
BLOCK = '''    # ---- per-draw output (as the PAIRED per-draw block of build_paired.py; the file name carries the arm label): writes the per-draw lists that rec() has already accumulated (cert, fpr, sens, clicks, selected candidate); nothing is computed here
    _stem = Path(__file__).stem
    for _arm in ARMS:
        _o = {}
        for _key, _e in merged.items():
            if _key[0] != _arm or _key[2] != 22: continue
            _tag = "|".join(str(_x) for _x in _key[3:])
            for _f in ("cert", "fpr", "sens", "clicks"): _o[_tag + "|" + _f] = np.array(_e[_f])
            if "kk" in _e and _tag.startswith("S3"): _o[_tag + "|sel"] = np.array([_k if _c else "none" for _k, _c in zip(_e["kk"], _e["cert"])])
            if "sel" in _e: _o[_tag + "|sel"] = np.array([CEN[_i] if _i >= 0 else "none" for _i in _e["sel"]])
        np.savez(OUT / (_stem + "_perdraw_" + _arm + "_m22.npz"), **_o)
'''
PD_OLD = '"_perdraw_" + _arm + "_m22.npz"'; PD_NEW = '"_perdraw_" + LAB + "_m22.npz"'
SITECHECK_BEGIN = "    # ---- PAIRED site check"; SITECHECK_END = "    tasks = [(ai, sd, si)"


def one(t, old, new, name):
    if t.count(old) != 1: raise SystemExit("edit not unique in %s: %r (%d)" % (name, old[:70], t.count(old)))
    return t.replace(old, new)


def build(key):
    path, sha = SRC[key]; raw = path.read_bytes(); assert hashlib.sha256(raw).hexdigest() == sha, "source changed: %s" % path
    t = raw.decode(); stem = "c5b_" + key; add_perdraw = key in ("bbbc_rev1", "mh_rev1", "bbbc_s11", "mh_s11")
    t = one(t, ARMS_OLD, ARMS_NEW, key)
    t = one(t, "for ai in range(2) for sd in SEEDS", "for ai in range(1) for sd in SEEDS", key)
    if key == "bbbc_rev1":
        t = one(t, 'SITES = ["bbbc041_all", "mpidb_full", "lacuna_B", "plasmo"]', 'SITES = ["bbbc041_all"]', key)
        t = one(t, 'np.load(INF/f"{arm}_s{sd}_scores.npz")[site+"_logits"]', 'np.load(LOGDIR/f"{LAB}_s{sd}_{site}_logits.npz")["logits"]', key)
        t = one(t, 'OUT/"rev1_results.json"', 'OUT/("c5b_bbbc_rev1_results_%s.json" % LAB)', key)
    if key == "mh_rev1":
        t = one(t, 'np.load(E8/f"{arm}_s{sd}_mh_logits.npz")["logits"]', 'np.load(LOGDIR/f"{LAB}_s{sd}_mh_logits.npz")["logits"]', key)
        t = one(t, 'OUT/"mh_rev1_results.json"', 'OUT/("c5b_mh_rev1_results_%s.json" % LAB)', key)
    if key == "lh_rev1":
        t = one(t, 'np.load(E8/f"{arm}_s{sd}_lh_logits.npz")["logits"]', 'np.load(LOGDIR/f"{LAB}_s{sd}_lh_logits.npz")["logits"]', key)
        t = one(t, 'OUT/"lh_rev1_results.json"', 'OUT/("c5b_lh_rev1_results_%s.json" % LAB)', key)
        i = t.index(SITECHECK_BEGIN); j = t.index(SITECHECK_END); assert t.count(SITECHECK_BEGIN) == 1 and t.count(SITECHECK_END) == 1 and i < j
        t = t[:i] + "    # (the OLD/NEW negative-score site check of lh_paired_rev1.py is removed: the C5b models are different models, so their scores differ by design; check_c5b.py compares rows, fields and labels instead)\n" + t[j:]
        t = one(t, PD_OLD, PD_NEW, key)
    if key == "lbsame_rev1":
        t = one(t, 'np.load(E8/f"{arm}_s{sd}_lb_same_logits.npz")["logits"]', 'np.load(LOGDIR/f"{LAB}_s{sd}_lb_same_logits.npz")["logits"]', key)
        t = one(t, 'OUT/"lb_same_rev1_results.json"', 'OUT/("c5b_lbsame_rev1_results_%s.json" % LAB)', key)
        t = one(t, PD_OLD, PD_NEW, key)
    if key == "bbbc_s11":
        t = one(t, 'SITES = ["bbbc041_all", "mpidb_full", "lacuna_B"]', 'SITES = ["bbbc041_all"]', key)
        t = one(t, "for si in range(3)]", "for si in range(len(SITES))]", key)
        t = one(t, '''    ids = np.load(INF/"ids_iter1.npz", allow_pickle=True); lid = np.load(E7/"ids_lacuna_B.npz", allow_pickle=True)
    if site == "lacuna_B":
        z = np.load(E7/f"{arm}_s{sd}_lacuna_B_logits.npz")["logits"]; s = z[:, 1]-z[:, 0]; y = lid["y"]; g = lid["field"]
    else:
        z = np.load(INF/f"{arm}_s{sd}_scores.npz")[site+"_logits"]; s = z[:, 1]-z[:, 0]; y = ids[site+"_y"]; g = ids[site+"_group"]
''', '''    assert site == "bbbc041_all"
    ids = np.load(INF/"ids_iter1.npz", allow_pickle=True)
    z = np.load(LOGDIR/f"{LAB}_s{sd}_{site}_logits.npz")["logits"]; s = z[:, 1]-z[:, 0]; y = ids[site+"_y"]; g = ids[site+"_group"]
''', key)
        t = one(t, "with Pool(8) as p:", "with Pool(4) as p:", key)                                            # worker count only (CPU budget of the C5b evaluation: 4 workers); no effect on any value
        t = one(t, 'OUT/"s11_results.json"', 'OUT/("c5b_bbbc_s11_results_%s.json" % LAB)', key)
    if key == "mh_s11":
        t = one(t, 'E7/f"{arm}_s{sd}_mh_logits.npz"', 'LOGDIR/f"{LAB}_s{sd}_mh_logits.npz"', key)
        t = one(t, 'OUT/"mh_s11_results.json"', 'OUT/("c5b_mh_s11_results_%s.json" % LAB)', key)
    if key == "lh_s11":
        t = one(t, 'E7/f"{arm}_s{sd}_lh_logits.npz"', 'LOGDIR/f"{LAB}_s{sd}_lh_logits.npz"', key)
        t = one(t, 'OUT/"lh_s11_results.json"', 'OUT/("c5b_lh_s11_results_%s.json" % LAB)', key)
        t = one(t, PD_OLD, PD_NEW, key)
    if add_perdraw:
        if "_rev1" in key:
            t = one(t, REC_OLD, REC_NEW, key); t = one(t, REC_INIT_OLD, REC_INIT_NEW, key)
        if not t.endswith("\n"): t += "\n"
        t += BLOCK.replace(PD_OLD, PD_NEW)
    hdr = "# C5b EVAL (PREREG_P5_C5B_FRAMING_AUG_20261026.md + errata 1-3): derived from %s (sha256 %s...) by build_c5b_eval.py. ARMS = ['a2'] keeps ai = 0 (the frozen A2's random streams); the arm label LAB is read from env C5B_ARM (c5b-f / c5b-fs), logits from env C5B_LOGDIR/<LAB>_s<seed>_<site>_logits.npz; output file names carry LAB; see the build script for the full list of edits.\n" % (path.name, sha[:12])
    return hdr + t, stem, raw.decode()


if __name__ == "__main__":
    for key in SRC:
        t, stem, src = build(key); dst = Path(stem + ".py")
        if dst.exists(): raise SystemExit("refusing to overwrite " + str(dst))
        dst.write_text(t); d = "".join(difflib.unified_diff(src.splitlines(True), t.splitlines(True), fromfile=SRC[key][0].name, tofile=dst.name)); Path("diffs").mkdir(exist_ok=True); Path("diffs", dst.name + ".diff").write_text(d)
        print(dst.name, hashlib.sha256(t.encode()).hexdigest()[:16], "changed lines", sum(1 for l in d.split("\n") if l[:1] in "+-" and l[:3] not in ("+++", "---")))

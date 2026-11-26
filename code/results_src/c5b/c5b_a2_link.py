"""Builds the unified-layout directory of the FROZEN A2 logits for the eval_a2 step of run_c5b_eval.sh (main-session instruction of 2026-10-27): <out>/a2_s<seed>_<site>_logits.npz for seeds 41-45 and the sites bbbc041_all, mh, lh, lb_same,
plus ids_mh.npz / ids_lh.npz / ids_lb_same.npz (the A2 ids the earlier runs read). BBBC: the key bbbc041_all_logits of <arch>/a2_s<seed>_scores.npz is written to a2_s<seed>_bbbc041_all_logits.npz (array copied unchanged, key 'logits');
mh / lh / lb_same: symbolic links to the existing files. MANIFEST.tsv records source path, source sha256, linked/written file sha256. Nothing is printed except the number of files; no score or label is read into any output other than the copied BBBC array.
Refuses to write into an existing directory. Usage: python3 c5b_a2_link.py --arch <iter1 infer dir> --mh-dir D --lh-dir D --lbs-dir D --ids-mh F --ids-lh F --ids-lbs F --out <new dir>"""
import argparse, hashlib, os, sys
from pathlib import Path
import numpy as np
ap = argparse.ArgumentParser()
for k in ("arch", "mh-dir", "lh-dir", "lbs-dir", "ids-mh", "ids-lh", "ids-lbs", "out"): ap.add_argument("--" + k, required=True)
a = ap.parse_args(); out = Path(a.out)
if out.exists(): raise SystemExit("output dir exists: %s" % out)
out.mkdir(parents=True)
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
rows = []; n = 0
for sd in (41, 42, 43, 44, 45):
    src = Path(a.arch) / ("a2_s%d_scores.npz" % sd); dst = out / ("a2_s%d_bbbc041_all_logits.npz" % sd)
    z = np.load(src)["bbbc041_all_logits"]; np.savez(dst, logits=z); rows.append((str(src) + "[bbbc041_all_logits]", sha(src), "written", sha(dst))); n += 1
    for site, d in (("mh", a.mh_dir), ("lh", a.lh_dir), ("lb_same", a.lbs_dir)):
        src = Path(d) / ("a2_s%d_%s_logits.npz" % (sd, site)); assert src.exists(), src
        dst = out / ("a2_s%d_%s_logits.npz" % (sd, site)); os.symlink(str(src.resolve()), str(dst)); rows.append((str(src), sha(src), "linked", sha(dst))); n += 1
for site, f in (("mh", a.ids_mh), ("lh", a.ids_lh), ("lb_same", a.ids_lbs)):
    src = Path(f); assert src.exists(), src; dst = out / ("ids_%s.npz" % site); os.symlink(str(src.resolve()), str(dst)); rows.append((str(src), sha(src), "linked", sha(dst))); n += 1
with open(out / "MANIFEST.tsv", "w") as fh:
    fh.write("source\tsource_sha256\tmode\tfile_sha256\n"); fh.write("".join("\t".join(r) + "\n" for r in rows))
print("a2 link dir: %d files" % n)

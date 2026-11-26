"""Synthetic smoke of check_c5b.py and the site / cache logic of c5b_infer.py (no torch, no real development-site data). Usage: python3 c5b_smoke_misc.py --work <new dir> [--nih-root <local NIH cell_images>] [--ids-iter1 <ids_iter1.npz>]
 (a) check_c5b.py: PASS on a consistent synthetic infer directory; each of five mutations makes it print exactly the expected failing key and exit 3;
 (b) c5b_infer.py --manifests-only on synthetic roots (C5B_SMOKE=1; sizes not asserted) for bbbc041_all, mh, lh, lb_same; cache() on the synthetic crops returns the expected shape and the rows in manifest order;
 (c) with --nih-root and --ids-iter1 (directory listing only, no image is opened): the NIH VALIDATION list of c5b_infer.py equals the archived ids_iter1 names and labels (assertion of the real code path, no SMOKE flag)."""
import argparse, csv, json, os, shutil, subprocess, sys
from pathlib import Path
import numpy as np
from PIL import Image
ap = argparse.ArgumentParser(); ap.add_argument("--work", required=True); ap.add_argument("--nih-root"); ap.add_argument("--ids-iter1"); a = ap.parse_args()
W = Path(a.work); assert not W.exists(); W.mkdir(parents=True); HERE = Path(__file__).resolve().parent; STAGE1 = HERE.parents[1] / "stage1"
rng = np.random.default_rng(3); fails = []
SITES = {"bbbc041_all": 40, "mh": 30, "lh": 36, "lb_same": 36, "nih_validation": 20}
ids = {}
for s, n in SITES.items():
    y = (rng.random(n) < 0.2).astype(int); field = np.array(["F%02d" % (i // 6) for i in range(n)]); ids[s] = dict(y=y, field=field, crop_file=np.array(["f%03d.png" % i for i in range(n)]), kind=np.where(y == 1, "pos", "neg"))
# ---- (a) check_c5b
def build(d, mut=None):
    d.mkdir(); iter1 = dict(bbbc041_all_y=ids["bbbc041_all"]["y"], bbbc041_all_group=ids["bbbc041_all"]["field"], bbbc041_all_name=ids["bbbc041_all"]["crop_file"], nih_validation_y=ids["nih_validation"]["y"], nih_validation_name=np.array(["n%02d.png" % i for i in range(20)]))
    np.savez(d / "ref_iter1.npz", **iter1)
    for s in ("mh", "lh", "lb_same"): np.savez(d / ("ref_%s.npz" % s), **ids[s])
    inf = d / "infer"; inf.mkdir()
    np.savez(inf / "ids_bbbc041_all.npz", **{k: ids["bbbc041_all"][k] for k in ("crop_file", "field", "y")}); np.savez(inf / "ids_nih_validation.npz", name=iter1["nih_validation_name"], y=ids["nih_validation"]["y"], field=ids["nih_validation"]["field"])
    for s in ("mh", "lh", "lb_same"): np.savez(inf / ("ids_%s.npz" % s), **ids[s])
    for g in ("c5b-f", "c5b-fs"):
        for sd in (41, 42, 43, 44, 45):
            for s, n in SITES.items(): np.savez(inf / ("%s_s%d_%s_logits.npz" % (g, sd, s)), logits=rng.normal(size=(n, 2)).astype(np.float32))
            json.dump({"x": 1}, open(inf / ("%s_s%d_eval.json" % (g, sd)), "w"))
    if mut == "y": z = dict(np.load(inf / "ids_lh.npz", allow_pickle=True)); z["y"] = 1 - z["y"]; np.savez(inf / "ids_lh.npz", **z)
    if mut == "missing": (inf / "c5b-fs_s43_mh_logits.npz").unlink()
    if mut == "nan": np.savez(inf / "c5b-f_s41_lb_same_logits.npz", logits=np.full((36, 2), np.nan, np.float32))
    if mut == "fieldorder": z = dict(np.load(inf / "ids_mh.npz", allow_pickle=True)); z["field"] = z["field"][::-1].copy(); np.savez(inf / "ids_mh.npz", **z)
    if mut == "shape": np.savez(inf / "c5b-f_s42_lh_logits.npz", logits=np.zeros((35, 2), np.float32))
def check(d):
    r = subprocess.run([sys.executable, str(HERE / "check_c5b.py"), "--new-dir", str(d / "infer"), "--ids-iter1", str(d / "ref_iter1.npz"), "--ids-mh", str(d / "ref_mh.npz"), "--ids-lh", str(d / "ref_lh.npz"), "--ids-lb-same", str(d / "ref_lb_same.npz")], capture_output=True, text=True)
    return r.returncode, r.stdout.strip()
build(W / "ok"); rc, out = check(W / "ok"); print("check_c5b consistent:", out, "rc", rc)
if (rc, out) != (0, "PASS"): fails.append("check_c5b did not PASS on a consistent directory")
for mut, key in (("y", "lh:y"), ("missing", "c5b-fs_s43:mh:missing"), ("nan", "c5b-f_s41:lb_same:nonfinite"), ("fieldorder", "mh:field_list_or_row_fields"), ("shape", "c5b-f_s42:lh:shape")):
    build(W / mut, mut); rc, out = check(W / mut); print("check_c5b mutation %-10s -> rc %d, printed %r" % (mut, rc, out))
    if rc != 3 or out != key: fails.append("mutation %s: expected %r, got rc %d %r" % (mut, key, rc, out))
# ---- (b) c5b_infer site logic on synthetic roots
R = W / "roots"; R.mkdir()
def png(p, i): Image.fromarray(np.random.default_rng(i).integers(0, 255, (40 + i % 7, 44, 3), dtype=np.uint8)).save(p)
b = R / "bbbc"; b.mkdir(); (b / "Uninfected").mkdir(); (b / "Parasitized").mkdir(); rows = []
for i in range(SITES["bbbc041_all"]):
    cls = "Parasitized" if ids["bbbc041_all"]["y"][i] else "Uninfected"; f = "%s/c%03d.png" % (cls, i); png(b / f, i); rows.append((f, ids["bbbc041_all"]["field"][i], int(ids["bbbc041_all"]["y"][i])))
with open(b / "manifest.csv", "w", newline="") as fh: w = csv.writer(fh); w.writerow(["file", "image_uuid", "label"]); w.writerows(rows)
m = R / "mh"; (m / "Parasitized").mkdir(parents=True); (m / "Uninfected").mkdir(); rows = []
for i in range(SITES["mh"]):
    cls = "Parasitized" if ids["mh"]["y"][i] else "Uninfected"; fld = ids["mh"]["field"][i]; png(m / cls / ("%s_c%d.png" % (fld, i)), i); rows.append((fld, i, cls))
with open(m / "manifest.csv", "w", newline="") as fh: w = csv.writer(fh); w.writerow(["field", "instance", "class"]); w.writerows(rows)
h = R / "lh"; h.mkdir(); rows = []
for i in range(SITES["lh"]):
    f = "h%03d.png" % i; png(h / f, i); rows.append((f, "B", ids["lh"]["field"][i], ids["lh"]["kind"][i], i, int(ids["lh"]["y"][i]), "0"))
with open(h / "manifest.csv", "w", newline="") as fh: w = csv.writer(fh); w.writerow(["file", "subsite", "field_id", "kind", "idx", "label", "tuning"]); w.writerows(rows)
o = R / "old"; o.mkdir()
for i in range(SITES["lb_same"]): png(o / ids["lb_same"]["crop_file"][i], i)
np.savez(R / "ids_lb_same.npz", **ids["lb_same"])
env = dict(os.environ, C5B_SMOKE="1", P5_CODE=str(STAGE1), C5B_BBBC_ROOT=str(b), C5B_MH_ROOT=str(m), C5B_LH_ROOT=str(h), C5B_OLD_ROOT=str(o))
r = subprocess.run([sys.executable, str(HERE / "c5b_infer.py"), "--out", str(W / "x"), "--manifests-only", "--sites", "bbbc041_all,mh,lh,lb_same", "--ids-iter1", str(W / "ok" / "ref_iter1.npz"), "--ids-lb-same", str(R / "ids_lb_same.npz")], env=env, capture_output=True, text=True)
print(r.stdout.strip()); 
if r.returncode != 0 or "manifests ok" not in r.stdout: fails.append("c5b_infer --manifests-only failed: " + r.stderr[-300:])
code = ("import os,sys,numpy as np; sys.path.insert(0,%r); os.environ['NPROC']='2'; import c5b_infer as M; M.ROOTS['bbbc']=%r\n"
        "a=type('A',(),{'ids_iter1':None,'ids_lb_same':None})(); d=M.site_rows('bbbc041_all',a); c=M.cache(d['paths'][:12],True); print(c.shape,c.dtype)\n"
        "from PIL import Image; import numpy as np; ref=np.asarray(Image.open(d['paths'][5]).convert('RGB').resize((128,128)),dtype=np.uint8)\n"
        "from p5stage.prod import legacy_float_pr_od; r=np.asarray(Image.fromarray(legacy_float_pr_od(np.asarray(Image.open(d['paths'][5]).convert('RGB'),dtype=np.uint8))).resize((128,128),Image.BILINEAR),dtype=np.uint8).transpose(2,0,1)\n"
        "print('row 5 equals the single-image pipeline:', bool((c[5]==r).all()))") % (str(HERE), str(b))
r = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True); print(r.stdout.strip(), r.stderr[-300:] if r.returncode else "")
if "(12, 3, 128, 128) uint8" not in r.stdout or "equals the single-image pipeline: True" not in r.stdout: fails.append("cache() shape/order check failed")
# ---- (c) NIH validation list against the archived ids (directory listing only)
if a.nih_root and a.ids_iter1:
    env2 = dict(os.environ, P5_CODE=str(STAGE1), C5B_NIH_ROOT=a.nih_root); env2.pop("C5B_SMOKE", None)
    r = subprocess.run([sys.executable, str(HERE / "c5b_infer.py"), "--out", str(W / "y"), "--manifests-only", "--sites", "nih_validation", "--ids-iter1", a.ids_iter1], env=env2, capture_output=True, text=True)
    print("NIH validation list vs ids_iter1 (names, labels, n=4081):", r.stdout.strip().replace("\n", " | "), "" if r.returncode == 0 else r.stderr[-300:])
    if r.returncode != 0 or "manifests ok" not in r.stdout: fails.append("NIH validation identity failed")
print("RESULT:", "ALL PASSED" if not fails else "FAILED"); [print("  ", x) for x in fails]; sys.exit(1 if fails else 0)

"""R8 step 0 (read-only): overlap of IML Malaria images with other datasets on the VM, by file name and md5 (md5 only for files whose size equals an IML file size).
Usage: overlap_vm.py IML_IMAGE_DIR OUT_JSON. Exit 3 if any name or md5 hit. Detects byte-identical copies only (re-encoded copies are not detected)."""
import hashlib, io, json, os, sys, zipfile, collections
IML = sys.argv[1]; OUT = sys.argv[2]; EXT = (".jpg", ".jpeg", ".png", ".tif", ".tiff")
B = "<DATA>/medprune/MedPrune/datasets"
DIRS = {"bbbc041": B+"/bbbc041", "mpidb_cells": B+"/mpidb_cells", "plasmo": "<DATA>/p5_e8_candidate/images", "lacuna": "<DATA>/p5_e7_lacuna/extracted/images", "nih": B+"/malaria/cell_images"}
ZIPS = {"mpidb_zip": B+"/mpidb.zip"}
md5 = lambda b: hashlib.md5(b).hexdigest()
iml = {}; sizes = collections.defaultdict(list)
for f in sorted(os.listdir(IML)):
    if not f.lower().endswith(".jpg"): continue
    d = open(os.path.join(IML, f), "rb").read(); iml[f] = md5(d); sizes[len(d)].append(f)
hmd5 = {v: k for k, v in iml.items()}; names = {k.lower(): k for k in iml}
res = {"iml_images": len(iml), "datasets": {}, "name_hits": [], "md5_hits": []}
def check(ds, rel, size, reader, st):
    st["scanned"] += 1; base = os.path.basename(rel).lower()
    if base in names: res["name_hits"].append({"dataset": ds, "file": rel, "iml": names[base]})
    if size in sizes:
        st["same_size"] += 1; h = md5(reader())
        if h in hmd5: res["md5_hits"].append({"dataset": ds, "file": rel, "iml": hmd5[h]})
for ds, root in DIRS.items():
    st = {"root": root, "scanned": 0, "same_size": 0}
    for r, _, fs in os.walk(root):
        for f in fs:
            if f.lower().endswith(EXT):
                p = os.path.join(r, f); check(ds, p, os.path.getsize(p), lambda p=p: open(p, "rb").read(), st)
    res["datasets"][ds] = st
for ds, zp in ZIPS.items():
    st = {"root": zp, "scanned": 0, "same_size": 0}
    with zipfile.ZipFile(zp) as z:
        for zi in z.infolist():
            if zi.filename.lower().endswith(EXT): check(ds, zi.filename, zi.file_size, lambda zi=zi, z=z: z.read(zi), st)
    res["datasets"][ds] = st
res["any_hit"] = bool(res["name_hits"] or res["md5_hits"])
res["limitation"] = "byte-identical copies and identical file names only; NIH and mpidb_cells hold cell crops, not fields"
json.dump(res, open(OUT, "w"), indent=1); print(json.dumps({k: v for k, v in res.items() if k != "limitation"}, indent=1))
sys.exit(3 if res["any_hit"] else 0)

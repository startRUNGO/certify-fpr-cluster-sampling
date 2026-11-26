"""R8 IML Malaria crop per PREREG_R8_IML_20261023.md (sha256 in .sha256 file). Adapted from e8_crop.py: only JSON reading, paths and class mapping changed.
positives = gametocyte/ring/trophozoite/schizont; negatives = red blood cell; difficult excluded. PAD=5, clipped to image. Stops on non-RGB/non-8-bit images."""
import csv, json, os, sys, collections
import numpy as np
from PIL import Image
D = "<DATA>/p5_e8_candidate/r8_iml"; OUT = sys.argv[1]; PAD = 5
POS = ("gametocyte", "ring", "trophozoite", "schizont"); NEG = ("red blood cell",)
ann = json.load(open(D+"/annotations.json")); have = set(os.listdir(D+"/images")); items = [a for a in ann if a["image_name"] in have]
os.makedirs(OUT+"/pos", exist_ok=True); os.makedirs(OUT+"/neg", exist_ok=True)
man = []; cnt = collections.Counter(); modes = collections.Counter(); sizes = collections.Counter(); per_unin = []; per_inf = []; oob = 0
for it in items:
    f = it["image_name"]; im = Image.open(D+"/images/"+f)
    if im.mode != "RGB": print("STOP non-RGB mode", im.mode, f); sys.exit(2)
    modes[im.mode] += 1; sizes[im.size] += 1; W, H = im.size; arr = np.asarray(im)
    nu = ni = 0
    for j, o in enumerate(it["objects"]):
        cat = o["type"]; cnt[cat] += 1
        if cat not in POS and cat not in NEG: continue
        b = {k: float(v) for k, v in o["bbox"].items()}
        c0 = int(max(b["x"]-PAD, 0)); c1 = int(min(b["x"]+b["w"]+PAD, W)); r0 = int(max(b["y"]-PAD, 0)); r1 = int(min(b["y"]+b["h"]+PAD, H))
        if r1-r0 < 2 or c1-c0 < 2: oob += 1; continue
        kind = "pos" if cat in POS else "neg"; name = f"{kind}/{os.path.splitext(f)[0]}_{j:04d}.png"
        Image.fromarray(arr[r0:r1, c0:c1]).save(f"{OUT}/{name}", compress_level=1); man.append((name, os.path.splitext(f)[0], kind, int(kind == "pos"), c0, r0, c1, r1, cat))
        if kind == "pos": ni += 1
        else: nu += 1
    per_unin.append(nu); per_inf.append(ni)
with open(OUT+"/manifest.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["file", "field_id", "kind", "label", "x0", "y0", "x1", "y1", "orig_class"]); w.writerows(man)
q = lambda v: [float(x) for x in np.quantile(v, [0, .25, .5, .75, 1])]
rep = {"images": len(items), "images_expected_in_json": len(ann), "modes": dict(modes), "sizes": {str(k): v for k, v in sizes.items()}, "json_box_counts": dict(cnt), "crops_pos": sum(1 for m in man if m[2] == "pos"), "crops_neg": sum(1 for m in man if m[2] == "neg"), "degenerate_boxes_skipped": oob,
       "uninfected_per_image_min_q1_med_q3_max": q(per_unin), "infected_per_image_min_q1_med_q3_max": q(per_inf), "images_without_infected": sum(1 for x in per_inf if x == 0)}
json.dump(rep, open(OUT+"/counts.json", "w"), indent=1); print(json.dumps(rep, indent=1))

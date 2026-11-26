"""E8 PlasmoCount crop per PREREG_E8_PLASMOCOUNT_20261023.md (sha256 2ce88dfd...). infected=pos, uninfected=neg (annotated), not_sure excluded. PAD=5 as in E7. Stops on non-8-bit images."""
import csv, json, os, sys, collections, statistics as st
import numpy as np
from PIL import Image
D = "<DATA>/p5_e8_candidate"; OUT = sys.argv[1]; PAD = 5
test = json.load(open(D+"/test.json")); have = set(os.listdir(D+"/images")); items = [it for it in test if it["image"]["pathname"] in have]
os.makedirs(OUT+"/pos", exist_ok=True); os.makedirs(OUT+"/neg", exist_ok=True)
man = []; cnt = collections.Counter(); modes = collections.Counter(); sizes = collections.Counter(); per_unin = []; per_inf = []; oob = 0
for it in items:
    f = it["image"]["pathname"]; im = Image.open(D+"/images/"+f)
    if im.mode not in ("RGB", "RGBA"): print("STOP non-8-bit/unsupported mode", im.mode, f); sys.exit(2)
    modes[im.mode] += 1; sizes[im.size] += 1; W, H = im.size; arr = np.asarray(im.convert("RGB"))
    nu = ni = 0
    for j, o in enumerate(it["objects"]):
        cat = o["category"]; cnt[cat] += 1
        if cat not in ("infected", "uninfected"): continue
        b = o["bounding_box"]; r0 = max(b["minimum"]["r"]-PAD, 0); r1 = min(b["maximum"]["r"]+PAD, H); c0 = max(b["minimum"]["c"]-PAD, 0); c1 = min(b["maximum"]["c"]+PAD, W)
        if r1-r0 < 2 or c1-c0 < 2: oob += 1; continue
        kind = "pos" if cat == "infected" else "neg"; name = f"{kind}/{os.path.splitext(f)[0]}_{j:04d}.png"
        Image.fromarray(arr[r0:r1, c0:c1]).save(f"{OUT}/{name}", compress_level=1); man.append((name, os.path.splitext(f)[0], kind, int(kind == "pos"), c0, r0, c1, r1))
        if kind == "pos": ni += 1
        else: nu += 1
    per_unin.append(nu); per_inf.append(ni)
with open(OUT+"/manifest.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["file", "field_id", "kind", "label", "x0", "y0", "x1", "y1"]); w.writerows(man)
q = lambda v: [float(x) for x in np.quantile(v, [0, .25, .5, .75, 1])]
rep = {"images": len(items), "images_expected_in_json": len(test), "modes": dict(modes), "sizes": {str(k): v for k, v in sizes.items()}, "json_box_counts": dict(cnt), "crops_pos": sum(1 for m in man if m[2] == "pos"), "crops_neg": sum(1 for m in man if m[2] == "neg"), "degenerate_boxes_skipped": oob,
       "uninfected_per_image_min_q1_med_q3_max": q(per_unin), "infected_per_image_min_q1_med_q3_max": q(per_inf), "images_without_infected": sum(1 for x in per_inf if x == 0)}
json.dump(rep, open(OUT+"/counts.json", "w"), indent=1); print(json.dumps(rep, indent=1))

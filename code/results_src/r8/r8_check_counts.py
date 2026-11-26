"""R8 step 2 (pre-scoring report check): compare crops/counts.json with the counts read from annotations.json before the preregistration (PREREG_R8_IML_20261023.md sec. 1). Exit 4 on any mismatch."""
import json, sys
c = json.load(open(sys.argv[1])); bad = []
exp = {"images": 345, "images_expected_in_json": 345, "modes": {"RGB": 345}, "sizes": {"(1280, 960)": 345}, "crops_pos": 529, "crops_neg": 37899, "degenerate_boxes_skipped": 0,
       "json_box_counts": {"red blood cell": 37899, "gametocyte": 261, "ring": 164, "trophozoite": 77, "schizont": 27, "difficult": 21},
       "uninfected_per_image_min_q1_med_q3_max_minmedmax": [14.0, 110.0, 164.0], "images_without_infected": 32}
for k, v in exp.items():
    if k.endswith("_minmedmax"): u = c["uninfected_per_image_min_q1_med_q3_max"]; got = [u[0], u[2], u[4]]
    else: got = c.get(k)
    if got != v: bad.append((k, got, v))
print("PRE-SCORING REPORT", json.dumps(c, indent=1)); print("MISMATCHES", bad)
sys.exit(4 if bad else 0)

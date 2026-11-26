"""C5b check (replaces the 'negatives unchanged' assertion of paired_check*.py, which does not apply: the C5b models are different models). Compares, for every site, the ids written by c5b_infer.py with the ids the frozen A2 run used:
row count, number of fields, labels (y), the field of every row (hence the field list and its order), and where the frozen file has them the crop file names and kinds; and checks that all 10 x 5 logits files exist with shape (n, 2) and only finite values.
Prints PASS or, on failure, only the failing keys (one per line) and exits with code 3. Prints no value of any score or label.
Usage: python3 check_c5b.py --new-dir <c5b_infer output> --ids-iter1 <ids_iter1.npz> --ids-mh <A2 ids_mh.npz> --ids-lh <A2 ids_lh.npz> --ids-lb-same <ids_lb_same.npz>"""
import argparse, sys
from pathlib import Path
import numpy as np
ap = argparse.ArgumentParser()
for k in ("new-dir", "ids-iter1", "ids-mh", "ids-lh", "ids-lb-same"): ap.add_argument("--" + k, required=True)
a = ap.parse_args(); new = Path(a.new_dir); bad = []
STEMS = ["%s_s%d" % (g, s) for g in ("c5b-f", "c5b-fs") for s in (41, 42, 43, 44, 45)]
def S(x): return [str(v) for v in x]
def cmp(site, got, ref_y, ref_field, ref_extra=None):
    ok = len(got["y"]) == len(ref_y)
    if not ok: bad.append(site + ":rows"); return
    if len(set(S(got["field"]))) != len(set(S(ref_field))): bad.append(site + ":n_fields")
    if not np.array_equal(np.asarray(got["y"]), np.asarray(ref_y)): bad.append(site + ":y")
    if S(got["field"]) != S(ref_field): bad.append(site + ":field_list_or_row_fields")
    for k, v in (ref_extra or {}).items():
        if k not in got.files or S(got[k]) != S(v): bad.append(site + ":" + k)
iter1 = np.load(a.ids_iter1, allow_pickle=True)
g = np.load(new / "ids_bbbc041_all.npz", allow_pickle=True)
cmp("bbbc041_all", g, iter1["bbbc041_all_y"], iter1["bbbc041_all_group"])
if S(Path(str(n)).name for n in g["crop_file"]) != S(iter1["bbbc041_all_name"]): bad.append("bbbc041_all:crop_file_names")
g = np.load(new / "ids_nih_validation.npz", allow_pickle=True)
if S(g["name"]) != S(iter1["nih_validation_name"]): bad.append("nih_validation:names")
if not np.array_equal(np.asarray(g["y"]), iter1["nih_validation_y"]): bad.append("nih_validation:y")
for site, ref in (("mh", a.ids_mh), ("lh", a.ids_lh), ("lb_same", a.ids_lb_same)):
    r = np.load(ref, allow_pickle=True); g = np.load(new / ("ids_%s.npz" % site), allow_pickle=True)
    cmp(site, g, r["y"], r["field"], {"crop_file": r["crop_file"], "kind": r["kind"]})
NS = {"nih_validation": len(iter1["nih_validation_y"]), "bbbc041_all": len(iter1["bbbc041_all_y"]), "mh": len(np.load(a.ids_mh, allow_pickle=True)["y"]), "lh": len(np.load(a.ids_lh, allow_pickle=True)["y"]), "lb_same": len(np.load(a.ids_lb_same, allow_pickle=True)["y"])}
for st in STEMS:
    for site, n in NS.items():
        f = new / ("%s_%s_logits.npz" % (st, site))
        if not f.exists(): bad.append("%s:%s:missing" % (st, site)); continue
        z = np.load(f)["logits"]
        if z.shape != (n, 2): bad.append("%s:%s:shape" % (st, site))
        elif not np.isfinite(z).all(): bad.append("%s:%s:nonfinite" % (st, site))
    if not (new / ("%s_eval.json" % st)).exists(): bad.append("%s:eval_json_missing" % st)
if bad: print("\n".join(sorted(set(bad)))); sys.exit(3)
print("PASS")

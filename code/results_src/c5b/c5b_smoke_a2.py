"""Smoke of c5b_a2_link.py and c5b_a2_compare.py on the synthetic world written by c5b_smoke_eval.py (--world = its --work directory; nothing real is read).
 link: the BBBC array written equals the A2 scores array, the mh / lh / lb_same files and the ids are symbolic links that resolve to the sources, the MANIFEST is complete, an existing output directory is refused;
 compare: PASS on the genuine pairs (c5b scripts on A2 logits vs the reviewed originals, as the wrapper builds them); a perturbation of 5e-10 still passes, one of 1e-6 prints exactly the perturbed key; a changed per-draw array and a dropped key are caught.
Usage: python3 c5b_smoke_a2.py --world <dir> --work <new dir>"""
import argparse, json, subprocess, sys, shutil
from pathlib import Path
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--world", required=True); ap.add_argument("--work", required=True); a = ap.parse_args()
Wd = Path(a.world); W = Path(a.work); W.mkdir(parents=True); HERE = Path(__file__).resolve().parent; fails = []
# ---- link
out = W / "a2_link"
r = subprocess.run([sys.executable, str(HERE / "c5b_a2_link.py"), "--arch", str(Wd / "orig_inf"), "--mh-dir", str(Wd / "orig_e"), "--lh-dir", str(Wd / "orig_e"), "--lbs-dir", str(Wd / "orig_e"), "--ids-mh", str(Wd / "orig_e/ids_mh.npz"),
                    "--ids-lh", str(Wd / "orig_e/ids_lh.npz"), "--ids-lbs", str(Wd / "orig_e/ids_lb_same.npz"), "--out", str(out)], capture_output=True, text=True)
print("link:", r.stdout.strip(), r.stderr[-200:])
if r.returncode != 0: fails.append("link failed")
for sd in (41, 42, 43, 44, 45):
    for site in ("bbbc041_all", "mh", "lh", "lb_same"):
        f = out / ("a2_s%d_%s_logits.npz" % (sd, site)); ref = np.load(Wd / "uni_a2" / ("a2_s%d_%s_logits.npz" % (sd, site)))["logits"]
        if not f.exists() or not np.array_equal(np.load(f)["logits"], ref): fails.append("link array differs: %s" % f.name)
        if site != "bbbc041_all" and not (f.is_symlink() and f.resolve().exists()): fails.append("not a resolving symlink: %s" % f.name)
for s in ("mh", "lh", "lb_same"):
    if not ((out / ("ids_%s.npz" % s)).is_symlink() and (out / ("ids_%s.npz" % s)).resolve().exists()): fails.append("ids link " + s)
if len((out / "MANIFEST.tsv").read_text().strip().split("\n")) != 1 + 20 + 3: fails.append("manifest rows")
r2 = subprocess.run([sys.executable, str(HERE / "c5b_a2_link.py"), "--arch", str(Wd / "orig_inf"), "--mh-dir", str(Wd / "orig_e"), "--lh-dir", str(Wd / "orig_e"), "--lbs-dir", str(Wd / "orig_e"), "--ids-mh", str(Wd / "orig_e/ids_mh.npz"),
                     "--ids-lh", str(Wd / "orig_e/ids_lh.npz"), "--ids-lbs", str(Wd / "orig_e/ids_lb_same.npz"), "--out", str(out)], capture_output=True, text=True)
if r2.returncode == 0: fails.append("existing output dir was not refused")
# ---- compare
NEWJ = {"bbbc_rev1": ("c5b_bbbc_rev1_results_a2.json", "rev1_results.json", "a2|bbbc041_all|", "a2|bbbc041_all"), "bbbc_s11": ("c5b_bbbc_s11_results_a2.json", "s11_results.json", "a2|bbbc041_all|", None),
        "mh_rev1": ("c5b_mh_rev1_results_a2.json", "mh_rev1_results.json", "a2|mh|", "a2|mh"), "mh_s11": ("c5b_mh_s11_results_a2.json", "mh_s11_results.json", "a2|mh|", None),
        "lh_rev1": ("c5b_lh_rev1_results_a2.json", "lh_rev1_results.json", "a2|lh|", "a2|lh"), "lh_s11": ("c5b_lh_s11_results_a2.json", "lh_s11_results.json", "a2|lh|", None),
        "lbsame_rev1": ("c5b_lbsame_rev1_results_a2.json", "lb_same_rev1_results.json", "a2|lb_same|", "a2|lb_same")}
def pairs(newroot):
    P = []
    for k, (nf, rf, pre, p1) in NEWJ.items(): P.append(dict(label=k, new=str(newroot / k / nf), ref=str(Wd / "out_ref" / k / rf), prefix=pre, prop1_prefix=p1))
    for k in ("lh_rev1", "lh_s11", "lbsame_rev1"): P.append(dict(label=k + "_perdraw", new=str(newroot / k / ("c5b_%s_perdraw_a2_m22.npz" % k)), ref=str(Wd / "out_ref" / k / ("ref_%s_perdraw_a2_m22.npz" % k)), perdraw=True))
    return P
def cmp(P, name):
    f = W / (name + ".json"); json.dump(P, open(f, "w")); r = subprocess.run([sys.executable, str(HERE / "c5b_a2_compare.py"), "--pairs", str(f)], capture_output=True, text=True); return r.returncode, r.stdout.strip()
base = Wd / "out_new_a2"; rc, out_ = cmp(pairs(base), "genuine"); print("compare genuine pairs:", out_, "rc", rc)
if (rc, out_) != (0, "PASS"): fails.append("genuine pairs did not PASS")
# perturbed copy of the new results
mut = W / "mut"; shutil.copytree(base, mut)
f = mut / "lh_rev1" / "c5b_lh_rev1_results_a2.json"; d = json.load(open(f)); k0 = sorted(k for k in d["agg"] if k.startswith("a2|lh|"))[0]; f0 = sorted(d["agg"][k0])[0]
d["agg"][k0][f0] = d["agg"][k0][f0] + 5e-10; json.dump(d, open(f, "w")); rc, o = cmp(pairs(mut), "m5e10"); print("compare +5e-10:", o, "rc", rc)
if (rc, o) != (0, "PASS"): fails.append("5e-10 perturbation should pass")
d["agg"][k0][f0] = d["agg"][k0][f0] + 1e-6; json.dump(d, open(f, "w")); rc, o = cmp(pairs(mut), "m1e6"); print("compare +1e-6:", o, "rc", rc)
if (rc, o) != (3, "lh_rev1:" + k0): fails.append("1e-6 perturbation: expected key %s, got %r" % (k0, o))
d["agg"][k0][f0] -= 1e-6; json.dump(d, open(f, "w")); del d["agg"][k0]; json.dump(d, open(f, "w")); rc, o = cmp(pairs(mut), "dropped"); print("compare dropped key:", o, "rc", rc)
if rc != 3 or "lh_rev1:key_set" not in o: fails.append("dropped key not caught")
z = dict(np.load(mut / "lh_s11" / "c5b_lh_s11_perdraw_a2_m22.npz", allow_pickle=True)); kk = sorted(z)[0]; z[kk] = z[kk].copy(); z[kk].flat[0] = (not z[kk].flat[0]) if z[kk].dtype == bool else z[kk].flat[0] + 1
np.savez(mut / "lh_s11" / "c5b_lh_s11_perdraw_a2_m22.npz", **z); shutil.copy(base / "lh_rev1" / "c5b_lh_rev1_results_a2.json", f); rc, o = cmp(pairs(mut), "perdraw"); print("compare changed per-draw array:", o, "rc", rc)
if (rc, o) != (3, "lh_s11_perdraw:" + kk): fails.append("per-draw change: got %r" % o)
print("RESULT:", "ALL PASSED" if not fails else "FAILED"); [print("  ", x) for x in fails]; sys.exit(1 if fails else 0)

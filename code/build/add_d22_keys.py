# P5a ADD-ONLY results_master builder for draft 22 (实验管理 2026-11-26): (H1) M5 positives, negatives and prevalence for Table 3, and (H7) the Pos-Census simulation grid values for S4.
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  Counts only: no score is read.  Same mechanics as add_sim50_derived3_keys.py (reviewed).
# Input: ids_m5hcm.npz pulled with the M5 inference outputs (array y, 1 = infected, one row per scored crop, manifest order), sha pinned below.
#   This is the label file the M5 confirmation run evaluated (m5_confirm.json n_crops 183105, n_fields 1217); the other sites' Table 3 counts come from their ids / counts files the same way.
# H7: the grid values of the sealed F2 prediction table results_src/p5a2/zpred/predictions.json (already in files_sha256, sha checked against that entry): zplus_pred/pi|0,1,
#   zplus_pred/delta_pos|0..4 (sorted distinct values over the cells), zplus_pred/xi_mid, zplus_pred/n_cells (number of cells); gates: 20 cells = 2 pi x 2 eps x 5 delta_pos, list sizes 12/7/5/4.
# Writes (all 'x', nothing overwritten): <tree>/results_src/m5/m5_counts.json (counts + the ids sha), registered in files_sha256; the master gains
#   dataset/m5|positives, dataset/m5|negatives (source: that file, path [key]) and p5aprev/m5 = positives / cells (as p5aprev/<site> for the other sites).
# Gates (exit 3, nothing written): master sha 168e49de; ids sha 32abb862; y is int64 with values in {0, 1}; len(y) == m5/n_crops; distinct fields == m5/n_fields;
#   new keys absent; write targets absent; after writing, every old key, source and files_sha256 entry unchanged, key count = old + 12, new values re-read equal.
# Backup: results_master_before_m5_counts_168e49de.json ('x').
# Usage: python3 -E -B add_m5_count_keys.py --tree <tree_bj_v1_...> --ids <.../m5_infer_20261031T180153Z/ids_m5hcm.npz>
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
MASTER_SHA = "168e49de2a90082b389f92b6dbb2287e9a1030e38a268bb5b3ea59f8aa4b9019"
IDS_SHA = "32abb8620e8d4ea42c2730bd94b1a1037b04e73019e0659587298f388e445723"
REL = "results_src/m5/m5_counts.json"
NEW = ("dataset/m5|positives", "dataset/m5|negatives", "p5aprev/m5", "zplus_pred/pi|0", "zplus_pred/pi|1", "zplus_pred/xi_mid", "zplus_pred/n_cells") + tuple(f"zplus_pred/delta_pos|{i}" for i in range(5))
ZPRED = "results_src/p5a2/zpred/predictions.json"


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); ap.add_argument("--ids", required=True); a = ap.parse_args()
    tree = Path(a.tree); rm = tree / "results_master.json"; out = tree / REL
    try:
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not 168e49de")
        if sha(a.ids) != IDS_SHA: raise ValueError("ids_m5hcm.npz is not 32abb862")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        clash = [k for k in NEW if k in vals] + ([REL] if REL in fsha else [])
        if clash: raise ValueError(f"keys or file entry exist already: {clash}")
        z = np.load(a.ids, allow_pickle=False); y = z["y"]; fld = z["field"]
        if y.dtype != np.int64 or not set(np.unique(y).tolist()) <= {0, 1}: raise ValueError("y is not int64 in {0, 1}")
        if len(y) != vals["m5/n_crops"] or len(fld) != len(y): raise ValueError(f"rows {len(y)} differ from m5/n_crops {vals['m5/n_crops']}")
        if len(set(fld.tolist())) != vals["m5/n_fields"]: raise ValueError("distinct fields differ from m5/n_fields")
        pos = int(y.sum()); neg = int(len(y) - pos)
        if sha(tree / ZPRED) != fsha.get(ZPRED): raise ValueError("predictions.json differs from its files_sha256 entry")
        zp = json.load(open(tree / ZPRED)); cells = zp["cells"]
        pis = sorted({c_["pi"] for c_ in cells}); dps = sorted({c_["delta_pos"] for c_ in cells}); eps_ = sorted({c_["eps"] for c_ in cells})
        if len(cells) != 20 or len(pis) != 2 or len(dps) != 5 or len(eps_) != 2 or len(cells) != len(pis) * len(dps) * len(eps_): raise ValueError("F2 grid is not 2 x 2 x 5 = 20 cells")
        if [len(zp["lists"][k_]) for k_ in ("target", "clear", "middle", "control")] != [12, 7, 5, 4]: raise ValueError("F2 lists are not 12 / 7 / 5 / 4")
        zsrc = lambda w_: dict(file=ZPRED, derived=w_)
        znew = {"zplus_pred/pi|0": (pis[0], zsrc("smallest distinct pi over the F2 cells")), "zplus_pred/pi|1": (pis[1], zsrc("largest distinct pi over the F2 cells")),
                "zplus_pred/xi_mid": (zp["xi_mid"], dict(file=ZPRED, path=["xi_mid"])), "zplus_pred/n_cells": (len(cells), zsrc("number of F2 cells"))}
        znew.update({f"zplus_pred/delta_pos|{i}": (v_, zsrc(f"distinct delta_pos values over the F2 cells, sorted, index {i}")) for i, v_ in enumerate(dps)})
        counts = {"ids_file": "ids_m5hcm.npz (M5 inference outputs, m5_infer_20261031T180153Z)", "ids_sha256": IDS_SHA, "cells": int(len(y)), "fields": int(len(set(fld.tolist()))),
                  "positives": pos, "negatives": neg}
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    bak = tree / "results_master_before_m5_counts_168e49de.json"; tmp = tree / "results_master.json.m5_counts_tmp"
    existing = [str(x) for x in (bak, tmp, out) if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}; old_fsha = dict(fsha)
    with open(rm, "rb") as s_, open(bak, "xb") as d_: d_.write(s_.read())
    with open(out, "x") as fh: json.dump(counts, fh, indent=1)
    new = {"dataset/m5|positives": (pos, dict(file=REL, path=["positives"])), "dataset/m5|negatives": (neg, dict(file=REL, path=["negatives"])),
           "p5aprev/m5": (pos / len(y), dict(file=REL, derived="positives / cells from m5_counts.json, as p5aprev/<site> for the other sites"))}
    new.update(znew)
    if sorted(new) != sorted(NEW): print("GATE (exit 3): new key set differs from NEW", file=sys.stderr); return 3
    for k, (v, s) in new.items(): vals[k] = v; srcs[k] = s
    fsha[REL] = sha(out); d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp))
    for k, (v, _) in new.items():
        if chk["values"].get(k) != v: print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    for k, v in old_fsha.items():
        if chk["files_sha256"].get(k) != v: print("GATE (exit 3): an old files_sha256 entry changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new) or chk["files_sha256"].get(REL) != sha(out): print("GATE (exit 3): key count or file entry mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    print(json.dumps(dict(added=len(new), n_values=d["n_values"], positives=pos, negatives=neg, prevalence=pos / len(y), new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

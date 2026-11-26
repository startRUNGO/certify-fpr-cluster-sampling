# P5 E10 derived summary keys (pooled fixed-threshold rates, pooled SE and limit, counts).  LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 reviews it.  Reads the E10 output json; the json is not changed.
# Input: --e10 <e10.json> (sha checked against the SHA256SUMS next to it).  Output: a NEW directory with e10_derived.json (+ SHA256SUMS), to be registered in the master as the file stem ebatch_e10_derived (FILES entry in
#   build_results_master.py, added by 实验管理 when the file exists).  The flattened master keys are:
#   ebatch_e10_derived/cells|<cell>|fixed_threshold_pooled|0.1|<kappa>|certify_rate|<arm>     arms B1, <EMAIL>, <EMAIL>/3, <EMAIL>; kappa 1.0 (boundary, type-I) and 0.5 (power); eps = 0.10
#   ebatch_e10_derived/pooled_n|0.1, pooled_se|0.1 = sqrt(0.1 * 0.9 / n), pooled_limit|0.1 = 0.1 + 3 * pooled_se
#   ebatch_e10_derived/count_B1_above_limit|0.1|1.0, count_B1u_above_alpha|0.1|1.0     (strict '>' on the unrounded pooled rates)
# Pooling: every realisation has the same number of draws (asserted), so the equal-weight mean of the realisations' certify rates equals the pooled rate over all draws; n = realisations x draws (asserted = 10000, the number
#   the review text uses).  The counts are over the cells of the design (10 cells at eps = 0.10, kappa = 1.0): B1 above the limit alpha + 3 SE; <EMAIL> above alpha = 0.10.  A rate exactly at the limit is not counted.
# Usage: python e10_derive.py --e10 runs/ebatch_e10_.../e10.json --out <new dir>/e10_derived.json
import argparse, hashlib, json, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "e_batch")); import ebatch_lib as EB
ARMS = ("B1", "<EMAIL>", "<EMAIL>/3", "<EMAIL>"); KAPPAS = ("1.0", "0.5"); EPS = "0.1"; ALPHA = 0.10; N_EXPECT = 10000


def derive(e10):
    cells = {}; ns = set()
    for ci, c in e10["cells"].items():
        reps = c["fixed_threshold"][EPS]; assert len(reps) == e10["reps"], ("realisations per cell", ci, len(reps)); out = {}
        for kappa in KAPPAS:
            per = {}
            for arm in ARMS:
                vals = []
                for sd, r in reps.items():
                    assert set(r[kappa]["certify_rate"]) == set(ARMS), ("arms differ", ci, sd, sorted(r[kappa]["certify_rate"])); vals.append(float(r[kappa]["certify_rate"][arm])); ns.add(int(r[kappa]["n_draws"]))
                per[arm] = float(np.mean(vals))
            out[kappa] = {"certify_rate": per}
        cells[str(ci)] = {"fixed_threshold_pooled": {EPS: out}}
    assert len(ns) == 1, ("draw counts differ between realisations: equal-weight mean is not the pooled rate", sorted(ns)); n = e10["reps"] * ns.pop(); assert n == N_EXPECT, ("pooled n", n)
    se = math.sqrt(ALPHA * (1 - ALPHA) / n); limit = ALPHA + 3 * se
    b1 = [cells[ci]["fixed_threshold_pooled"][EPS]["1.0"]["certify_rate"]["B1"] for ci in cells]; b1u = [cells[ci]["fixed_threshold_pooled"][EPS]["1.0"]["certify_rate"]["<EMAIL>"] for ci in cells]
    return dict(source_sha256=None, cells=cells, pooled_n={EPS: n}, pooled_se={EPS: se}, pooled_limit={EPS: limit}, count_B1_above_limit={EPS: {"1.0": int(sum(v > limit for v in b1))}},
                count_B1u_above_alpha={EPS: {"1.0": int(sum(v > ALPHA for v in b1u))}}, n_cells=len(cells))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--e10", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(); p = Path(a.e10); out = Path(a.out); assert not out.parent.exists(), "output directory exists: use a new directory"
    sums = {l.split()[1].lstrip("*"): l.split()[0] for l in (p.parent / "SHA256SUMS").read_text().splitlines() if l.strip()}; h = EB.sha(p); assert sums.get(p.name) == h, "e10.json differs from its SHA256SUMS"
    e10 = json.load(open(p)); assert e10["version"] == "e10_r1" and e10["spec"].startswith("P5_规格_E10"), "not an E10 output"; res = derive(e10); res["source_sha256"] = h; EB.write_out(out, res)


if __name__ == "__main__":
    main()

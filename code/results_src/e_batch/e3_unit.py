# P5 E3 certification unit: one random negative cell per field (R1-1).  Spec E-batch draft3 section E3.  LOCAL, read-only.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# NATURE: descriptive, development sites only (A2, seeds 41..45).  The other two units (field cluster = S3/SPLIT6, iid cells = R1-all) are reused from the master, not recomputed.
#   M5, IML, PlasmoCount, NIH test, any image, any classifier: never read.
# R1-1: m fields (m in {22, 40}; a site with fewer fields with negatives gives "infeasible", not an extrapolation) are drawn as perm[:m] of a permutation of the fields with >= 1 negative, using the E1 stream
#   default_rng([20261103, 81, site index, seed, eps index, draw]) (so m = 22 uses the same fields as E1 and the nested m = 40 set contains them); from each selected field ONE negative cell is drawn with
#   default_rng([20261103, 81, site index, seed, eps index, draw, 3]); the m cells are certified by the exact rule CP (ebatch_lib.tau_cp) at n = m.  NDRAW = 1000.
# Expected (analytic): m = 22 certifies eps = 0.10 only with zero sample failures (n_min = 22); eps = 0.05 is vacuous at m = 22 and at m = 40.  These cells are reported as vacuous, not as coverage.
# Truth, coverage, ds_own, anchors: as E2-n (master pin, AUROC anchor).
import argparse, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ebatch_lib as E

MS = (22, 40); NDRAW = 1000


def e3_cell(pop, si, sd, ei, eps, m, ndraw=NDRAW):
    fl = pop["fl"]
    if len(fl) < m: return None
    orac = E.oracle_own(pop, eps); recs = []
    for d in range(ndraw):
        key = [20261103, 81, si, sd, ei, d]; F = np.random.default_rng(key).permutation(len(fl))[:m]; r2 = np.random.default_rng(key + [3])
        cells = np.sort(np.array([pop["negf"][fl[i]][r2.integers(len(pop["negf"][fl[i]]))] for i in F]))
        recs.append(E.evaluate(pop, E.tau_cp(cells, eps), eps, orac))
    return recs


def main():
    ap = E.common_args(argparse.ArgumentParser()); ap.add_argument("--ndraw", type=int, default=NDRAW); a = ap.parse_args()
    L, inp = E.open_inputs(a); master = E.load_master(a.master); out = Path(a.out); assert not out.parent.exists(), "output directory exists: use a new directory"
    seeds = [int(x) for x in a.seeds.split(",")]; by = {}; aur = {s: [] for s in E.SITES}; infeasible = []
    for si, site in enumerate(E.SITES):
        for sd in seeds:
            s, y, g = inp.load_site(sd, site); aur[site].append(E.auroc_bases(s, y)); pop = E.make_pop(s, y, g)
            for ei, eps in enumerate(E.EPSS):
                for m in MS:
                    recs = e3_cell(pop, si, sd, ei, eps, m, a.ndraw)
                    if recs is None: infeasible.append((site, m)); continue
                    by.setdefault((site, m, eps), {})[sd] = recs
    pooled = {}
    for (site, m, eps), rr in by.items(): pooled.setdefault(site, {}).setdefault(str(m), {})[str(eps)] = E.pooled_summary(rr)
    cells = [(s, m, e, v) for s, d in pooled.items() for m, dd in d.items() for e, v in dd.items()]
    reading = dict(E3_a=[dict(site=s, m=int(m), eps=float(e), coverage_nonvacuous=v["coverage_nonvacuous"], vacuous_share=v["vacuous_share"]) for s, m, e, v in cells if v["vacuous_share"] < 1.0 and v["coverage_nonvacuous"] >= 0.90],
                   E3_b=[dict(site=s, m=int(m), eps=float(e), coverage_nonvacuous=v["coverage_nonvacuous"]) for s, m, e, v in cells if v["vacuous_share"] < 1.0 and v["coverage_nonvacuous"] < 0.90],
                   all_vacuous=[dict(site=s, m=int(m), eps=float(e)) for s, m, e, v in cells if v["vacuous_share"] == 1.0])
    basis = {site: E.anchor_auroc(master, site, aur[site]) for site in E.SITES} if seeds == E.SEEDS else None
    E.write_out(out, dict(version="e3_r1", anchor_auroc_basis=basis, spec="E-batch draft3", ndraw=a.ndraw, seeds=seeds, ms=list(MS), pooled=pooled, infeasible=sorted(set(infeasible)), reading=reading, vacuous_convention=E.VAC_CONVENTION, ci_definition=E.CI_DEFINITION, cp_definition="R1 exact (order statistic) on one cell per field, no 200-point grid", inputs_used=dict(inp.used)))


if __name__ == "__main__":
    try: main()
    except E.AnchorError as e: print("ANCHOR MISMATCH (exit 3, no output):", e, file=sys.stderr); sys.exit(3)

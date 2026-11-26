# P5 E5 clinical conversion of the cell-level FPR (spec E-batch draft3 section E5).  LOCAL, read-only, pure arithmetic on stored numbers.  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Inputs: the E1 output json (checked against the SHA256SUMS next to it; N0, Q, MC realised FPR mean and 90th percentile over labeling draws; convention: a vacuous draw has FPR 0), the stored REV-1
#   aggregates (S3 and B1u_c1s, fpr_mean / fpr_p90 at m = 22; convention: manuscript, cannot-certify draws are evaluated at the fallback threshold), the pinned master (SPLIT6 on the development sites:
#   keys ssrcps/cells|SPLIT6|<site>|<eps>|fpr_mean and |fpr_median, same manuscript convention; there is NO fpr_p90 for SPLIT6 in the master or in ssrcps.json, so its p90 is null and the p90 conversion
#   waits for E5b with per-draw thresholds), and per site the mean number of negative cells per field.  That mean uses y and field only, over fields with >= 1 negative (inp.load_site loads the scores too, unused here).
# Every output row carries `convention` ("E1: vacuous->FPR 0" or "manuscript: fallback threshold"); the two conventions must not be mixed in one table.
# Output rows per site, eps, rule:  fp_per_100_fields_mean = fpr_mean x neg_per_field x 100;  fp_per_100_fields_p90 = fpr_p90 x neg_per_field x 100 (p90 over labeling draws, i.e. over the
#   threshold the rule returns, not over fields);  scenario only: P(at least one false positive among N cells) = 1 - (1 - fpr)^N for N in {500, 2000, 5000}, independent-cell approximation.
# There are no slide or patient identifiers: no slide-level statistic is measured.  The cluster-corrected version (needs the per-draw thresholds) is NOT computed here.
# Usage: python -I -B e5_clinical.py <common args without --nihval_dir> --e1 <e1.json> --out <new dir>/e5.json
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ebatch_lib as E

NCELLS = (500, 2000, 5000)


def convert(fpr, neg_per_field):
    return dict(fp_per_100_fields=fpr * neg_per_field * 100.0, p_any_fp={str(n): 1.0 - (1.0 - fpr) ** n for n in NCELLS})


def main():
    ap = E.common_args(argparse.ArgumentParser()); ap.add_argument("--e1", required=True); a = ap.parse_args()
    L, inp = E.open_inputs(a); master = E.load_master(a.master); e1p = Path(a.e1); out = Path(a.out); assert not out.parent.exists(), "output directory exists: use a new directory"
    sums = [l.split()[0] for l in (e1p.parent / "SHA256SUMS").read_text().splitlines() if l.split()[1].lstrip("*") == e1p.name]
    assert sums == [E.sha(e1p)], "E1 json does not match its SHA256SUMS"
    e1 = json.load(open(e1p)); assert e1["version"] == "e1_r1" and e1["seeds"] == E.SEEDS, "E1 json is not the full five-seed run"
    rows = []; npfs = {}
    for site in E.SITES:
        s, y, g = inp.load_site(41, site); pop = E.make_pop(s, y, g); npf = npfs[site] = pop["N"] / len(pop["fl"]); agg = json.load(open(inp.rev1[site]))["agg"]
        for eps in E.EPSS:
            src = {r: e1["pooled"][site][str(eps)][r] for r in ("N0", "Q", "MC")}
            for r, v in src.items():
                rows.append(dict(site=site, eps=eps, rule=r, source="E1", convention="E1: vacuous draw -> FPR 0", neg_per_field=npf, fpr_mean=v["fpr_mean"], fpr_p90=v["fpr_p90"], vacuous_share=v["vacuous_share"],
                                 mean=convert(v["fpr_mean"], npf), p90=convert(v["fpr_p90"], npf)))
            for r in ("S3", "B1u_c1s"):
                v = agg[L.STORED_KEY[site].format(rule=r, eps=eps)]
                rows.append(dict(site=site, eps=eps, rule=r, source="REV-1 stored", convention="manuscript: cannot-certify evaluated at the fallback threshold", neg_per_field=npf, fpr_mean=v["fpr_mean"], fpr_p90=v["fpr_p90"], vacuous_share=v["cannot_certify"],
                                 mean=convert(v["fpr_mean"], npf), p90=convert(v["fpr_p90"], npf)))
    for site in E.SITES:
        npf = npfs[site]
        for eps in E.EPSS:
            k = "ssrcps/cells|SPLIT6|%s|%s|" % (site, eps); fm, fmed, cc = master.get(k + "fpr_mean"), master.get(k + "fpr_median"), master.get(k + "cannot_certify")
            if fm is None: continue
            rows.append(dict(site=site, eps=eps, rule="SPLIT6", source="master " + k + "fpr_mean", convention="manuscript: cannot-certify evaluated at the fallback threshold", neg_per_field=npf, fpr_mean=fm, fpr_median=fmed, fpr_p90=None,
                             vacuous_share=cc, mean=convert(fm, npf), p90=None))
    E.write_out(out, dict(version="e5_r1", spec="E-batch draft3", e1_sha=E.sha(e1p), note="cell-level FPR converted; no slide or patient data; scenario only; not for clinical claims", rows=rows, inputs_used=dict(inp.used)))


if __name__ == "__main__":
    main()

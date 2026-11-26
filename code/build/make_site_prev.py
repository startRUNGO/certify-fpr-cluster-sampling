# P5a Table 4: site prevalence = positives / cells, from the registered dataset counts, as a small JSON for add-only key registration.  WRITTEN, NOT RUN, until 实验管理 reviews it.  Reads ONLY results_master.json (sha pinned).
# Sites with counts in the master: bbbc041_all, mpidb_full, lacuna_B, plasmo, iml (dataset/<site>|positives and |cells).  M5 has no registered positive / negative counts (only m5/n_crops, m5/n_fields), so no M5 prevalence here.
# Cross-check (exit 3 on mismatch): for the four sites that also have rev5_results/agg|a2|<site>|source_threshold|pi_site, that value must equal positives / cells to 1e-12.
# Output: <out>/site_prev.json = {"<site>": prevalence}.  Usage: python3 -B make_site_prev.py --master <tree>/results_master.json --out <new dir>
import argparse, hashlib, json, sys
from pathlib import Path
MASTER_SHA = "f03855aa6994563ab05364ffb329d75a78e0c37c1528867949c4635984ad5301"
SITES = ("bbbc041_all", "mpidb_full", "lacuna_B", "plasmo", "iml"); XCHECK = ("bbbc041_all", "mpidb_full", "lacuna_B", "plasmo")


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--master", required=True); ap.add_argument("--out", required=True); a = ap.parse_args(argv); out = Path(a.out)
    if out.exists(): print("output directory exists", file=sys.stderr); return 2
    if hashlib.sha256(Path(a.master).read_bytes()).hexdigest() != MASTER_SHA: print("GATE (exit 3): master sha differs", file=sys.stderr); return 3
    v = json.loads(Path(a.master).read_text())["values"]; res = {}
    for s in SITES:
        p, c = v[f"dataset/{s}|positives"], v[f"dataset/{s}|cells"]
        # v1.1 (Gemini): bool excluded from the integer check
        if not (type(p) is int and type(c) is int and 0 < p < c): print("GATE (exit 3): counts", s, p, c, file=sys.stderr); return 3
        res[s] = p / c
    for s in XCHECK:
        # v1.1 (Gemini): written as not (<=) so that NaN fails the check
        if not (abs(v[f"rev5_results/agg|a2|{s}|source_threshold|pi_site"] - res[s]) <= 1e-12): print("GATE (exit 3): rev5 pi_site differs", s, file=sys.stderr); return 3
    out.mkdir(parents=True)
    with open(out / "site_prev.json", "x") as fh: json.dump(res, fh, indent=1)
    print(res); return 0


if __name__ == "__main__": sys.exit(main())

# P5a results master r2 (实验管理 2026-11-28): ADD-ONLY builder for the S2 pre-submission robustness results.  LOCAL assembly over pulled, sha-checked summaries; nothing runs.
# Base: the frozen master 131b06d0 (284,763 values).  Every base value, source and file pin is copied unchanged; the builder refuses to overwrite any existing key.
# Keys added (value = number, string or bool; source = {file, path, sha256}):
#   s2A/<cell>|m<m>|<rule>|<eps>|<metric>|<stat>   from the part A summary  (metric in coverage, cannot_certify, sens, ds_eps, ds_own, fpr_mean; stat mean / mcse; coverage also mcse_binomial)
#   s2B/<cell>|<rule>|<eps>|<metric>|<stat>        from the part B simulation summary
#   s2C/<cell>|<variant>|<rule>|<eps>|<metric>|<stat>  from the part C simulation summary
#   s2Bdev/<rule>|<site>|<eps>|<field> and s2Cdev/<variant>|<rule>|<site>|<eps>|<field>  from the development-site summaries (fields n, coverage, fpr_mean, fpr_median,
#     cannot_certify, sens, ds_eps, ds_own, coverage_mcse_seed)
#   <prefix>/n_realisations|<label>|<rule>|<eps>, <prefix>/gate|<name> (the gate statement string), s2A/not_applicable|e12:19|m60 (reason string)
# Per-realisation lists are not copied (they stay in the pulled summaries, pinned by sha).
# Output: a NEW file results_master_p5a_r2_<UTC>.json ('x' mode, allow_nan=False) plus a keys digest TSV of the added keys only.
# Usage: python3 -B add_s2_keys_r2.py --base <results_master.json 131b06d0> --pull <dir with s2A_/s2B_/s2C_/s2devB_/s2devC_ subdirs> --out <dir>
import argparse, hashlib, json, sys, time
from pathlib import Path

BASE_SHA = "131b06d0d0aa8ebeaba6f2d8c3fb102494106ca0182fbdb4d4b40c7913547489"
TAG = "20261128T055603Z"
SRC = {"s2A": (f"s2A_{TAG}/s2_summary.json", "c76470d95fe3ff79273aac8ca20abafc774ec49830a58b0b0f1f40c1e038a06c"),
       "s2B": (f"s2B_{TAG}/s2_summary.json", "6436b5eab1cb0068c2d617fb6593c5579b824297e632e25eed7a6d5bc05a4960"),
       "s2C": (f"s2C_{TAG}/s2_summary.json", "444978b85b91f4ec4650808f11d875853387b573055d8049fd40ad5c39da5863"),
       "s2Bdev": (f"s2devB_{TAG}/s2dev_summary.json", "d287c041d8215b7390eb5504802d608dff293dfb0c8f24e6e051257229817cab"),
       "s2Cdev": (f"s2devC_{TAG}/s2dev_summary.json", None)}
METRICS = ("coverage", "cannot_certify", "sens", "ds_eps", "ds_own", "fpr_mean")
DEV_FIELDS = ("n", "coverage", "fpr_mean", "fpr_median", "cannot_certify", "sens", "ds_eps", "ds_own", "coverage_mcse_seed")


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--base", required=True); ap.add_argument("--pull", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--s2Cdev_sha", required=True)
    a = ap.parse_args(); pull = Path(a.pull); out = Path(a.out); SRC["s2Cdev"] = (SRC["s2Cdev"][0], a.s2Cdev_sha)
    if sha(a.base) != BASE_SHA: sys.exit("GATE base master pin")
    for k, (f, h) in SRC.items():
        if sha(pull / f) != h: sys.exit(f"GATE source pin {k}")
    M = json.load(open(a.base)); V, S = M["values"], M["sources"]; n0 = len(V); added = []

    def add(key, val, k, path):
        if key in V or key in S: sys.exit(f"GATE key exists {key}")
        if not (val is None or isinstance(val, (bool, int, float, str))): sys.exit(f"GATE non-scalar {key}")
        V[key] = val; S[key] = {"file": f"s2/{SRC[k][0]}", "path": path, "sha256": SRC[k][1], "origin": "P5a S2 pre-submission robustness (prereg addendum v1r2 766e7e54)"}; added.append(key)
    for k in ("s2A", "s2B", "s2C"):
        d = json.load(open(pull / SRC[k][0]))
        for g, txt in d["gates"].items(): add(f"{k}/gate|{g}", txt, k, ["gates", g])
        for lab, b in d["results"].items():
            add(f"{k}/n_realisations|{lab}", b["n_realisations"], k, ["results", lab, "n_realisations"])
            for mn in METRICS:
                add(f"{k}/{lab}|{mn}|mean", b[mn]["mean"], k, ["results", lab, mn, "mean"])
                add(f"{k}/{lab}|{mn}|mcse", b[mn]["mcse"], k, ["results", lab, mn, "mcse"])
            add(f"{k}/{lab}|coverage|mcse_binomial", b["coverage"]["mcse_binomial"], k, ["results", lab, "coverage", "mcse_binomial"])
        for na in d.get("not_applicable", []): add(f"{k}/not_applicable|{na['cell']}|m{na['m']}", na["reason"], k, ["not_applicable"])
    for k in ("s2Bdev", "s2Cdev"):
        d = json.load(open(pull / SRC[k][0]))
        for g, txt in d["gates"].items(): add(f"{k}/gate|{g}", txt, k, ["gates", g])
        for lab, b in d["results"].items():
            for f in DEV_FIELDS: add(f"{k}/{lab}|{f}", b[f], k, ["results", lab, f])
    if "r2_added" in M: sys.exit("GATE base already has r2_added")
    M["n_values"] = len(V)
    M["r2_added"] = dict(builder="add_s2_keys_r2.py", builder_sha256=sha(__file__), base_sha256=BASE_SHA, n_base=n0, n_added=len(added), sources={k: {"file": f, "sha256": h} for k, (f, h) in SRC.items()})
    ts = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()); out.mkdir(exist_ok=True)
    mp = out / f"results_master_p5a_r2_{ts}.json"
    with open(mp, "x", encoding="utf-8") as fh: json.dump(M, fh, ensure_ascii=False, allow_nan=False)
    dp = out / f"results_master_p5a_r2_added_keys_{ts}.tsv"
    with open(dp, "x", encoding="utf-8") as fh:
        fh.write("key\tvalue\tsource\tpath\n")
        for key in added: fh.write(f"{key}\t{json.dumps(V[key], ensure_ascii=False)}\t{S[key]['file']} ({S[key]['sha256'][:8]})\t{json.dumps(S[key]['path'], ensure_ascii=False)}\n")
    print(n0, "base +", len(added), "added =", len(V)); print(sha(mp), mp); print(sha(dp), dp)


if __name__ == "__main__":
    main()

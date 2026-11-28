# P5a results master r3 (实验管理 2026-11-28): ADD-ONLY builder for the S3 addendum (prereg addendum v1r2 713ebe2b). LOCAL assembly over pulled, sha-checked summaries; nothing runs.
# Base: master r2 7e56c5c9 (292,241 values). Every base value, source and pin is copied unchanged. The builder refuses to overwrite any existing key.
# Keys added: every scalar leaf (number, string, bool, null) of
#   s3sim_summary.json: results/<label>/... -> s3sim/<label>|<path...>, gates/<g> -> s3sim/gate|<g>, settings/<scalar> -> s3sim/settings|<name>, calib/<f> -> s3sim/calib|<f>
#   s3dev_summary.json: results/<label>/... -> s3dev/<label>|<path...>, gates/<g> -> s3dev/gate|<g>
# Lists (per_realisation, coverage_per_seed, settings/realisations, variants) are not copied. They stay in the pulled summaries, pinned by sha.
# Output: a NEW file results_master_p5a_r3_<UTC>.json ('x' mode, allow_nan=False) plus a TSV digest of the added keys only.
import argparse, hashlib, json, sys, time
from pathlib import Path
BASE_SHA = "7e56c5c99bdcda85109aa9d56790ab52a555af3af13264786db86777c36c86d8"
SRC = {"s3sim": ("s3sim_20261128T094156Z/s3sim_summary.json", "a21b37ff"), "s3dev": ("s3dev_20261128T094156Z/s3dev_summary.json", "3513614d0cd35cfe2f9a088bffde4871cf0fef0c0950e90733b8b636cb6eeb24")}
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--base", required=True); ap.add_argument("--pull", required=True); ap.add_argument("--out", required=True); ap.add_argument("--s3sim_sha", required=True)
    a = ap.parse_args(); pull = Path(a.pull); out = Path(a.out)
    if not a.s3sim_sha.startswith(SRC["s3sim"][1]) or len(a.s3sim_sha) != 64: sys.exit("GATE s3sim sha argument")
    SRC["s3sim"] = (SRC["s3sim"][0], a.s3sim_sha)
    if sha(a.base) != BASE_SHA: sys.exit("GATE base master pin")
    for k, (f, h) in SRC.items():
        if sha(pull / f) != h: sys.exit(f"GATE source pin {k}")
    M = json.load(open(a.base)); V, S = M["values"], M["sources"]; n0 = len(V); added = []
    def add(key, val, k, path):
        if key in V or key in S: sys.exit(f"GATE key exists {key}")
        if not (val is None or isinstance(val, (bool, int, float, str))): sys.exit(f"GATE non-scalar {key}")
        V[key] = val; S[key] = {"file": f"s3/{SRC[k][0]}", "path": path, "sha256": SRC[k][1], "origin": "P5a S3 addendum (prereg addendum v1r2 713ebe2b)"}; added.append(key)
    def walk(o, key, k, path):
        if isinstance(o, dict):
            for c, v in o.items(): walk(v, f"{key}|{c}", k, path + [c])
        elif isinstance(o, list): return
        else: add(key, o, k, path)
    for k in ("s3sim", "s3dev"):
        d = json.load(open(pull / SRC[k][0]))
        for g, txt in d["gates"].items(): add(f"{k}/gate|{g}", txt, k, ["gates", g])
        for lab, b in d["results"].items(): walk(b, f"{k}/{lab}", k, ["results", lab])
        if k == "s3sim":
            for n, v in d["settings"].items():
                if not isinstance(v, (list, dict)): add(f"s3sim/settings|{n}", v, k, ["settings", n])
            for n, v in d["calib"].items(): add(f"s3sim/calib|{n}", v, k, ["calib", n])
            if d["exceptions"]: sys.exit("GATE s3sim exceptions not empty, decide by hand")
    if "r3_added" in M: sys.exit("GATE base already has r3_added")
    M["n_values"] = len(V)
    M["r3_added"] = dict(builder="add_s3_keys_r3.py", builder_sha256=sha(__file__), base_sha256=BASE_SHA, n_base=n0, n_added=len(added), sources={k: {"file": f, "sha256": h} for k, (f, h) in SRC.items()})
    ts = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()); out.mkdir(exist_ok=True)
    mp = out / f"results_master_p5a_r3_{ts}.json"
    with open(mp, "x", encoding="utf-8") as fh: json.dump(M, fh, ensure_ascii=False, allow_nan=False)
    dp = out / f"results_master_p5a_r3_added_keys_{ts}.tsv"
    with open(dp, "x", encoding="utf-8") as fh:
        fh.write("key\tvalue\tsource\tpath\n")
        for key in added: fh.write(f"{key}\t{json.dumps(V[key], ensure_ascii=False)}\t{S[key]['file']} ({S[key]['sha256'][:8]})\t{json.dumps(S[key]['path'], ensure_ascii=False)}\n")
    B = json.load(open(a.base))["values"]
    bad = [x for x in B if x not in V or V[x] != B[x]]
    print("add-only check: base keys", len(B), "changed or missing", len(bad))
    print(n0, "base +", len(added), "added =", len(V)); print(sha(mp), mp); print(sha(dp), dp)
if __name__ == "__main__":
    main()

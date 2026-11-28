# P5a results master r3b (实验管理 2026-11-28): ADD-ONLY calibration keys of the S3 site-calibrated generators. LOCAL assembly, nothing runs.
# Base: master r3a 52de2781. Adds every scalar field of s3calib.json (ad11c571) calibration/<site>|<seed> as s3calib/<site>|<seed>|<field>.
# Fields: F, mu0, tau_b, s_w, prevalence, n_fields_dropped, n_pos_dropped, n_res, n_off. Lists and dicts inside an entry, score_files_sha256 (VM paths) and the note are not copied.
import hashlib, json, sys, time
from pathlib import Path
BASE_SHA = "52de27812200a9827053a36bf52ae696fa030568bda1e0a97a5ed1f6d133c5b1"
CAL_SHA = "ad11c5712a1915b1aa11c4ae45bdbac4de7721fb8a0795bd89ef1e0623fcca9d"
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
base, cal, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
if sha(base) != BASE_SHA: sys.exit("GATE base pin")
if sha(cal) != CAL_SHA: sys.exit("GATE calib pin")
M = json.load(open(base)); V, S = M["values"], M["sources"]; n0 = len(V); added = []
C = json.load(open(cal))["calibration"]
for lab, e in C.items():
    for f, v in e.items():
        if isinstance(v, (list, dict)): continue
        k = f"s3calib/{lab}|{f}"
        if k in V or k in S: sys.exit(f"GATE key exists {k}")
        if not (v is None or isinstance(v, (bool, int, float, str))): sys.exit(f"GATE non-scalar {k}")
        if isinstance(v, str) and ("/" in v or "<USER>" in v): sys.exit(f"GATE path-like string {k}")
        V[k] = v; S[k] = {"file": "s3/s3calib_20261128T094156Z/s3calib.json", "path": ["calibration", lab, f], "sha256": CAL_SHA, "origin": "P5a S3 addendum calibration (prereg addendum v1r2 713ebe2b)"}; added.append(k)
if "r3b_added" in M: sys.exit("GATE r3b exists")
M["n_values"] = len(V); M["r3b_added"] = dict(builder="add_s3calib_keys_r3b.py", builder_sha256=sha(__file__), base_sha256=BASE_SHA, n_base=n0, n_added=len(added), source={"file": "s3calib.json", "sha256": CAL_SHA})
ts = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()); out.mkdir(exist_ok=True)
mp = out / f"results_master_p5a_r3b_{ts}.json"
with open(mp, "x", encoding="utf-8") as fh: json.dump(M, fh, ensure_ascii=False, allow_nan=False)
dp = out / f"results_master_p5a_r3b_added_keys_{ts}.tsv"
with open(dp, "x", encoding="utf-8") as fh:
    fh.write("key\tvalue\n")
    for k in added: fh.write(f"{k}\t{json.dumps(V[k])}\n")
B = json.load(open(base))["values"]
print("add-only check: changed or missing", sum(1 for x in B if x not in V or V[x] != B[x]))
print(n0, "base +", len(added), "added =", len(V)); print(sha(mp), mp); print(sha(dp), dp)

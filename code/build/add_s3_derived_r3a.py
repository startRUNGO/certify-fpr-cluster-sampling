# P5a results master r3a (实验管理 2026-11-28): ADD-ONLY derived keys for the S3 text of v6.32. LOCAL assembly, nothing runs.
# Base: master r3 ee23a22a. Every base value is copied unchanged and no existing key is overwritten.
# Derived keys (source = {"derived": definition, "from_keys": [...]}):
#   s3sim_derived/min_coverage|0.1  = min of s3sim/<gen>|<site>|<rule>|0.1|coverage|mean over gen {A,B} x site {bbbc041_all, lacuna_B} x rule {S3, SPLIT6}
#   s3sim_derived/n_checks|0.1      = number of those cells (8)
#   s3sim_derived/n_pass|0.1        = number of those cells with mean >= 0.90 - 2 mcse
#   s3dev_derived/split6_min_coverage|0.1 = min of s3dev/<variant>|SPLIT6|<site>|0.1|coverage over 7 variants x 3 sites
#   s3dev_derived/n_cells|0.1       = number of those cells (21)
#   s3dev_derived/n_pass|0.1        = number of those cells with coverage >= 0.90 - 2 coverage_mcse_seed
import hashlib, json, sys, time
from pathlib import Path
BASE_SHA = "ee23a22a343a7a62a3aae5af1dd1f22f489987fc52eca0d5efbd127d4303f684"
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
base, out = Path(sys.argv[1]), Path(sys.argv[2])
if sha(base) != BASE_SHA: sys.exit("GATE base pin")
M = json.load(open(base)); V, S = M["values"], M["sources"]; n0 = len(V); added = []
def add(k, v, d, fk):
    if k in V or k in S: sys.exit(f"GATE key exists {k}")
    V[k] = v; S[k] = {"derived": d, "from_keys": fk}; added.append(k)
cells = [f"s3sim/{g}|{s}|{r}|0.1|coverage|" for g in "AB" for s in ("bbbc041_all", "lacuna_B") for r in ("S3", "SPLIT6")]
for c in cells:
    if c + "mean" not in V or c + "mcse" not in V: sys.exit(f"GATE missing {c}")
mn = min(V[c + "mean"] for c in cells); npass = sum(V[c + "mean"] >= 0.90 - 2 * V[c + "mcse"] for c in cells)
fk = [c + "mean" for c in cells] + [c + "mcse" for c in cells]
add("s3sim_derived/min_coverage|0.1", mn, "min over generator x site x rule of the mean coverage at eps 0.10", [c + "mean" for c in cells])
add("s3sim_derived/n_checks|0.1", len(cells), "number of generator x site x rule cells at eps 0.10", [c + "mean" for c in cells])
add("s3sim_derived/n_pass|0.1", int(npass), "cells with mean coverage >= 0.90 - 2 MCSE at eps 0.10 (prediction 1)", fk)
dv = [f"s3dev/{v}|SPLIT6|{s}|0.1|" for v in ("default", "ns10", "ns40", "rlow0.6", "rlow0.8", "G100", "G400") for s in ("bbbc041_all", "mpidb_full", "lacuna_B")]
for c in dv:
    if c + "coverage" not in V or c + "coverage_mcse_seed" not in V: sys.exit(f"GATE missing {c}")
add("s3dev_derived/split6_min_coverage|0.1", min(V[c + "coverage"] for c in dv), "min of Split-Select coverage at eps 0.10 over 7 tuning variants x 3 development sites", [c + "coverage" for c in dv])
add("s3dev_derived/n_cells|0.1", len(dv), "number of tuning variant x development site cells", [c + "coverage" for c in dv])
add("s3dev_derived/n_pass|0.1", int(sum(V[c + "coverage"] >= 0.90 - 2 * V[c + "coverage_mcse_seed"] for c in dv)), "cells with coverage >= 0.90 - 2 seed MCSE (prediction 2)", [c + x for c in dv for x in ("coverage", "coverage_mcse_seed")])
if "r3a_added" in M: sys.exit("GATE r3a exists")
M["n_values"] = len(V); M["r3a_added"] = dict(builder="add_s3_derived_r3a.py", builder_sha256=sha(__file__), base_sha256=BASE_SHA, n_base=n0, n_added=len(added))
ts = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()); out.mkdir(exist_ok=True)
mp = out / f"results_master_p5a_r3a_{ts}.json"
with open(mp, "x", encoding="utf-8") as fh: json.dump(M, fh, ensure_ascii=False, allow_nan=False)
dp = out / f"results_master_p5a_r3a_added_keys_{ts}.tsv"
with open(dp, "x", encoding="utf-8") as fh:
    fh.write("key\tvalue\tdefinition\n")
    for k in added: fh.write(f"{k}\t{json.dumps(V[k])}\t{S[k]['derived']}\n")
B = json.load(open(base))["values"]
print("add-only check: changed or missing", sum(1 for x in B if x not in V or V[x] != B[x]))
for k in added: print(k, V[k])
print(sha(mp), mp); print(sha(dp), dp)

# P5a ADD-ONLY results_master builder r3: largest all50 MCSE of Delta Sens_eps and of the cannot-certify rate over the 340 sim50 settings (实验管理 2026-11-26, draft 17 item D6).
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  No input file: every new key is computed from sim50/* keys already in the master.
# Same mechanics as add_sim50_derived2_keys.py (reviewed): write targets checked absent first; indent=1; every old key, source and files_sha256 entry unchanged;
#   key count checked; every new key recomputed from the written file and compared.
# Keys (none may exist before the run):
#   sim50_derived/max_mcse_ds_eps|value, |argmax and sim50_derived/max_mcse_cannot_certify|value, |argmax ('<d>|<ci>|<rule>|<eps>') over the 340 settings
#   (34 cells x 5 rules x 2 eps) of sim50/<d>|<ci>|<rule>|<eps>|all50|<stat>|mcse, as sim50_derived/max_mcse|value is for coverage.
# Sources: {'derived': <rule in words>, 'from_keys': <keys used>}.
# Gates (exit 3, nothing written): master sha 5dee0985; every new key absent; every input key present and finite; 340 values per stat; all write targets absent.
# Backup: results_master_before_sim50_derived3_5dee0985.json ('x').
# Usage: python3 -E -B add_sim50_derived3_keys.py --tree <tree_bj_v1_...>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "5dee0985a20b96e9995e6ae5f2adaeb98ff2f11bf7dda491dc0a2895c61200a6"
PRE = "sim50_derived"
RULES_ALL = ("S3", "SPLIT6", "CB_c1s", "GEE_c1s", "BB_c1s")
EPS = ("0.1", "0.05")
CELLS = {"e10": range(10), "e12": range(24)}
STATS = ("ds_eps", "cannot_certify")


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def derive(V):
    """Returns {key: (value, source)} computed from V, or raises ValueError."""
    def num(k):
        if k not in V or not isinstance(V[k], (int, float)) or isinstance(V[k], bool) or not math.isfinite(V[k]): raise ValueError(f"missing or non-finite {k}")
        return V[k]
    new = {}
    for st in STATS:
        rows = [(num(f"sim50/{d}|{ci}|{r}|{e}|all50|{st}|mcse"), f"{d}|{ci}|{r}|{e}", f"sim50/{d}|{ci}|{r}|{e}|all50|{st}|mcse") for d, cs in CELLS.items() for ci in cs for r in RULES_ALL for e in EPS]
        if len(rows) != 340: raise ValueError(f"{st}: {len(rows)} MCSE values, expected 340")
        mx = max(rows, key=lambda x: x[0]); used = sorted(k for _, _, k in rows)
        w = f"over the 340 settings (34 cells x 5 rules x 2 eps) of the all50 {st} MCSE"
        for k, v, words in ((f"max_mcse_{st}|value", mx[0], "maximum "), (f"max_mcse_{st}|argmax", mx[1], "setting of the maximum ")):
            k = f"{PRE}/{k}"
            if k in new: raise ValueError(f"duplicate {k}")
            new[k] = (v, dict(derived=words + w, from_keys=used))
    return new


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); a = ap.parse_args()
    tree = Path(a.tree); rm = tree / "results_master.json"
    try:
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not 5dee0985")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        new = derive(vals)
        clash = [k for k in new if k in vals]
        if clash: raise ValueError(f"keys exist already: {clash}")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    bak = tree / "results_master_before_sim50_derived3_5dee0985.json"; tmp = tree / "results_master.json.sim50_derived3_tmp"
    existing = [str(x) for x in (bak, tmp) if x.exists()]
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}; old_fsha = dict(fsha)
    with open(rm, "rb") as s_, open(bak, "xb") as d_: d_.write(s_.read())
    for k, (v, s) in new.items(): vals[k] = v; srcs[k] = s
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp))
    try: again = derive({k: v for k, v in chk["values"].items() if k not in new})
    except ValueError as e: print("GATE (exit 3): recompute failed", e, file=sys.stderr); return 3
    for k, (v, _) in new.items():
        if chk["values"].get(k) != v or again[k][0] != v: print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    for k, v in old_fsha.items():
        if chk["files_sha256"].get(k) != v: print("GATE (exit 3): an old files_sha256 entry changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new): print("GATE (exit 3): key count mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    print(json.dumps(dict(added=len(new), n_values=d["n_values"], new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

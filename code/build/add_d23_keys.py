# P5a ADD-ONLY results_master builder for draft 23 item O2 (实验管理 2026-11-26): the cost of certification by simulation block, for the new main-text table in Section 4.6.
#   WRITTEN, NOT RUN, until 实验管理 reviews it.  No input file: every new key is the mean of sim50 / sizefpr keys already in the master.
# Same mechanics as add_sim50_derived3_keys.py (reviewed): write targets checked absent first; indent=1; every old key, source and files_sha256 entry unchanged;
#   key count checked; every new key recomputed from the written file and compared.
# Blocks (scenario indices; e12 membership gated by the scenario names in the master, e12/cells|<i>|cell|name):
#   single = sim50/e10|0..9 (the single-factor design); rare = e12 1..6 (names A_*); size = e12 7..10 (B_*); slide = e12 11..16 (C_*); twostage = e12 17..20 (C2_*);
#   smallsite = e12 21..23 (D_*); sizefpr = sizefpr/sc|0..3.  The stress base scenario e12 0 (name 'base') belongs to no block and is left out.
# Keys (none may exist before the run), at eps = 0.1, rules S3 (DE-Grid) and SPLIT6 (Split-Select):
#   sim50_derived/block|<block>|<rule>|0.1|<stat>|mean = mean over the block's scenarios of <prefix>|<i>|<rule>|0.1|all50|<stat>|mean, stat in {cannot_certify, ds_eps};
#   sim50_derived/block|<block>|n_scenarios = number of scenarios in the block.  7 x 2 x 2 + 7 = 35 keys.
# Sources: {'derived': <rule in words>, 'from_keys': <keys used>}.
# Gates (exit 3, nothing written): master sha 0f4068e4; every new key absent; every input key present and finite; block membership by name; all write targets absent.
# Backup: results_master_before_d23_0f4068e4.json ('x').
# Usage: python3 -E -B add_d23_keys.py --tree <tree_bj_v1_...>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "0f4068e4035bcae1254ebc7fae4108867f1bb83143da9d5ac724227bc32c7d82"
PRE = "sim50_derived/block"
BLOCKS = {"single": ("sim50/e10", range(0, 10), None), "rare": ("sim50/e12", range(1, 7), "A_"), "size": ("sim50/e12", range(7, 11), "B_"),
          "slide": ("sim50/e12", range(11, 17), "C_"), "twostage": ("sim50/e12", range(17, 21), "C2_"), "smallsite": ("sim50/e12", range(21, 24), "D_"),
          "sizefpr": ("sizefpr/sc", range(0, 4), None)}
RULES = ("S3", "SPLIT6")
STATS = ("cannot_certify", "ds_eps")


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return hashlib.sha256(fh.read()).hexdigest()


def derive(V):
    """Returns {key: (value, source)} computed from V, or raises ValueError."""
    def num(k):
        if k not in V or not isinstance(V[k], (int, float)) or isinstance(V[k], bool) or not math.isfinite(V[k]): raise ValueError(f"missing or non-finite {k}")
        return V[k]
    new = {}
    for b, (pre, idx, name) in BLOCKS.items():
        if name is not None:
            for i in idx:
                nm = V.get(f"e12/cells|{i}|cell|name")
                if not isinstance(nm, str) or not nm.startswith(name) or (name == "C_" and nm.startswith("C2_")): raise ValueError(f"block {b}: e12 scenario {i} is {nm}")
        if V.get("e12/cells|0|cell|name") != "base": raise ValueError("e12 scenario 0 is not the base")
        new[f"{PRE}|{b}|n_scenarios"] = (len(idx), dict(derived=f"number of scenarios in block {b} ({pre} {idx.start}..{idx.stop - 1})", from_keys=[]))
        for r in RULES:
            for st in STATS:
                ks = [f"{pre}|{i}|{r}|0.1|all50|{st}|mean" for i in idx]
                vals = [num(k) for k in ks]
                new[f"{PRE}|{b}|{r}|0.1|{st}|mean"] = (sum(vals) / len(vals), dict(derived=f"mean over the {len(idx)} scenarios of block {b} of the all50 {st} mean, rule {r}, eps 0.1", from_keys=ks))
    if len(new) != 35: raise ValueError(f"{len(new)} new keys, expected 35")
    return new


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); a = ap.parse_args()
    tree = Path(a.tree); rm = tree / "results_master.json"
    try:
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master is not 0f4068e4")
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        new = derive(vals)
        clash = [k for k in new if k in vals]
        if clash: raise ValueError(f"keys exist already: {clash}")
    except (ValueError, KeyError, OSError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    bak = tree / "results_master_before_d23_0f4068e4.json"; tmp = tree / "results_master.json.d23_tmp"
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

# P5a ADD-ONLY results_master builder: size-FPR check (sizefpr_check.py 7d4727be, sizecorr_gen.py 3c3a7b1e, VM run sizefpr_20261125T214523Z).  WRITTEN, NOT RUN, until 实验管理 reviews it.
# Same write mechanics as add_sim50_keys.py and add_sim50_derived2_keys.py (reviewed): write targets checked absent first; backup 'x'; tmp master re-read; every new key re-read
#   from its source; every old key, source and files_sha256 entry unchanged; key count checked; indent=1.
# Input: the run directory pulled via the scratchpad (sizefpr_summary.json + SHA256SUMS).  The sha of sizefpr_summary.json is asserted == SUMMARY_SHA on the raw bytes before
#   the file is parsed.
# Keys, prefix sizefpr/ (results key as in the file, e.g. sc|1|S3|0.05; cells 0..3, rules S3 SPLIT6 CB_c1s GEE_c1s BB_c1s, eps 0.1 0.05):
#   sizefpr/<rk>|all50|n_realisations, sizefpr/<rk>|all50|<stat>|mean and |mcse for stat in coverage cannot_certify sens ds_eps ds_own fpr_mean,
#     sizefpr/<rk>|all50|coverage|mcse_binomial, sizefpr/<rk>|all50|coverage|ci90_seedcluster|0 and |1  (per-realisation values and seed lists are not registered);
#   sizefpr/cell|<i>|size_sigma and sizefpr/cell|<i>|rho_c_nominal (the generator parameters; rho_c is the NOMINAL value and is never cited as achieved);
#   sizefpr/ndraw, sizefpr/m, sizefpr/n_seeds;
#   sizefpr/corr_realised|<i>|mean (from the file) and the derived sizefpr/corr_realised|<i>|min, |max (over the 50 per-realisation values of the file).
# Derived keys, prefix sizefpr_derived/ (sources: {'derived': <rule in words>, 'from_keys': [...]} or {'file': ..., 'derived': ...}):
#   <rule>|n_settings, n_hold, n_between, n_fail, min_coverage, argmin, min_mcse for rule in S3, SPLIT6 over the 8 size-FPR settings (4 cells x 2 eps), delta = protocol/delta;
#   <rule>|n_hold_56 = sim50_derived/e12|<rule>|n_hold + sizefpr_derived/<rule>|n_hold and <rule>|n_settings_56 = sim50_derived/e12|<rule>|n_settings + 8, for S3 and SPLIT6.
#   The old keys sim50_derived/e12|<rule>|n_hold and |n_settings are read only and stay unchanged.
# Gates (exit 3, nothing written): summary sha == SUMMARY_SHA (before parsing); SHA256SUMS of the run dir; master sha == MASTER_SHA; version sizefpr_r1; script_sha256 == SCRIPT_SHA;
#   pins sizecorr_gen 3c3a7b1e and sim50_extend ea3ea087; anchors g, r, c start with "passed"; platt_fallback empty; settings cells == CELLS, seeds 41..90, ndraw 200, m 22;
#   result keys == 4 cells x 5 rules x 2 eps = 40; n_realisations 50 and |mean - mean(per_realisation)| <= 1e-12 for every stat registered; ci90 lo <= coverage mean <= hi;
#   s3_hold_1_minus_delta of the file equals the recomputed hold (coverage >= 1 - delta) for each S3 setting and s3_holds_everywhere equals all of them;
#   50 correlation values per cell, mean equal to their mean; protocol/delta == 0.1; old count keys present and n_settings == 48;
#   no existing key starts with sizefpr/ or sizefpr_derived/; no non-finite value; all write targets absent.
# Backup: results_master_before_sizefpr_<first 8 of MASTER_SHA>.json ('x').  Copy: <tree>/results_src/sizefpr/<run dir name>/sizefpr_summary.json (+ files_sha256).
# Usage: python3 -E -B add_sizefpr_keys.py --tree <tree_bj_v1_...> --src <run dir>
import argparse, hashlib, json, math, sys
from pathlib import Path
MASTER_SHA = "43717cbe47d37a402cc81258ad350cd4c0482680fab3a13575b7f3c4901f142d"
SUMMARY_SHA = "9bc82641fb6c6468fbe0a5f24cea2e08ca54ead394fbbe16f4b4ff8161a24140"
SCRIPT_SHA = "7d4727bee33654da23384fc8712cf452530421ef9c67f4c289efbf1965fe5bf3"
PIN_GEN = "3c3a7b1edcc6e9ff003321a890420a73d5e58660e49098ff3b8d07152bdf8220"
PIN_SIM50 = "ea3ea087620e8bcda8d08a95a80db6897fd26622ede8c5842262699eb8b602a6"
RUN_DIR = "sizefpr_20261125T214523Z"
PRE, PRED, FNAME, DST = "sizefpr", "sizefpr_derived", "sizefpr_summary.json", "results_src/sizefpr"
CELLS = [(0, 1.0, 0.4), (1, 1.0, 0.8), (2, 0.5, 0.8), (3, 1.0, -0.8)]
RULES = ("S3", "SPLIT6", "CB_c1s", "GEE_c1s", "BB_c1s")
COUNT_RULES = ("S3", "SPLIT6")
EPS = ("0.1", "0.05")
STATS = ("coverage", "cannot_certify", "sens", "ds_eps", "ds_own", "fpr_mean")
TOL = 1e-12


def sha_bytes(b): return hashlib.sha256(b).hexdigest()


def sha(p):
    if not Path(p).is_file(): return None
    with open(p, "rb") as fh: return sha_bytes(fh.read())


def finite(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def file_keys(o):
    """Returns {key: (value, path)} read from the summary object o, or raises ValueError."""
    new = {}

    def put(k, path):
        v = o
        for p in path: v = v[p]
        if not finite(v): raise ValueError(f"non-finite or non-numeric {path}")
        k = f"{PRE}/{k}"
        if k in new: raise ValueError(f"duplicate {k}")
        new[k] = (v, list(path))
    for rk in o["results"]:
        g = o["results"][rk]["all50"]; b = ["results", rk, "all50"]
        put(f"{rk}|all50|n_realisations", b + ["n_realisations"])
        for st in STATS:
            put(f"{rk}|all50|{st}|mean", b + [st, "mean"]); put(f"{rk}|all50|{st}|mcse", b + [st, "mcse"])
        put(f"{rk}|all50|coverage|mcse_binomial", b + ["coverage", "mcse_binomial"])
        put(f"{rk}|all50|coverage|ci90_seedcluster|0", b + ["coverage", "ci90_seedcluster", 0])
        put(f"{rk}|all50|coverage|ci90_seedcluster|1", b + ["coverage", "ci90_seedcluster", 1])
    for i, c in enumerate(o["settings"]["cells"]):
        put(f"cell|{i}|size_sigma", ["settings", "cells", i, "size_sigma"]); put(f"cell|{i}|rho_c_nominal", ["settings", "cells", i, "rho_c"])
    put("ndraw", ["settings", "ndraw"]); put("m", ["settings", "m"])
    for i in o["size_fpr_correlation"]: put(f"corr_realised|{i}|mean", ["size_fpr_correlation", i, "mean"])
    return new


def derived_keys(o, V):
    """Returns {key: (value, source)} computed from o and the old master values V, or raises ValueError."""
    delta = V["protocol/delta"]
    if delta != 0.1: raise ValueError("protocol/delta differs from 0.1")
    lo1, lo2 = 1 - delta, 1 - 2 * delta
    rel_file = f"{DST}/{RUN_DIR}/{FNAME}"; new = {}

    def put(k, v, src):
        if not (finite(v) or isinstance(v, str)): raise ValueError(f"bad derived value {k}")
        k = f"{PRED}/{k}" if not k.startswith(PRE + "/") else k
        if k in new: raise ValueError(f"duplicate {k}")
        new[k] = (v, src)
    put(f"{PRE}/n_seeds", len(o["settings"]["seeds"]), dict(file=rel_file, derived="number of entries of settings.seeds"))
    for i, v in o["size_fpr_correlation"].items():
        pr = v["per_realisation"]
        put(f"{PRE}/corr_realised|{i}|min", min(pr), dict(file=rel_file, derived=f"minimum over the 50 realisations of the realised size-FPR correlation, cell {i}"))
        put(f"{PRE}/corr_realised|{i}|max", max(pr), dict(file=rel_file, derived=f"maximum over the 50 realisations of the realised size-FPR correlation, cell {i}"))
    for r in COUNT_RULES:
        rows = []
        for ci, _, _ in CELLS:
            for e in EPS:
                b = f"{PRE}/sc|{ci}|{r}|{e}|all50|coverage|"
                g = o["results"][f"sc|{ci}|{r}|{e}"]["all50"]["coverage"]
                rows.append((g["mean"], g["mcse"], f"sc|{ci}|{e}", [b + "mean", b + "mcse"]))
        used = sorted(k for *_, ks in rows for k in ks)
        cmin, mmin, arg, _ = min(rows, key=lambda x: x[0])
        nh = sum(c >= lo1 for c, *_ in rows); nb = sum(lo2 <= c < lo1 for c, *_ in rows); nf = sum(c < lo2 for c, *_ in rows)
        if nh + nb + nf != len(rows) or len(rows) != 8: raise ValueError(f"{r}: counts do not add up")
        w = f"over the 8 size-FPR settings (4 cells x eps 0.1, 0.05) of all50 coverage of rule {r}, delta = protocol/delta"
        src = dict(derived=None, from_keys=used + ["protocol/delta"])
        for k, v, words in ((f"{r}|n_settings", len(rows), "number of settings "), (f"{r}|n_hold", nh, "settings with coverage >= 1 - delta "),
                            (f"{r}|n_between", nb, "settings with 1 - 2 delta <= coverage < 1 - delta "), (f"{r}|n_fail", nf, "settings with coverage < 1 - 2 delta "),
                            (f"{r}|min_coverage", cmin, "minimum "), (f"{r}|argmin", arg, "setting of the minimum "), (f"{r}|min_mcse", mmin, "all50 coverage MCSE at the minimum ")):
            put(k, v, dict(src, derived=words + w))
        oh, on = f"sim50_derived/e12|{r}|n_hold", f"sim50_derived/e12|{r}|n_settings"
        if not (finite(V.get(oh)) and finite(V.get(on))) or V[on] != 48: raise ValueError(f"old count keys of {r} missing or n_settings differs from 48")
        put(f"{r}|n_hold_56", V[oh] + nh, dict(derived=f"{oh} + {PRED}/{r}|n_hold (the 48 E12 stress settings and the 8 size-FPR settings)", from_keys=[oh, f"{PRED}/{r}|n_hold"]))
        put(f"{r}|n_settings_56", V[on] + len(rows), dict(derived=f"{on} + {PRED}/{r}|n_settings", from_keys=[on, f"{PRED}/{r}|n_settings"]))
        if V[on] + len(rows) != 56: raise ValueError("combined setting count differs from 56")
    return new


def check_object(o):
    if o["version"] != "sizefpr_r1" or o["script_sha256"] != SCRIPT_SHA: raise ValueError("version or script_sha256")
    if o["pins"].get("sizecorr_gen") != PIN_GEN or o["pins"].get("sim50_extend") != PIN_SIM50: raise ValueError("pins")
    if not all(o["anchors"][a].startswith("passed") for a in ("g", "r", "c")): raise ValueError(f"anchors {o['anchors']}")
    if o["platt_fallback"]: raise ValueError("platt_fallback not empty")
    s = o["settings"]
    if [(c["index"], c["size_sigma"], c["rho_c"]) for c in s["cells"]] != CELLS: raise ValueError("cells differ")
    if s["seeds"] != list(range(41, 91)) or s["ndraw"] != 200 or s["m"] != 22: raise ValueError("seeds / ndraw / m")
    want = {f"sc|{c}|{r}|{e}" for c, _, _ in CELLS for r in RULES for e in EPS}
    if set(o["results"]) != want: raise ValueError(f"result keys differ: missing {sorted(want - set(o['results']))[:5]}, extra {sorted(set(o['results']) - want)[:5]}")
    for rk, v in o["results"].items():
        g = v["all50"]
        if g["n_realisations"] != 50 or g["seeds"] != list(range(41, 91)): raise ValueError(f"realisations {rk}")
        for st in STATS:
            pr = g[st]["per_realisation"]
            if len(pr) != 50 or abs(g[st]["mean"] - sum(pr) / 50) > TOL: raise ValueError(f"mean differs from its realisations {rk} {st}")
        lo, hi = g["coverage"]["ci90_seedcluster"]
        if not lo <= g["coverage"]["mean"] <= hi: raise ValueError(f"coverage outside its interval {rk}")
    hold = {f"sc|{c}|S3|{e}": o["results"][f"sc|{c}|S3|{e}"]["all50"]["coverage"]["mean"] >= 0.9 for c, _, _ in CELLS for e in EPS}
    if o["s3_hold_1_minus_delta"] != hold or o["s3_holds_everywhere"] != all(hold.values()): raise ValueError("hold dict differs from the recomputed hold")
    if set(o["size_fpr_correlation"]) != {str(c) for c, _, _ in CELLS}: raise ValueError("correlation cells")
    for i, v in o["size_fpr_correlation"].items():
        pr = v["per_realisation"]
        if len(pr) != 50 or abs(v["mean"] - sum(pr) / 50) > TOL: raise ValueError(f"correlation mean cell {i}")


def get(o, path):
    for p in path: o = o[p]
    return o


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tree", required=True); ap.add_argument("--src", required=True); a = ap.parse_args()
    tree = Path(a.tree); rm = tree / "results_master.json"; src = Path(a.src)
    try:
        if src.name != RUN_DIR: raise ValueError(f"--src must be the run dir {RUN_DIR}")
        raw = (src / FNAME).read_bytes()
        if sha_bytes(raw) != SUMMARY_SHA: raise ValueError("summary sha differs from SUMMARY_SHA (checked before parsing)")
        rows = [l.split("  ", 1) for l in (src / "SHA256SUMS").read_text().split("\n") if l.strip()]
        if any(len(r) != 2 or sha(src / r[1]) != r[0] for r in rows) or FNAME not in {r[1] for r in rows}: raise ValueError(f"SHA256SUMS check failed in {src}")
        if sha(rm) != MASTER_SHA: raise ValueError("tree results_master differs from MASTER_SHA")
        o = json.loads(raw.decode("utf-8"))
        check_object(o)
        d = json.load(open(rm)); vals, srcs, fsha = d["values"], d["sources"], d["files_sha256"]
        if any(k.startswith(PRE + "/") or k.startswith(PRED + "/") for k in vals): raise ValueError(f"existing keys already use prefix {PRE}/ or {PRED}/")
        new = file_keys(o); der = derived_keys(o, vals)
        clash = [k for k in list(new) + list(der) if k in vals] + [k for k in der if k in new]
        if clash: raise ValueError(f"keys exist already: {clash[:5]}")
    except (ValueError, KeyError, OSError, TypeError, IndexError) as e: print("GATE (exit 3):", e, file=sys.stderr); return 3
    rel_file = f"{DST}/{RUN_DIR}/{FNAME}"
    bak = tree / f"results_master_before_sizefpr_{MASTER_SHA[:8]}.json"; tmp = tree / "results_master.json.sizefpr_tmp"; dst = tree / rel_file
    existing = [str(x) for x in (bak, tmp, dst) if x.exists()] + ([rel_file] if rel_file in fsha else [])
    if existing: print("GATE (exit 3): write targets already exist, nothing written:", existing, file=sys.stderr); return 3
    old_vals = dict(vals); old_srcs = {k: dict(v) if isinstance(v, dict) else v for k, v in srcs.items()}; old_fsha = dict(fsha)
    with open(rm, "rb") as s_, open(bak, "xb") as d_: d_.write(s_.read())
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "xb") as d_: d_.write(raw)
    if sha(dst) != SUMMARY_SHA: print("GATE (exit 3): copy sha differs", file=sys.stderr); return 3
    fsha[rel_file] = SUMMARY_SHA
    for k, (v, path) in new.items(): vals[k] = v; srcs[k] = dict(file=rel_file, path=path)
    for k, (v, s) in der.items(): vals[k] = v; srcs[k] = s
    d["n_values"] = len(vals)
    with open(tmp, "x") as fh: json.dump(d, fh, indent=1, allow_nan=False)
    chk = json.load(open(tmp)); o2 = json.load(open(dst))
    for k, (v, path) in new.items():
        if get(o2, path) != chk["values"].get(k): print("GATE (exit 3): re-read differs", k, file=sys.stderr); return 3
    try: again = derived_keys(o2, {k: v for k, v in chk["values"].items() if k not in new and k not in der})
    except (ValueError, KeyError) as e: print("GATE (exit 3): recompute failed", e, file=sys.stderr); return 3
    for k, (v, _) in der.items():
        if chk["values"].get(k) != v or again[k][0] != v: print("GATE (exit 3): derived re-read differs", k, file=sys.stderr); return 3
    for k, v in old_vals.items():
        if k not in chk["values"] or chk["values"][k] != v or chk["sources"].get(k) != old_srcs.get(k): print("GATE (exit 3): an old key changed", k, file=sys.stderr); return 3
    for k, v in old_fsha.items():
        if chk["files_sha256"].get(k) != v: print("GATE (exit 3): an old files_sha256 entry changed", k, file=sys.stderr); return 3
    if len(chk["values"]) != len(old_vals) + len(new) + len(der): print("GATE (exit 3): key count mismatch", file=sys.stderr); return 3
    tmp.replace(rm)
    print(json.dumps(dict(added=len(new) + len(der), from_file=len(new), derived=len(der), n_values=d["n_values"], new_master_sha256=sha(rm)), indent=1)); return 0


if __name__ == "__main__": sys.exit(main())

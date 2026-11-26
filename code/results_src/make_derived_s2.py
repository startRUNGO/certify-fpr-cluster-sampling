"""Derived numbers for segment 2 (per-field workload and the zero-flag floor). Reads results_master.json (ids dataset/<site>|cells, |fields_all; rev1_results / r8_rev1_results n_labeled_neg). No new analysis.
Run from Paper5_StainRef/: python3 results_src/make_derived_s2.py ; then rebuild results_master.json (derived_s2.json is registered there)."""
import json
from pathlib import Path
from scipy.stats import beta
ROOT = Path(__file__).resolve().parents[1]; V = json.load(open(ROOT/"paper_en/results_master.json"))["values"]
SITES = ["bbbc041_all", "mpidb_full", "lacuna_B", "plasmo", "iml"]; out = {}
for s in SITES:
    cpf = V[f"dataset/{s}|cells"]/V[f"dataset/{s}|fields_all"]
    src = f"r8_rev1_results/agg|a2|iml|22|S3|0.1|n_labeled_neg" if s == "iml" else f"rev1_results/agg|a2|{s}|22|S3|0.1|n_labeled_neg"
    N = V[src]
    out[f"{s}|cells_per_field"] = cpf
    out[f"{s}|false_flags_per_field_eps005"] = 0.05*cpf
    out[f"{s}|zero_flag_floor"] = 1-(0.10/3)**(1.0/N)   # Clopper-Pearson upper limit at level delta/3 with zero exceedances in N cells, no clustering: 1-(delta/3)^(1/N)
    out[f"{s}|zero_flag_floor_N"] = N
json.dump(out, open(ROOT/"results_src/derived_s2.json", "w"), indent=1)
for k, v in out.items(): print(k, v)
# check: closed form equals the Beta quantile with x_eff=0
for s in SITES: assert abs(out[f"{s}|zero_flag_floor"]-float(beta.ppf(1-0.10/3, 1, out[f"{s}|zero_flag_floor_N"]))) < 1e-12

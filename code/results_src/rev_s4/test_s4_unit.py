"""Unit test of the S4 upper bounds on synthetic beta-binomial data only (no site scores). Usage: python3 test_s4_unit.py s4.py TMP_OUT_DIR (the directory must not exist; s4.py creates it on import)."""
import importlib.util, sys, numpy as np
from scipy.stats import beta
s4py, tmp = sys.argv[1], sys.argv[2]
sys.argv = ["s4.py", "/nonexistent_inf", "/nonexistent_e7", tmp, "/nonexistent_rev1.json"]
spec = importlib.util.spec_from_file_location("s4", s4py); S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
rng = np.random.default_rng(7); m = 22; d = 0.10/3; ok = True
for rho, p in ((0.02, 0.06), (0.10, 0.06), (0.25, 0.06), (0.10, 0.15)):
    a = p*(1-rho)/rho; b = (1-p)*(1-rho)/rho
    n = rng.integers(60, 160, m); pf = rng.beta(a, b, m); x = rng.binomial(n, pf)
    u3 = S.b1u(x, n, d); ua = S.ub_s4("S4a", x, n, d); ub = S.ub_s4("S4b", x, n, d); ph = x.sum()/n.sum()
    good = all(np.isfinite(v) and 0 <= v <= 1 for v in (u3, ua, ub)) and ua <= u3 + 1e-12
    print(f"rho={rho} p={p}: p_hat={ph:.4f}  UB S3={u3:.4f}  S4a={ua:.4f}  S4b={ub:.4f}  (S4a<=S3: {ua <= u3 + 1e-12})"); ok &= good
n = np.full(m, 100); x0 = np.zeros(m); N = n.sum()
ub0 = S.ub_s4("S4b", x0, n, d); ref0 = float(beta.ppf(1-d, 1, N)); print("S4b p_hat=0:  UB", ub0, " CP(n_eff=N,x_eff=0)", ref0, " DEFF=1 ->", abs(ub0-ref0) < 1e-12); ok &= abs(ub0-ref0) < 1e-12
ua0 = S.ub_s4("S4a", x0, n, d); print("S4a p_hat=0:  UB", ua0, abs(ua0-ref0) < 1e-12); ok &= abs(ua0-ref0) < 1e-12
n1 = np.full(m, 100); x1 = n1.copy(); u1 = S.ub_s4("S4b", x1, n1, d); u1s3 = S.b1u(x1, n1, d)
print("p_hat=1: S4b UB", u1, " S3 UB", u1s3, "-> both nan (Beta shape n_eff-x_eff=0), the comparison UB<=eps is False, i.e. not certified; same as S3:", bool(np.isnan(u1) and np.isnan(u1s3)) and not (u1 <= 0.5)); ok &= bool(np.isnan(u1) and np.isnan(u1s3)) and not (u1 <= 0.5)
# DEFF recomputed by hand for a clustered draw: n_eff must be < N
rho, p = 0.25, 0.10; a = p*(1-rho)/rho; b = (1-p)*(1-rho)/rho; n = np.full(m, 100); x = rng.binomial(n, rng.beta(a, b, m)); ph = x.sum()/n.sum(); N = n.sum()
v = m/(m-1)*((x-ph*n)**2).sum()/N**2; deff = max(1.0, v/(ph*(1-ph)/N)); print(f"clustered draw: DEFF(S4b)={deff:.2f}, n_eff={N/deff:.0f} of N={N}"); ok &= deff >= 1
print("ALL OK" if ok else "FAILED"); sys.exit(0 if ok else 1)

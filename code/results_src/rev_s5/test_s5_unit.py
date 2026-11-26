"""Unit test of the S5 upper bounds on synthetic beta-binomial data only (no site scores). Usage: python3 test_s5_unit.py s5.py S4_PY TMP_OUT_DIR  (TMP_OUT_DIR must not exist; s5.py creates it on import)."""
import importlib.util, sys, numpy as np
from scipy.stats import beta, norm, t
s5py, s4py, tmp = sys.argv[1:4]
def load(path, name, outdir):
    sys.argv = [path, "/nonexistent_inf", "/nonexistent_e7", outdir, "/nonexistent_rev1.json"]
    sp = importlib.util.spec_from_file_location(name, path); M = importlib.util.module_from_spec(sp); sp.loader.exec_module(M); return M
S5 = load(s5py, "s5", tmp); S4 = load(s4py, "s4", tmp+"_s4")
rng = np.random.default_rng(11); m = 22; d = 0.10/3; ok = True
fac = (norm.ppf(1-d)/t.ppf(1-d, m-1))**2; print(f"(z/t)^2 at delta/3={d:.4f}, df={m-1}: {fac:.4f}")
def neff_s4b(x, n):
    N = n.sum(); ph = x.sum()/N
    if ph <= 0 or ph >= 1: return float(N)
    v = len(n)/(len(n)-1)*((x-ph*n)**2).sum()/N**2; return N/max(1.0, v/(ph*(1-ph)/N))
for rho, p in ((0.02, 0.06), (0.10, 0.06), (0.25, 0.06), (0.10, 0.15)):
    a = p*(1-rho)/rho; b = (1-p)*(1-rho)/rho; n = rng.integers(60, 160, m).astype(float); x = rng.binomial(n.astype(int), rng.beta(a, b, m)).astype(float)
    u3 = S5.b1u(x, n, d); u4b = S4.ub_s4("S4b", x, n, d); u4a = S4.ub_s4("S4a", x, n, d); u5a = S5.ub_s5("S5a", x, n, d); u5b = S5.ub_s5("S5b", x, n, d); ph = x.sum()/n.sum()
    # S5a n_eff must be exactly S4b n_eff * (z/t)^2: recompute the bound from that n_eff with CP at delta
    ne = neff_s4b(x, n)*fac; ref = float(beta.ppf(1-d, ph*ne+1, ne-ph*ne)); c1 = abs(ref-u5a) < 1e-12; c2 = u4b <= u5a + 1e-12; c3 = u5b <= u3 + 1e-12
    print(f"rho={rho} p={p}: p_hat={ph:.4f}  S3={u3:.4f}  S5b={u5b:.4f}  S4a={u4a:.4f}  S5a={u5a:.4f}  S4b={u4b:.4f} | S5a n_eff exact: {c1}  S4b<=S5a: {c2}  S5b<=S3: {c3}  (S5b>=S4a: {u5b >= u4a - 1e-12})")
    ok &= c1 and c2 and c3
n = np.full(m, 100.0); x0 = np.zeros(m); N = n.sum()
u5a0 = S5.ub_s5("S5a", x0, n, d); ref0 = float(beta.ppf(1-d, 1, N*fac)); print(f"p_hat=0: S5a UB {u5a0:.6f} == CP(n_eff=N*(z/t)^2={N*fac:.1f}, x_eff=0) {ref0:.6f}: {abs(u5a0-ref0) < 1e-12}  (S3 UB {S5.b1u(x0, n, d):.6f}, S5b UB {S5.ub_s5('S5b', x0, n, d):.6f})"); ok &= abs(u5a0-ref0) < 1e-12
u5b0 = S5.ub_s5("S5b", x0, n, d); u30 = S5.b1u(x0, n, d); ok &= abs(u5b0-u30) < 1e-12 and abs(u30-float(beta.ppf(1-d, 1, N))) < 1e-12; print("p_hat=0: S5b == S3 == CP(N, 0):", abs(u5b0-u30) < 1e-12)
n1 = np.full(m, 100.0); x1 = n1.copy(); a1 = S5.ub_s5("S5a", x1, n1, d); b1 = S5.ub_s5("S5b", x1, n1, d); print("p_hat=1: S5a", a1, "S5b", b1, "S3", S5.b1u(x1, n1, d), "(nan -> comparison UB<=eps False -> not certified, as in S3)"); ok &= np.isnan(a1) and np.isnan(b1)
print("ALL OK" if ok else "FAILED"); sys.exit(0 if ok else 1)

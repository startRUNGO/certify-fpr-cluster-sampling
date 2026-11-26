# P5a size-FPR check: NEW generator module (实验管理 2026-11-25).  Pinned by sizefpr_check.py; never edited into e12_stress.
# make_site_sizecorr is the generator of P5b p5b_g0g1_sim.py (98d2e898) with rho_c and size_sigma as arguments: E12 cell B field sizes (log-normal, dispersion size_sigma at mean r)
#   and a field shift correlated with the size latent, mu_f = tau (rho_c z_f + sqrt(1 - rho_c^2) w_f).  With rho_c = 0 the draws and arithmetic equal e12_stress.make_site2 for the
#   same B cell bit for bit (with rho_c = 0 the shift is tau * w_f, w_f drawn at the same stream position as make_site2's rng.normal(0, tau, F)); sizefpr_check anchors this.
import numpy as np


def make_site_sizecorr(E10, rng, tau, r, pi, tail, F, size_sigma, rho_c, sep=None):
    """Scores s, labels y, field ids g of one synthetic site."""
    sep = E10.SEP if sep is None else sep
    z = rng.normal(size=F); nf = np.maximum(5, np.round(r * np.exp(size_sigma * z - size_sigma ** 2 / 2)).astype(int))
    npos = rng.poisson(pi / (1 - pi) * nf)
    w = rng.normal(size=F); mu = tau * (rho_c * z + np.sqrt(1 - rho_c ** 2) * w)
    s, y, g = [], [], []
    for f in range(F):
        s.append(mu[f] + E10.noise(rng, nf[f], tail)); y.append(np.zeros(nf[f], int)); g.append(np.full(nf[f], f))
        if npos[f]: s.append(mu[f] + sep + E10.noise(rng, npos[f], tail)); y.append(np.ones(npos[f], int)); g.append(np.full(npos[f], f))
    return np.concatenate(s), np.concatenate(y), np.concatenate(g)

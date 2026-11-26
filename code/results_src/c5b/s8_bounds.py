"""S8 bounds per PREREG_ITER8_S8_20261022.md (sha256 c8037e9a...). B1 Kish+CP, B2 beta-binomial profile LR, B3 cluster-robust ratio."""
import numpy as np
from scipy.special import gammaln
from scipy.stats import beta as betad, norm, t as tdist

def b1(x, n, delta):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); X = x.sum(); ph = X/N
    n0 = (N-(n**2).sum()/N)/(m-1); pi = x/n
    msb = (n*(pi-ph)**2).sum()/(m-1); msw = (n*pi*(1-pi)).sum()/(N-m); den = msb+(n0-1)*msw
    icc = 0.0 if den <= 0 else min(1.0, max(0.0, (msb-msw)/den))
    nbw = (n**2).sum()/N; neff = N/(1+(nbw-1)*icc); xeff = X*neff/N
    return float(betad.ppf(1-delta, xeff+1, neff-xeff)), icc

def b3(x, n, delta):
    x = np.asarray(x, float); n = np.asarray(n, float); m = len(n); N = n.sum(); ph = x.sum()/N; nb = N/m
    se = np.sqrt(((x-ph*n)**2).sum()/(m*(m-1)*nb**2))
    return float(ph+tdist.ppf(1-delta, m-1)*se)

RHO = 1/(1+np.exp(-np.linspace(np.log(1e-4/(1-1e-4)), np.log((1-1e-4)/1e-4), 48)))
PHI = ((1-RHO)/RHO)[:, None, None]   # (R,1,1)
PG = np.concatenate([[1e-6], np.geomspace(1e-4, 0.5, 100), np.linspace(0.52, 0.99, 25)])
def _prof(p, x, n):
    """profile log-lik over rho for each p in array p; x,n (m,)"""
    p = np.asarray(p, float)[None, :, None]; a = p*PHI; b = (1-p)*PHI
    ll = gammaln(x+a)+gammaln(n-x+b)-gammaln(n+PHI)+gammaln(PHI)-gammaln(a)-gammaln(b)
    return ll.sum(2).max(0)   # (P,)
def b2(x, n, delta):
    x = np.asarray(x, float); n = np.asarray(n, float); cut = norm.ppf(1-delta)**2/2
    l = _prof(PG, x, n); lmax = l.max(); ph = x.sum()/n.sum(); k0 = int(np.argmax(l))
    # first grid index at/above the argmax (and above p-hat) where lmax - l > cut
    idx = np.where((PG >= min(PG[k0], max(ph, 1e-6)) - 1e-12) & (lmax-l > cut))[0]
    idx = idx[PG[idx] >= PG[k0]]
    if len(idx) == 0: return 1.0
    k = int(idx[0]); lo, hi = PG[k-1], PG[k]
    fine = np.linspace(lo, hi, 25); lf = _prof(fine, x, n); ok = np.where(lmax-lf > cut)[0]
    if len(ok) == 0: return float(hi)
    j = int(ok[0]); return float(fine[j-1] if j > 0 else lo)

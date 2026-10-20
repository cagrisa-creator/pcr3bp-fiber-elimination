#!/usr/bin/env python3
"""
Global Levi--Civita convexity reconnaissance for the PCR3BP, using the reduced
(fiber-eliminated) gates of the accompanying paper.

For a grid of mass ratios and subcritical energy offsets, this module evaluates:
  * calT_min       -- minimum of the reduced tangential surrogate calT
                      (Lemma 1: exactly min_phi T);
  * calD_bound_min -- minimum of the rigorous determinant lower bound calD
                      (Corollary to Lemma 4: A0 - sqrt(A1^2+B1^2) - 128 F^2 s);
  * D_true_min     -- minimum over the fiber of the exact degree-two trigonometric
                      polynomial for D (Lemma 3), evaluated on the worst base points UNION the
                      entire collar band (v1.3.0 fix; see D_true evaluation below).

The first critical energy c_{L1}(mu) is located by robust bisection of the
effective-potential derivative between the primaries. As a sanity check,
c_{L1}(equal masses) = 2.000000 is recovered.

This is a floating-point reconnaissance, not an interval-certified computation.

Convention (v1.3.0 correction): in the defining function F of Eq. (2.1), mu is
the mass fraction of the DISTANT (non-regularized) primary, located at q=(1,0);
the REGULARIZED primary, at the origin, has mass fraction m_r = 1 - mu. This
follows from the exact identity F = |q|(Omega - c) with Omega the effective
potential of the primary pair (1-mu at the origin, mu at (1,0), barycentre at
(mu,0)); see the symbolic audit in audit/. Earlier versions (<= 1.2.0)
mislabeled m_r as mu itself. c is one half the Jacobi constant; s = |z|^2,
a = z1^2 - z2^2. Requires sympy, numpy, scipy.
Running as a script regenerates results/phase_map_data.csv (now with an
explicit m_r column).

Author: Cesar Augusto Grisa (UFF, Brazil).
"""
import sympy as sp
import numpy as np
import csv
import warnings
warnings.filterwarnings('ignore')
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Symbolic setup: defining function F, gradient, Hessian, structure matrices,
# reduced tangential surrogate calT, and the harmonic coefficients of D.
# ---------------------------------------------------------------------------
# v1.3.2: expression cache. If _gates_cache.py is present (shipped with the
# package and verified by audit/verify_gates_cache.py), load the six gate
# expressions from it -- import then costs seconds on every tested Python.
# Deleting the cache file falls back to the full symbolic derivation below.
_CACHE = None
try:
    from _gates_cache import CACHE as _CACHE
except Exception:
    _CACHE = None

z1, z2, mu, c = sp.symbols('z1 z2 mu c', real=True)
s_ = z1**2 + z2**2
a_ = z1**2 - z2**2
R_ = sp.sqrt(1 - 4*a_ + 4*s_**2)
F = 1 - mu + 2*mu*s_/R_ - 2*c*s_ + mu**2*s_ + 4*s_**3 - 4*mu*s_*a_
if _CACHE is None:
    # Full symbolic derivation (slow path; exercised by
    # audit/verify_gates_cache.py and whenever _gates_cache.py is absent).
    Fz = sp.Matrix([sp.diff(F, z1), sp.diff(F, z2)])
    Fzz = sp.hessian(F, (z1, z2))
    M1 = sp.Matrix([[-4*z2, -4*z1], [-4*z1, -12*z2]])
    M2 = sp.Matrix([[12*z1, 4*z2], [4*z2, 4*z1]])
    L0 = -sp.Rational(1, 2)*Fzz
    adj = lambda M: sp.Matrix([[M[1, 1], -M[0, 1]], [-M[1, 0], M[0, 0]]])
    sig = lambda X, Y: sp.trace(X)*sp.trace(Y) - sp.trace(X*Y)

    # Reduced tangential surrogate (Lemma 1).
    calT = (Fz.T*Fz)[0] - 2*F*sp.trace(Fzz) - 64*F*sp.sqrt(F)*sp.sqrt(s_)

    # Harmonic coefficients of D (Lemma 3, specialized to the LC structure).
    A0e = (Fz.T*adj(L0)*Fz)[0] + 4*F*L0.det() + 64*F**2*s_
    A1e = sp.sqrt(F)*((Fz.T*adj(M1)*Fz)[0] + 4*F*sig(L0, M1))
    B1e = sp.sqrt(F)*((Fz.T*adj(M2)*Fz)[0] + 4*F*sig(L0, M2))
    A2e = -128*F**2*a_
    B2e = -256*F**2*z1*z2
else:
    calT = sp.sympify(_CACHE['calT'])
    A0e = sp.sympify(_CACHE['A0e']); A1e = sp.sympify(_CACHE['A1e'])
    B1e = sp.sympify(_CACHE['B1e']); A2e = sp.sympify(_CACHE['A2e'])
    B2e = sp.sympify(_CACHE['B2e'])

args = (z1, z2, mu, c)
_L = sp.lambdify
fT = _L(args, calT, 'numpy')
fA0 = _L(args, A0e, 'numpy')
fA1 = _L(args, A1e, 'numpy')
fB1 = _L(args, B1e, 'numpy')
fA2 = _L(args, A2e, 'numpy')
fB2 = _L(args, B2e, 'numpy')

# F as a function of (S, Q=cos 2theta) for radial root-finding on each ray.
S, Q = sp.symbols('S Q', real=True)
fF = _L((S, Q, mu, c),
        1 - mu + 2*mu*S/sp.sqrt(1 - 4*S*Q + 4*S**2) - 2*c*S + mu**2*S + 4*S**3 - 4*mu*S**2*Q,
        'numpy')

# Angular sample for the exact fiber polynomial of D.
PH = np.linspace(0, 2*np.pi, 256, endpoint=False)
CP, SP_, C2, S2 = np.cos(PH), np.sin(PH), np.cos(2*PH), np.sin(2*PH)


def cL1_of(muf):
    """First critical (L1) energy for mass parameter muf, via bisection of the
    effective-potential derivative between the primaries. Sanity: cL1_of(0.5)=2."""
    eps = float(min(muf, 1 - muf))
    dOm = lambda x: x - (1 - eps)/(x + eps)**2 + eps/((1 - eps) - x)**2
    xL1 = brentq(dOm, -eps + 1e-9, (1 - eps) - 1e-9, xtol=1e-15)
    return float((xL1**2)/2 + (1 - eps)/(xL1 + eps) + eps/((1 - eps) - xL1))


def scan(muf, cc, n_th=96, n_s=40, n_pre=2000):
    """Scan the base at (muf, cc). Returns (calT_min, calD_bound_min, D_true_min)
    or None if the collar root is not bracketed on some ray."""
    th = np.linspace(0, np.pi/2, n_th)
    q = np.cos(2*th)
    hi = np.where(np.abs(q - 1) < 0.02, 0.499, 1.0)
    Sg = np.linspace(1e-9, 1, n_pre)[:, None]*hi[None, :]
    V = fF(Sg, q[None, :], muf, cc)
    ch = np.diff(np.sign(V), axis=0) != 0
    if not ch.any(axis=0).all():
        return None
    i0 = ch.argmax(axis=0)
    ar = np.arange(n_th)
    lo, hi2 = Sg[i0, ar], Sg[i0 + 1, ar]
    for _ in range(45):  # vectorized bisection to the collar root on each ray
        mid = (lo + hi2)/2
        sw = np.sign(fF(mid, q, muf, cc)) == np.sign(fF(lo, q, muf, cc))
        lo = np.where(sw, mid, lo)
        hi2 = np.where(sw, hi2, mid)
    r0 = (lo + hi2)/2
    sv = (np.arange(1, n_s + 1)/n_s)[:, None]*r0[None, :]
    r = np.sqrt(sv)
    X = r*np.cos(th)[None, :]
    Y = r*np.sin(th)[None, :]
    Tm = float(np.nanmin(fT(X, Y, muf, cc)))
    A0, A1, B1, A2, B2 = (f(X, Y, muf, cc) for f in (fA0, fA1, fB1, fA2, fB2))
    Db = A0 - np.sqrt(A1**2 + B1**2) - np.sqrt(A2**2 + B2**2)
    Dbm = float(np.nanmin(Db))
    # True fiber minimum of D, evaluated exactly on the union of
    #   (i) the 80 worst base points of the lower bound Db, and
    #   (ii) the ENTIRE collar band (the three outermost radial rows on every
    #        ray, s/s_boundary >= (n_s-2)/n_s), where the bound is sharp.
    # Versions <= 1.2.0 used (i) alone; because Db is loose in the interior,
    # its worst points can all be interior while the true minimum of D sits on
    # the collar, producing false 'convex' entries. Adding (ii) removes this
    # failure mode at negligible cost.
    flat = np.argsort(np.nan_to_num(Db, nan=np.inf), axis=None)[:80]
    ii, jj = np.unravel_index(flat, Db.shape)
    n_s_eff, n_th_eff = Db.shape
    rows = np.arange(max(0, n_s_eff - 3), n_s_eff)
    ii_c = np.repeat(rows, n_th_eff)
    jj_c = np.tile(np.arange(n_th_eff), rows.size)
    ii = np.concatenate([ii, ii_c]); jj = np.concatenate([jj, jj_c])
    Dtrue = (A0[ii, jj][:, None] + A1[ii, jj][:, None]*CP + B1[ii, jj][:, None]*SP_
             + A2[ii, jj][:, None]*C2 + B2[ii, jj][:, None]*S2)
    Dtm = float(np.nanmin(Dtrue))
    return Tm, Dbm, Dtm


if __name__ == '__main__':
    print("sanity: c_L1(0.5)=%.6f (expected 2.0), c_L1(Earth-Moon)=%.9f"
          % (cL1_of(0.5), cL1_of(0.987849414390376)))
    mus = [0.001, 0.005, 0.012150585609624, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80,
           0.90, 0.941, 0.97, 0.987849414390376, 0.999]
    dcs = [1e-4, 1e-3, 0.01, 0.03, 0.1, 0.3, 1.0]
    rows = []
    for muf in mus:
        cl = cL1_of(muf)
        for d in dcs:
            out = scan(muf, cl + d, n_th=192, n_s=80)
            rows.append((muf, cl, d) + (out if out else (None, None, None)))
        print("done mu=%.4f" % muf, flush=True)
    with open('../results/phase_map_data.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['mu_formula', 'm_r', 'c_L1', 'delta_c',
                    'calT_min', 'calD_bound_min', 'D_true_min'])
        w.writerows([(r[0], 1 - r[0]) + r[1:] for r in rows])
    print("wrote ../results/phase_map_data.csv")

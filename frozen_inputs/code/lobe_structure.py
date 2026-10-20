#!/usr/bin/env python3
"""
Fine structure of the non-convexity lobes (Observation 6.2 of the paper).

Two quantitative facts about the off-axis lobes where the determinant gate is
negative, across the comparable-mass band (c = c_L1 + 1e-3):

  (F1) Simple-fold structure. At every base point inside a lobe (true fiber
       minimum of D negative), the reduced tangential gate calT stays strictly
       positive, with a uniform margin. The tangential gate controls the trace
       of the convexity pencil and the determinant gate its determinant, so a
       point with D<0<calT carries exactly ONE negative principal direction:
       the non-convexity is a single simple fold, never a double concavity.

  (F2) Radial transversality (fiberwise star-shapedness of the base).
       Along every ray the defining function F crosses the collar with
       dF/dr < 0, uniformly, including on rays through the lobes. The base
       component is radially star-shaped and the LC non-convexity does not
       degrade the star-shape margin.

Floating-point observation at a finite grid, not a proof. Conventions match
convexity_map.py / dynamical_relevance.py: mu = DISTANT-primary mass fraction
(v1.3.0 correction; the regularized primary has m_r = 1-mu),
c = one half the Jacobi constant, s = |z|^2, a = z1^2 - z2^2.

To keep evaluation fast we assemble the harmonic coefficients A0..B2 and the
reduced tangential gate calT NUMERICALLY from cheap lambdified blocks
(F, grad F, Hess F), avoiding the costly symbolic build of the bordered
determinant. The formulas are exactly those of Lemmas 3.1-3.5 of the paper.

Author: Cesar Augusto Grisa (UFF, Brazil).
"""
import numpy as np
import sympy as sp
from scipy.optimize import brentq


def cL1_of(muf):
    """First critical (L1) energy for mass parameter muf, via bisection of the
    effective-potential derivative between the primaries. Sanity: cL1_of(0.5)=2.
    Identical to convexity_map.cL1_of; inlined here so importing this module
    does not trigger the (slow) symbolic build in convexity_map's top level."""
    eps = float(min(muf, 1 - muf))
    dOm = lambda x: x - (1 - eps)/(x + eps)**2 + eps/((1 - eps) - x)**2
    xL1 = brentq(dOm, -eps + 1e-9, (1 - eps) - 1e-9, xtol=1e-15)
    return float((xL1**2)/2 + (1 - eps)/(xL1 + eps) + eps/((1 - eps) - xL1))

# ------------------------------------------------ cheap symbolic blocks -----
z1, z2, mu, c = sp.symbols('z1 z2 mu c', real=True)
s_ = z1**2 + z2**2
a_ = z1**2 - z2**2
R_ = sp.sqrt(1 - 4*a_ + 4*s_**2)
F = 1 - mu + 2*mu*s_/R_ - 2*c*s_ + mu**2*s_ + 4*s_**3 - 4*mu*s_*a_
args = (z1, z2, mu, c)
_L = sp.lambdify
fF = _L(args, F, 'numpy')
fFx = _L(args, sp.diff(F, z1), 'numpy')
fFy = _L(args, sp.diff(F, z2), 'numpy')
fFxx = _L(args, sp.diff(F, z1, 2), 'numpy')
fFxy = _L(args, sp.diff(F, z1, z2), 'numpy')
fFyy = _L(args, sp.diff(F, z2, 2), 'numpy')

# F on a ray, F(S, Q=cos 2theta), for collar root-finding and dF/dr.
S, Q = sp.symbols('S Q', real=True)
F_ray = (1 - mu + 2*mu*S/sp.sqrt(1 - 4*S*Q + 4*S**2)
         - 2*c*S + mu**2*S + 4*S**3 - 4*mu*S**2*Q)
fFray = _L((S, Q, mu, c), F_ray, 'numpy')
fdFds = _L((S, Q, mu, c), sp.diff(F_ray, S), 'numpy')

PH = np.linspace(0, 2*np.pi, 192, endpoint=False)
CP, SP_, C2, S2 = np.cos(PH), np.sin(PH), np.cos(2*PH), np.sin(2*PH)


def _blocks(x, y, mf, cc):
    Fv = fF(x, y, mf, cc)
    gx, gy = fFx(x, y, mf, cc), fFy(x, y, mf, cc)
    hxx, hxy, hyy = fFxx(x, y, mf, cc), fFxy(x, y, mf, cc), fFyy(x, y, mf, cc)
    return Fv, gx, gy, hxx, hxy, hyy


def gate_values(x, y, mf, cc):
    """Return (calT, and the harmonic coeffs A0,A1,B1,A2,B2 of D) at (x,y).
    All assembled numerically from the cheap blocks, per Lemmas 3.1-3.5."""
    Fv, gx, gy, hxx, hxy, hyy = _blocks(x, y, mf, cc)
    s = x*x + y*y
    a = x*x - y*y
    # L0 = -1/2 Hess F
    l11, l12, l22 = -0.5*hxx, -0.5*hxy, -0.5*hyy
    # structure matrices
    m1_11, m1_12, m1_22 = -4*y, -4*x, -12*y
    m2_11, m2_12, m2_22 = 12*x, 4*y, 4*x
    # g^T adj(X) g for symmetric X = [[a11,a12],[a12,a22]]:
    #   = a22 gx^2 - 2 a12 gx gy + a11 gy^2
    gAg = lambda a11, a12, a22: a22*gx*gx - 2*a12*gx*gy + a11*gy*gy
    detL0 = l11*l22 - l12*l12
    trL0 = l11 + l22
    sig = lambda a11, a12, a22: trL0*(a11 + a22) - (l11*a11 + 2*l12*a12 + l22*a22)
    # tangential gate calT (Lemma 3.1): |Fz|^2 - 2 F Lap F - 64 F^{3/2} sqrt(s)
    lapF = hxx + hyy
    calT = (gx*gx + gy*gy) - 2*Fv*lapF - 64*Fv*np.sqrt(max(Fv, 0.0))*np.sqrt(s)
    # harmonic coeffs of D (Lemma 3.4/3.5 specialization)
    A0 = gAg(l11, l12, l22) + 4*Fv*detL0 + 64*Fv*Fv*s
    rt = np.sqrt(max(Fv, 0.0))
    A1 = rt*(gAg(m1_11, m1_12, m1_22) + 4*Fv*sig(m1_11, m1_12, m1_22))
    B1 = rt*(gAg(m2_11, m2_12, m2_22) + 4*Fv*sig(m2_11, m2_12, m2_22))
    A2 = -128*Fv*Fv*a
    B2 = -256*Fv*Fv*x*y
    return calT, A0, A1, B1, A2, B2


def Dmin(x, y, mf, cc):
    _, A0, A1, B1, A2, B2 = gate_values(x, y, mf, cc)
    return np.min(A0 + A1*CP + B1*SP_ + A2*C2 + B2*S2)


def calT_at(x, y, mf, cc):
    return gate_values(x, y, mf, cc)[0]


def sstar(th, mf, cc, n=4000):
    q = np.cos(2*th)
    hi = 0.499 if abs(q - 1) < 0.02 else 1.0
    ss = np.linspace(1e-9, hi, n)
    v = fFray(ss, q, mf, cc)
    idx = np.where(np.diff(np.sign(v)) != 0)[0]
    if not len(idx):
        return None
    return brentq(lambda x: fFray(x, q, mf, cc), ss[idx[0]], ss[idx[0] + 1], xtol=1e-15)


def analyze(mf, dc=1e-3, n_theta=200, n_rad=80):
    cc = cL1_of(mf) + dc
    lobe_angles = []
    minT_lobe = np.inf
    minD_lobe = np.inf
    max_dFdr_collar = -np.inf
    for th in np.linspace(1e-4, np.pi/2 - 1e-4, n_theta):
        s0 = sstar(th, mf, cc)
        if s0 is None:
            continue
        q = np.cos(2*th)
        dFdr = 2*np.sqrt(s0)*fdFds(s0, q, mf, cc)   # dF/dr = 2 r F_s at r=sqrt(s0)
        max_dFdr_collar = max(max_dFdr_collar, dFdr)
        ray_has_lobe = False
        for fr in np.linspace(0.5, 0.9999, n_rad):
            r = np.sqrt(fr*s0)
            x, y = r*np.cos(th), r*np.sin(th)
            if Dmin(x, y, mf, cc) < 0:
                ray_has_lobe = True
                minD_lobe = min(minD_lobe, Dmin(x, y, mf, cc))
                minT_lobe = min(minT_lobe, calT_at(x, y, mf, cc))
        if ray_has_lobe:
            lobe_angles.append(np.degrees(th))
    band = (min(lobe_angles), max(lobe_angles)) if lobe_angles else None
    return dict(mu=mf, c=cc, band=band,
                minT=minT_lobe if lobe_angles else None,
                minD=minD_lobe if lobe_angles else None,
                max_dFdr=max_dFdr_collar)


def main():
    print("Fine structure of the non-convexity lobes  (c = c_L1 + 1e-3)")
    print("=" * 84)
    print("%-6s %-20s %-15s %-16s %-16s"
          % ("mu", "lobe band (deg)", "min D in lobe", "min calT in lobe", "max dF/dr collar"))
    ok_fold = True
    ok_star = True
    for mf in (0.20, 0.35, 0.50, 0.65, 0.80):
        r = analyze(mf)
        band = ("[%.1f, %.1f]" % r['band']) if r['band'] else "none"
        minT = ("%+.3e" % r['minT']) if r['minT'] is not None else "-"
        minD = ("%+.3e" % r['minD']) if r['minD'] is not None else "-"
        print("%-6.2f %-20s %-15s %-16s %+.3e"
              % (r['mu'], band, minD, minT, r['max_dFdr']))
        if r['minT'] is not None and r['minT'] <= 0:
            ok_fold = False
        if r['max_dFdr'] >= 0:
            ok_star = False
    print("-" * 84)
    print("simple fold (calT>0 at every lobe point) : %s" % ("PASS" if ok_fold else "FAIL"))
    print("radial transversality (dF/dr<0 at collar): %s" % ("PASS" if ok_star else "FAIL"))
    if ok_fold and ok_star:
        print("PASS_LOBE_STRUCTURE")


if __name__ == '__main__':
    main()

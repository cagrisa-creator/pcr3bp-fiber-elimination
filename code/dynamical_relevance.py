#!/usr/bin/env python3
"""
Dynamical-relevance observation for the PCR3BP non-convexity window.

In the comparable-mass regime the Levi--Civita determinant gate D becomes
negative in a pair of lateral lobes hugging the collar. This script localizes
those lobes in base polar coordinates and measures the true fiber minimum of D
along the two coordinate axes, which carry the retrograde and direct collision
orbits of Birkhoff. The finding (see the accompanying paper, Observation 6.1):
the lobes sit off-axis (angular offset >= a few degrees), while D is strictly
positive along both axes.

Outputs paper/figures/dynamical_relevance.png and prints the on-axis minima.
This is a floating-point observation at a finite grid, not a proof.

Convention: mu is the large-primary mass parameter; c is one half the Jacobi
constant; s = |z|^2. Requires sympy, numpy, scipy, matplotlib.

Author: Cesar Augusto Grisa (UFF, Brazil).
"""
import sys
import numpy as np
import sympy as sp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '.')
from convexity_map import cL1_of  # noqa: E402

# Harmonic coefficients of D (Lemma 3), built locally to keep this script standalone.
z1, z2, mu, c = sp.symbols('z1 z2 mu c', real=True)
s_ = z1**2 + z2**2
a_ = z1**2 - z2**2
R_ = sp.sqrt(1 - 4*a_ + 4*s_**2)
F = 1 - mu + 2*mu*s_/R_ - 2*c*s_ + mu**2*s_ + 4*s_**3 - 4*mu*s_*a_
Fz = sp.Matrix([sp.diff(F, z1), sp.diff(F, z2)])
Fzz = sp.hessian(F, (z1, z2))
M1 = sp.Matrix([[-4*z2, -4*z1], [-4*z1, -12*z2]])
M2 = sp.Matrix([[12*z1, 4*z2], [4*z2, 4*z1]])
L0 = -sp.Rational(1, 2)*Fzz
adj = lambda M: sp.Matrix([[M[1, 1], -M[0, 1]], [-M[1, 0], M[0, 0]]])
sig = lambda X, Y: sp.trace(X)*sp.trace(Y) - sp.trace(X*Y)
A0e = (Fz.T*adj(L0)*Fz)[0] + 4*F*L0.det() + 64*F**2*s_
A1e = sp.sqrt(F)*((Fz.T*adj(M1)*Fz)[0] + 4*F*sig(L0, M1))
B1e = sp.sqrt(F)*((Fz.T*adj(M2)*Fz)[0] + 4*F*sig(L0, M2))
A2e = -128*F**2*a_
B2e = -256*F**2*z1*z2
args = (z1, z2, mu, c)
_L = sp.lambdify
fA0, fA1, fB1, fA2, fB2 = (_L(args, e, 'numpy') for e in (A0e, A1e, B1e, A2e, B2e))
Sy, Qy = sp.symbols('S Q')
fF = _L((Sy, Qy, mu, c),
        1 - mu + 2*mu*Sy/sp.sqrt(1 - 4*Sy*Qy + 4*Sy**2) - 2*c*Sy + mu**2*Sy + 4*Sy**3 - 4*mu*Sy**2*Qy,
        'numpy')

PH = np.linspace(0, 2*np.pi, 192, endpoint=False)
CP, SP_, C2, S2 = np.cos(PH), np.sin(PH), np.cos(2*PH), np.sin(2*PH)


def Dmin(x, y, mf, cc):
    """True fiber minimum of D at base point (x, y) via the exact trig polynomial."""
    A0, A1, B1, A2, B2 = (f(x, y, mf, cc) for f in (fA0, fA1, fB1, fA2, fB2))
    return np.min(A0 + A1*CP + B1*SP_ + A2*C2 + B2*S2)


def sstar(th, mf, cc, n=4000):
    """Collar radius s* along the ray at angle th."""
    q = np.cos(2*th)
    hi = 0.499 if abs(q - 1) < 0.02 else 1.0
    ss = np.linspace(1e-9, hi, n)
    v = fF(ss, q, mf, cc)
    idx = np.where(np.diff(np.sign(v)) != 0)[0]
    return brentq(lambda x: fF(x, q, mf, cc), ss[idx[0]], ss[idx[0] + 1], xtol=1e-15) if len(idx) else None


def main():
    mf, cc = 0.5, cL1_of(0.5) + 1e-3
    print("Equal masses, c = c_L1 + 1e-3 = %.6f" % cc)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5.4))

    # Panel 1: where D<0 lives (the lobes), with the two orbit axes.
    for t in np.linspace(0, 2*np.pi, 720, endpoint=False):
        tq = abs((t) % (np.pi/2))
        r0 = sstar(min(tq, np.pi/2 - 1e-6) if tq < np.pi/2 else np.pi/4, mf, cc)
        if r0 is None:
            continue
        for fr in np.linspace(0.5, 1.0, 60):
            r = np.sqrt(fr*r0)
            x, y = r*np.cos(t), r*np.sin(t)
            if Dmin(x, y, mf, cc) < 0:
                ax[0].plot(x, y, 's', color='#d62728', ms=2.5)
    tt = np.linspace(0, 2*np.pi, 400)
    rb = [np.sqrt(sstar(min(abs((t) % (np.pi/2)), np.pi/2 - 1e-6), mf, cc) or 0) for t in tt]
    ax[0].plot(np.array(rb)*np.cos(tt), np.array(rb)*np.sin(tt), '-', color='#888', lw=1, label='collar (F=0)')
    ax[0].plot([-0.75, 0.75], [0, 0], '-', color='#1f77b4', lw=2, label='retrograde-orbit axis')
    ax[0].plot([0, 0], [-0.75, 0.75], '-', color='#2ca02c', lw=1.5, label='direct-orbit axis')
    ax[0].plot([], [], 's', color='#d62728', ms=6, label='D<0 (non-convex)')
    ax[0].set_aspect('equal')
    ax[0].set_title('Where D<0 lives (equal masses, $c=c_{L_1}+10^{-3}$)')
    ax[0].legend(fontsize=8, loc='upper right')
    ax[0].set_xlabel('$z_1$'); ax[0].set_ylabel('$z_2$')

    # Panel 2: D along on-axis vs off-axis rays.
    frs = np.linspace(0.1, 0.999, 60)
    on_axis_min = None
    for ang, lbl, col in [(0.0, r'on-axis ($\theta=0°$)', '#1f77b4'),
                          (np.radians(7), r'lobe ($\theta=7°$)', '#d62728'),
                          (np.radians(45), r'diagonal ($\theta=45°$)', '#2ca02c')]:
        r0 = sstar(min(ang, np.pi/2 - 1e-6), mf, cc)
        D = [Dmin(np.sqrt(fr*r0)*np.cos(ang), np.sqrt(fr*r0)*np.sin(ang), mf, cc) for fr in frs]
        ax[1].plot(frs, D, '-', color=col, label=lbl)
        if ang == 0.0:
            on_axis_min = min(D)
    ax[1].axhline(0, color='k', lw=0.8, ls='--')
    ax[1].set_xlabel('radial fraction $s/s^*$ (1.0 = collar)')
    ax[1].set_ylabel('true determinant minimum over fiber')
    ax[1].set_title('D stays positive on-axis, dips below zero only off-axis')
    ax[1].legend(fontsize=9)
    ax[1].grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('../paper/figures/dynamical_relevance.png', dpi=150)

    # On-axis minima for both orbit axes.
    for th, name in [(0.0, 'retrograde axis (z2=0)'), (np.pi/2, 'direct axis (z1=0)')]:
        r0 = sstar(th, mf, cc)
        vals = [Dmin(np.sqrt(fr*r0)*np.cos(th), np.sqrt(fr*r0)*np.sin(th), mf, cc)
                for fr in np.linspace(0.1, 0.999, 50)]
        print("  min D along %s = %+.4e  (positive: %s)" % (name, min(vals), min(vals) > 0))
    print("wrote ../paper/figures/dynamical_relevance.png")


if __name__ == '__main__':
    main()

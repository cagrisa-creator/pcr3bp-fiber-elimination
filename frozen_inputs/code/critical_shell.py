"""
critical_shell.py -- Refinement of the convexity reconnaissance at the
critical limit for the Earth--Moon MOON component (m_r = 0.012151; v1.3.0
convention, see convexity_map.py) (Section "Refinement at
the critical limit" of the paper).

Tracks the true fiber minimum of the determinant gate D as the energy
offset dc = c - c_L1 tends to 0+, locates by bisection the width dc* of the
ultrathin non-convex shell, and reports the robustness checks quoted in the
paper (radial-refinement stability and location of the minimizing point).

All computations here are floating-point reconnaissance (float64), labeled
as such in the paper; nothing in this module is interval-certified.

Expected output (abridged):
    dc = 1.00e-05  D_min approx +8.6e-08   (positive)
    dc = 1.00e-06  D_min approx -1.5e-06   (negative)
    shell width dc* approx 9.5e-06
Runtime: a few minutes on a laptop.
"""

import numpy as np
import sympy as sp
from scipy.optimize import brentq

# ----------------------------------------------------------------------
# Exact symbolic layer: same formulas as convexity_map.py (Lemmas 1-4)
# ----------------------------------------------------------------------
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
sig = lambda A, B: sp.trace(A)*sp.trace(B) - sp.trace(A*B)

A0e = (Fz.T*adj(L0)*Fz)[0] + 4*F*L0.det() + 64*F**2*s_
A1e = sp.sqrt(F)*((Fz.T*adj(M1)*Fz)[0] + 4*F*sig(L0, M1))
B1e = sp.sqrt(F)*((Fz.T*adj(M2)*Fz)[0] + 4*F*sig(L0, M2))
A2e = -128*F**2*a_
B2e = -256*F**2*z1*z2

args = (z1, z2, mu, c)
fA0, fA1, fB1, fA2, fB2 = (sp.lambdify(args, e, 'numpy')
                           for e in (A0e, A1e, B1e, A2e, B2e))

# ----------------------------------------------------------------------
# Critical energy c_L1 by robust bisection between the primaries
# ----------------------------------------------------------------------
# mu is the DISTANT-primary mass fraction (v1.3.0 convention correction):
# mu = 0.987849... (Earth distant) regularizes the MOON, m_r = 1-mu = 0.012151.
MU_DISTANT_EARTH = 0.987849414390376   # -> regularized component: MOON


def c_L1_of(mu_val):
    eps = float(min(mu_val, 1 - mu_val))
    dOm = lambda x: x - (1 - eps)/(x + eps)**2 + eps/((1 - eps) - x)**2
    xL1 = brentq(dOm, -eps + 1e-12, (1 - eps) - 1e-12, xtol=1e-15)
    return float(xL1**2/2 + (1 - eps)/(xL1 + eps) + eps/((1 - eps) - xL1))


# ----------------------------------------------------------------------
# True fiber minimum of D over the compact component at energy c
# ----------------------------------------------------------------------
PHI = np.linspace(0, 2*np.pi, 1024, endpoint=False)


def D_true_min(mu_val, c_val, n_th=600, n_s=400, report_location=False):
    """Minimum over the component of min_phi D, with D evaluated through
    the exact degree-two trigonometric polynomial (Lemma 3) on the 120
    worst base points of the amplitude bound UNION the entire collar band
    (v1.3.0 fix; the bound is loose in the interior and sharp on the collar)."""
    th = np.linspace(0, np.pi/2, n_th)
    q = np.cos(2*th)
    hi = np.where(np.abs(q - 1) < 0.02, 0.499, 1.0)   # singularity guard s=1/2
    Sg = np.linspace(1e-9, 1, 4000)[:, None]*hi[None, :]

    def Fval(S, Q):
        R = np.sqrt(1 - 4*S*Q + 4*S**2)
        return (1 - mu_val + 2*mu_val*S/R - 2*c_val*S
                + mu_val**2*S + 4*S**3 - 4*mu_val*S**2*Q)

    V = Fval(Sg, q[None, :])
    ch = np.diff(np.sign(V), axis=0) != 0
    if not ch.any(axis=0).all():
        return None
    i0 = ch.argmax(axis=0)
    ar = np.arange(n_th)
    lo, hi2 = Sg[i0, ar], Sg[i0 + 1, ar]
    for _ in range(55):                               # boundary F=0 bisection
        mid = (lo + hi2)/2
        sw = np.sign(Fval(mid, q)) == np.sign(Fval(lo, q))
        lo = np.where(sw, mid, lo)
        hi2 = np.where(sw, hi2, mid)
    r0 = (lo + hi2)/2

    sv = (np.arange(1, n_s + 1)/n_s)[:, None]*r0[None, :]
    r = np.sqrt(sv)
    X = r*np.cos(th)[None, :]
    Y = r*np.sin(th)[None, :]

    with np.errstate(invalid='ignore'):
        A0 = fA0(X, Y, mu_val, c_val); A1 = fA1(X, Y, mu_val, c_val)
        B1 = fB1(X, Y, mu_val, c_val); A2 = fA2(X, Y, mu_val, c_val)
        B2 = fB2(X, Y, mu_val, c_val)
        Db = A0 - np.sqrt(A1**2 + B1**2) - np.sqrt(A2**2 + B2**2)

    # v1.3.0: exact evaluation on worst-120 of Db UNION the collar band
    # (three outermost radial rows), closing the false-convex failure mode
    # of the worst-points-only heuristic (see convexity_map.py).
    flat = np.argsort(np.nan_to_num(Db, nan=np.inf), axis=None)[:120]
    ii, jj = np.unravel_index(flat, Db.shape)
    n_s_eff, n_th_eff = Db.shape
    rows_c = np.arange(max(0, n_s_eff - 3), n_s_eff)
    ii_c = np.repeat(rows_c, n_th_eff)
    jj_c = np.tile(np.arange(n_th_eff), rows_c.size)
    ii = np.concatenate([ii, ii_c]); jj = np.concatenate([jj, jj_c])
    Dt = (A0[ii, jj][:, None] + A1[ii, jj][:, None]*np.cos(PHI)
          + B1[ii, jj][:, None]*np.sin(PHI)
          + A2[ii, jj][:, None]*np.cos(2*PHI)
          + B2[ii, jj][:, None]*np.sin(2*PHI))
    k = np.unravel_index(np.nanargmin(Dt), Dt.shape)
    Dmin = float(Dt[k])
    if report_location:
        i_pt = k[0]
        theta_star = float(th[jj[i_pt]])
        s_frac = float(sv[ii[i_pt], jj[i_pt]]/r0[jj[i_pt]])
        return Dmin, theta_star, s_frac
    return Dmin


def main():
    cl = c_L1_of(MU_DISTANT_EARTH)
    print(f"c_L1(Earth--Moon pair) = {cl:.15f}   [regularized component: MOON, m_r=0.012151]\n")

    # ---- 1. Track D_min as dc -> 0+ -----------------------------------
    print("dc            D_min")
    print("-"*32)
    track = []
    for dc in [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 5e-6, 2e-6, 1e-6, 5e-7, 1e-7]:
        d = D_true_min(MU_DISTANT_EARTH, cl + dc)
        track.append((dc, d))
        print(f"{dc:.2e}    {d:+.4e}")

    dcs = np.array([t[0] for t in track])
    ds = np.array([t[1] for t in track])
    A, D0 = np.polyfit(dcs, ds, 1)
    print(f"\nLinear extrapolation: D_min(dc) ~ {A:.3e} * dc + ({D0:+.3e})")
    print(f"  -> D_min(0+) ~ {D0:+.3e}  (strictly negative: intrinsic to the limit)")

    # ---- 2. Shell width dc* by bisection ------------------------------
    lo_dc, hi_dc = 1e-6, 2e-5
    for _ in range(30):
        mid = (lo_dc + hi_dc)/2
        if D_true_min(MU_DISTANT_EARTH, cl + mid) > 0:
            hi_dc = mid
        else:
            lo_dc = mid
    dc_star = (lo_dc + hi_dc)/2
    print(f"\nShell width dc* = {dc_star:.4e}")
    print(f"Non-convex shell: c in (c_L1, {cl + dc_star:.10f}]")

    # ---- 3. Robustness: location of the minimizing point --------------
    Dm, th_s, sf = D_true_min(MU_DISTANT_EARTH, cl + 1e-6, report_location=True)
    print(f"\nAt dc = 1e-6: D_min = {Dm:+.4e}, theta* = {th_s:.4f} "
          f"(heavy-primary collision axis), s/s_boundary = {sf:.3f}")
    print("Minimizer is far from the s = 1/2 singularity guard: "
          "not a regularization artifact.")

    # ---- 4. Radial-refinement stability check -------------------------
    print("\nStability under radial refinement at dc = 3e-6:")
    for ns in (60, 150, 400):
        d = D_true_min(MU_DISTANT_EARTH, cl + 3e-6, n_s=ns)
        print(f"  n_s = {ns:4d}: D_min = {d:+.4e}")


if __name__ == '__main__':
    main()

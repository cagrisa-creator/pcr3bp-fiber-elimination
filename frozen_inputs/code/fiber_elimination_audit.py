#!/usr/bin/env python3
"""
Exact symbolic audit of the fiber-elimination lemmas for the PCR3BP
Levi--Civita convexity gates.

This script verifies, by exact computer algebra, the eight identities that
underlie Lemmas 1--4 and Theorem 1 of the accompanying paper. Each identity is
formed as a residual (difference of the two sides) and simplified to the symbol
0. If every residual is identically zero, the script prints
    PASS_C82_FIBER_ELIMINATION_AUDIT
and exits with status 0; otherwise it exits 1.

Residual map (see the paper for the corresponding equation numbers):

  R1  Gap formula (Lemma 1):
      T(z,phi) = calT(z) + 64 F^{3/2} ( sqrt(s) - z2 cos(phi) + z1 sin(phi) ),
      calT = |Fz|^2 - 2 F LapF - 64 F^{3/2} sqrt(s).
  R2  Perfect square (Lemma 1):
      s - (-z2 cos phi + z1 sin phi)^2 = (z1 cos phi + z2 sin phi)^2.
  R3  Attained angle (Lemma 1).
  R4  Gradient reduction (Lemma 2): |Fz|^2 = 4 s Fs^2 + 8 a Fs Fa + 4 s Fa^2.
  R5  Laplacian reduction (Lemma 2): LapF = 4 s Fss + 8 a Fsa + 4 s Faa + 4 Fs.
  R6  Generic harmonic decomposition (Lemma 3): pure multilinear identity for
      arbitrary symmetric 2x2 L0, M1, M2, vector g, and scalars rho, F.
  R7  Determinant sum (Lemma 4): det M1 + det M2 = 32 s.
  R8  Amplitude (Lemma 4): (det M1 - det M2)^2 + sigma(M1,M2)^2 = (64 s)^2.

Convention: mu is the large-primary mass parameter; c is one half the Jacobi
constant; s = |z|^2, a = z1^2 - z2^2. Runtime: a few seconds. Requires sympy.

Author: Cesar Augusto Grisa (UFF, Brazil).
"""
from __future__ import annotations
import sympy as sp

def _iszero(expr):
    """Exact zero-test with cheap exact pipelines tried before full simplify.
    Every step is symbolically exact; only the ORDER changes (v1.3.1 speedup,
    identical verdicts). Returns the residual in simplest achieved form."""
    import sympy as _sp
    for f in (_sp.expand, lambda e: _sp.cancel(_sp.together(e)),
              _sp.radsimp, _sp.expand_trig):
        try:
            expr = f(expr)
        except Exception:
            continue
        if expr == 0:
            return _sp.S.Zero
    return __import__('sympy').simplify(expr)




def main() -> int:
    ok = True
    z1, z2, mu, c, phi = sp.symbols('z1 z2 mu c phi', real=True)
    s_ = z1**2 + z2**2
    a_ = z1**2 - z2**2
    R_ = sp.sqrt(1 - 4*a_ + 4*s_**2)
    F = 1 - mu + 2*mu*s_/R_ - 2*c*s_ + mu**2*s_ + 4*s_**3 - 4*mu*s_*a_
    Fz = sp.Matrix([sp.diff(F, z1), sp.diff(F, z2)])
    Fzz = sp.hessian(F, (z1, z2))
    M1 = sp.Matrix([[-4*z2, -4*z1], [-4*z1, -12*z2]])
    M2 = sp.Matrix([[12*z1, 4*z2], [4*z2, 4*z1]])
    sqF = sp.sqrt(F)
    u1, u2 = sqF*sp.cos(phi), sqF*sp.sin(phi)
    L = -sp.Rational(1, 2)*Fzz + u1*M1 + u2*M2
    T = (Fz.T*Fz)[0] + 4*F*sp.trace(L)

    # --- Lemma 1 ---
    calT = (Fz.T*Fz)[0] - 2*F*sp.trace(Fzz) - 64*F*sqF*sp.sqrt(s_)
    r = _iszero(T - calT - 64*F*sqF*(sp.sqrt(s_) - z2*sp.cos(phi) + z1*sp.sin(phi)))
    print('R1_GAP_RESIDUAL            =', r); ok &= (r == 0)
    r = _iszero(s_ - (-z2*sp.cos(phi) + z1*sp.sin(phi))**2 - (z1*sp.cos(phi) + z2*sp.sin(phi))**2)
    print('R2_PERFECT_SQUARE_RESIDUAL =', r); ok &= (r == 0)
    r = _iszero(sp.sqrt(s_) - z2*(z2/sp.sqrt(s_)) + z1*(-z1/sp.sqrt(s_)))
    print('R3_ATTAINED_RESIDUAL       =', r); ok &= (r == 0)

    # --- Lemma 2: reduction to (s, a) ---
    S, A = sp.symbols('S A', real=True)
    Fsa = 1 - mu + 2*mu*S/sp.sqrt(1 - 4*A + 4*S**2) - 2*c*S + mu**2*S + 4*S**3 - 4*mu*S*A
    FS, FA = sp.diff(Fsa, S), sp.diff(Fsa, A)
    FSS, FSA, FAA = sp.diff(Fsa, S, 2), sp.diff(Fsa, S, A), sp.diff(Fsa, A, 2)
    sub = {S: s_, A: a_}
    r = _iszero((Fz.T*Fz)[0] - (4*S*FS**2 + 8*A*FS*FA + 4*S*FA**2).subs(sub))
    print('R4_GRAD_SA_RESIDUAL        =', r); ok &= (r == 0)
    r = _iszero(sp.trace(Fzz) - (4*S*FSS + 8*A*FSA + 4*S*FAA + 4*FS).subs(sub))
    print('R5_LAP_SA_RESIDUAL         =', r); ok &= (r == 0)

    # --- Lemma 3: generic harmonic decomposition (multilinear identity) ---
    g1, g2, rho = sp.symbols('g1 g2 rho', real=True)
    l11, l12, l22 = sp.symbols('l11 l12 l22', real=True)
    m11, m12, m22 = sp.symbols('m11 m12 m22', real=True)
    n11, n12, n22 = sp.symbols('n11 n12 n22', real=True)
    Fv = sp.symbols('Fv', real=True)
    L0g = sp.Matrix([[l11, l12], [l12, l22]])
    M1g = sp.Matrix([[m11, m12], [m12, m22]])
    M2g = sp.Matrix([[n11, n12], [n12, n22]])
    g = sp.Matrix([g1, g2])
    adj = lambda M: sp.Matrix([[M[1, 1], -M[0, 1]], [-M[1, 0], M[0, 0]]])
    sig = lambda X, Y: sp.trace(X)*sp.trace(Y) - sp.trace(X*Y)
    v1, v2 = rho*sp.cos(phi), rho*sp.sin(phi)
    Lg = L0g + v1*M1g + v2*M2g
    Dg = (g.T*adj(Lg)*g)[0] + 4*Fv*Lg.det()
    A0 = (g.T*adj(L0g)*g)[0] + 4*Fv*L0g.det() + 2*Fv*rho**2*(M1g.det() + M2g.det())
    A1 = rho*((g.T*adj(M1g)*g)[0] + 4*Fv*sig(L0g, M1g))
    B1 = rho*((g.T*adj(M2g)*g)[0] + 4*Fv*sig(L0g, M2g))
    A2 = 2*Fv*rho**2*(M1g.det() - M2g.det())
    B2 = 2*Fv*rho**2*sig(M1g, M2g)
    Dh = A0 + A1*sp.cos(phi) + B1*sp.sin(phi) + A2*sp.cos(2*phi) + B2*sp.sin(2*phi)
    r = _iszero(sp.expand_trig(sp.expand(Dg - Dh)))
    print('R6_HARMONIC_GENERIC_RES    =', r); ok &= (r == 0)

    # --- Lemma 4: structure-matrix amplitude ---
    r = _iszero(M1.det() + M2.det() - 32*s_)
    print('R7_detM_sum_RESIDUAL       =', r); ok &= (r == 0)
    r = _iszero((M1.det() - M2.det())**2 + sig(M1, M2)**2 - (64*s_)**2)
    print('R8_amp2_RESIDUAL           =', r); ok &= (r == 0)

    print('PASS_C82_FIBER_ELIMINATION_AUDIT' if ok else 'FAIL_C82_FIBER_ELIMINATION_AUDIT')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())

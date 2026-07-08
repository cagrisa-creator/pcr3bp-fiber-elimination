#!/usr/bin/env python3
"""
Symbolic + numerical audit of the closed form for the regularized rotation
term in the doubly-covered elliptic coordinates of Kim (arXiv:1701.07258),
equation (rotation) of the paper (Section "Limitations and outlook"):

    (cosh^2 lambda - cos^2 nu) * (q1 p2 - q2 p1)
        = (1/2) [ p_lambda sin(2 nu) + p_nu sinh(2 lambda) ].

The elliptic coordinates are (eq. (3) of Kim)
    q1 = (1/2) cosh(lambda) cos(nu),  q2 = (1/2) sinh(lambda) sin(nu),
with canonical momenta p_lambda, p_nu determined by
    p_lambda d(lambda) + p_nu d(nu) = p1 dq1 + p2 dq2.

This identity is the exact fact quoted in the outlook: the rotation term is
LINEAR in the momenta, hence Q_c^a = Q1(lambda,p_lambda) + Q2(nu,p_nu)
+ a*(rotation) is non-mechanical for any a>0, so Kim's curvature reduction
(valid only for the mechanical case a=0) does not carry over.

Prints PASS_ROTATION_TERM_AUDIT on success.

Author: Cesar Augusto Grisa (UFF, Brazil).
"""
import numpy as np
import sympy as sp


def symbolic_identity():
    lam, nu, p1, p2 = sp.symbols('lambda nu p1 p2', real=True)
    q1 = sp.Rational(1, 2)*sp.cosh(lam)*sp.cos(nu)
    q2 = sp.Rational(1, 2)*sp.sinh(lam)*sp.sin(nu)
    # canonical momenta p_lambda, p_nu in terms of p1, p2
    J = sp.Matrix([[sp.diff(q1, lam), sp.diff(q2, lam)],
                   [sp.diff(q1, nu),  sp.diff(q2, nu)]])
    plam, pnu = (J*sp.Matrix([p1, p2]))
    L = sp.expand_trig(sp.simplify(q1*p2 - q2*p1))
    denom = sp.cosh(lam)**2 - sp.cos(nu)**2
    lhs = sp.simplify(sp.expand_trig(L*denom))
    rhs = sp.Rational(1, 2)*(plam*sp.sin(2*nu) + pnu*sp.sinh(2*lam))
    residual = sp.simplify(sp.expand_trig(lhs - rhs))
    return residual


def numerical_audit(n=200000, seed=20260703):
    rng = np.random.default_rng(seed)
    max_err = 0.0
    for _ in range(n):
        lam = rng.uniform(-2, 2)
        nu = rng.uniform(-np.pi, np.pi)
        p1 = rng.uniform(-3, 3)
        p2 = rng.uniform(-3, 3)
        q1 = 0.5*np.cosh(lam)*np.cos(nu)
        q2 = 0.5*np.sinh(lam)*np.sin(nu)
        plam = p1*(0.5*np.sinh(lam)*np.cos(nu)) + p2*(0.5*np.cosh(lam)*np.sin(nu))
        pnu = p1*(-0.5*np.cosh(lam)*np.sin(nu)) + p2*(0.5*np.sinh(lam)*np.cos(nu))
        lhs = (np.cosh(lam)**2 - np.cos(nu)**2)*(q1*p2 - q2*p1)
        rhs = 0.5*(plam*np.sin(2*nu) + pnu*np.sinh(2*lam))
        max_err = max(max_err, abs(lhs - rhs))
    return max_err


def main():
    res = symbolic_identity()
    print("symbolic residual (should be 0):", res)
    err = numerical_audit()
    print("max |lhs - rhs| over 200000 random points: %.3e" % err)
    if res == 0 and err < 1e-9:
        print("PASS_ROTATION_TERM_AUDIT")
    else:
        print("FAIL")


if __name__ == '__main__':
    main()

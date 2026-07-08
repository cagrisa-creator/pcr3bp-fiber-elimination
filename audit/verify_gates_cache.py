#!/usr/bin/env python3
"""Verifies _gates_cache.py against the FULL symbolic derivation: rebuilds the
six gate expressions from scratch and checks exact zero residuals versus the
cached srepr forms. Run from audit/:  python3 verify_gates_cache.py
(This is the slow path by design; ~1-2 min depending on SymPy version.)"""
import sys, os, sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
from _gates_cache import CACHE
z1, z2, mu, c = sp.symbols('z1 z2 mu c', real=True)
s_ = z1**2 + z2**2; a_ = z1**2 - z2**2
R_ = sp.sqrt(1 - 4*a_ + 4*s_**2)
F = 1 - mu + 2*mu*s_/R_ - 2*c*s_ + mu**2*s_ + 4*s_**3 - 4*mu*s_*a_
Fz = sp.Matrix([sp.diff(F, z1), sp.diff(F, z2)])
Fzz = sp.hessian(F, (z1, z2))
M1 = sp.Matrix([[-4*z2, -4*z1], [-4*z1, -12*z2]])
M2 = sp.Matrix([[12*z1, 4*z2], [4*z2, 4*z1]])
L0 = -sp.Rational(1, 2)*Fzz
adj = lambda M: sp.Matrix([[M[1, 1], -M[0, 1]], [-M[1, 0], M[0, 0]]])
sig = lambda X, Y: sp.trace(X)*sp.trace(Y) - sp.trace(X*Y)
D = {
 'calT': (Fz.T*Fz)[0] - 2*F*sp.trace(Fzz) - 64*F*sp.sqrt(F)*sp.sqrt(s_),
 'A0e': (Fz.T*adj(L0)*Fz)[0] + 4*F*L0.det() + 64*F**2*s_,
 'A1e': sp.sqrt(F)*((Fz.T*adj(M1)*Fz)[0] + 4*F*sig(L0, M1)),
 'B1e': sp.sqrt(F)*((Fz.T*adj(M2)*Fz)[0] + 4*F*sig(L0, M2)),
 'A2e': -128*F**2*a_, 'B2e': -256*F**2*z1*z2}
ok = True
for k, e in D.items():
    r = sp.expand(sp.sympify(CACHE[k]) - e)
    if r != 0:
        r = sp.simplify(r)
    print("%-5s residual = %s" % (k, r))
    ok &= (r == 0)
print("PASS_GATES_CACHE" if ok else "FAIL_GATES_CACHE")

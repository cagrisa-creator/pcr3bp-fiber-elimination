#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v1.3.1 validation artifact: re-derives, by EXHAUSTIVE base x fiber scans
(every grid point, exact degree-two fiber polynomial, NO worst-point
heuristic), the sign and magnitude of every dedicated boundary and spot
value quoted in the paper. Run from the audit/ directory:
    python3 v130_validation.py
Resolution: n_th=400, n_s=160, n_phi=512 (documented; boundary claims are
stated at this resolution). Expected total runtime ~3 min.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
import numpy as np
from convexity_map import fA0, fA1, fB1, fA2, fB2, fF, cL1_of

def D_full_min(muf, cc, n_th=400, n_s=160, n_phi=512):
    th = np.linspace(0, np.pi/2, n_th); qq = np.cos(2*th)
    hi = np.where(np.abs(qq - 1) < 0.02, 0.499, 1.0)
    Sg = np.linspace(1e-9, 1, 4000)[:, None]*hi[None, :]
    V = fF(Sg, qq[None, :], muf, cc)
    ch = np.diff(np.sign(V), axis=0) != 0
    i0 = ch.argmax(axis=0); ar = np.arange(n_th)
    lo, hi2 = Sg[i0, ar], Sg[i0+1, ar]
    for _ in range(55):
        mid = (lo+hi2)/2
        sw = np.sign(fF(mid, qq, muf, cc)) == np.sign(fF(lo, qq, muf, cc))
        lo, hi2 = np.where(sw, mid, lo), np.where(sw, hi2, mid)
    r0 = (lo+hi2)/2
    sv = (np.arange(1, n_s+1)/n_s)[:, None]*r0[None, :]
    r = np.sqrt(sv); X = r*np.cos(th)[None, :]; Y = r*np.sin(th)[None, :]
    with np.errstate(invalid='ignore'):
        A0 = fA0(X, Y, muf, cc); A1 = fA1(X, Y, muf, cc); B1 = fB1(X, Y, muf, cc)
        A2 = fA2(X, Y, muf, cc); B2 = fB2(X, Y, muf, cc)
    ph = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
    c1, s1, c2, s2 = np.cos(ph), np.sin(ph), np.cos(2*ph), np.sin(2*ph)
    Dm = np.full(A0.shape, np.inf)
    for k in range(0, n_phi, 64):
        blk = (A0[..., None] + A1[..., None]*c1[k:k+64] + B1[..., None]*s1[k:k+64]
               + A2[..., None]*c2[k:k+64] + B2[..., None]*s2[k:k+64])
        Dm = np.minimum(Dm, np.nanmin(blk, axis=-1))
    return float(np.nanmin(Dm))

CHECKS = [
    # (label, mu_formula, delta_c, expected_sign)   sign: -1 non-convex, +1 convex
    ("CSV cross-check mu=0.05  dc=1e-4 (was false-convex in v1)", 0.05,      1e-4, -1),
    ("CSV cross-check mu=0.9   dc=1e-3",                          0.9,       1e-3, +1),
    ("CSV cross-check mu=0.97  dc=1e-4",                          0.97,      1e-4, +1),
    ("Light boundary: m_r=0.05 side (negative)",                  0.95,      1e-4, -1),
    ("Light boundary: m_r=0.04 side (positive)",                  0.96,      1e-4, +1),
    ("EARTH component m_r=0.9878, dc=1e-4 (macroscopic lobe)",    0.012150585609624, 1e-4,   -1),
    ("EARTH recovery bracket low  dc=1.1e-2 (negative)",          0.012150585609624, 1.1e-2, -1),
    ("EARTH recovery bracket high dc=1.2e-2 (positive)",          0.012150585609624, 1.2e-2, +1),
    ("MOON component m_r=0.0122, dc=1e-4 (convex; tiny margin)",  0.987849414390376, 1e-4,   +1),
]

ok_all = True
lines = ["v1.3.1 exhaustive validation (n_th=400, n_s=160, n_phi=512)", ""]
for label, mf, dc, expect in CHECKS:
    d = D_full_min(mf, cL1_of(mf) + dc)
    got = -1 if d < 0 else +1
    ok = (got == expect); ok_all &= ok
    line = "%-58s D_min=%+.4e  %s" % (label, d, "ok" if ok else "FAIL")
    print(line); lines.append(line)
verdict = "PASS_V130_VALIDATION" if ok_all else "FAIL_V130_VALIDATION"
print(verdict); lines += ["", verdict]
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'v130_validation_output.txt')
open(out, 'w').write("\n".join(lines) + "\n")
print("wrote", os.path.relpath(out))

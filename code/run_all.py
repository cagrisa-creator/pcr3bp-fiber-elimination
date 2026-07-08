#!/usr/bin/env python3
"""
Reproducibility driver.

Runs the full pipeline and checks the key reproducible claims:
  1. The exact symbolic audit (all eight residuals must be zero).
  2. The critical-energy sanity check c_L1(equal masses) = 2.000000.
  3. A spot reconnaissance at three (mu, c) points, checking the sign of the
     true determinant minimum against the paper's phase map:
       - Earth--Moon MOON component (mu=0.9878, m_r=0.0122): convex (D_true > 0)
       - Earth--Moon EARTH component (mu=0.0122, m_r=0.9878), near critical:
         non-convex (D_true < 0)  [new in v1.3.0 -- was a false 'convex' before]
       - equal masses, near critical  : non-convex (D_true < 0)
       - unequal masses, larger energy: convex   (D_true > 0)
  4. The fine structure of the non-convexity lobes (Observation 6.2): across
     the comparable-mass band the reduced tangential gate calT stays positive
     inside every lobe (single simple fold), and the collar radial derivative
     dF/dr stays negative (fiberwise star-shapedness).
  5. The closed form for the regularized rotation term in Kim's doubly-covered
     elliptic coordinates (eq. (rotation) of the paper), checked symbolically
     and numerically; the term is linear in the momenta, which is the exact
     obstruction to extending Kim's mechanical-case proof to the rotating
     problem.

Exits 0 iff every check passes. Run from the code/ directory:
    python run_all.py

Full map regeneration (slower, ~2 min) is a separate step:
    python convexity_map.py

Author: Cesar Augusto Grisa (UFF, Brazil).
"""
import sys
import subprocess


def check_audit():
    print("[1/5] Exact symbolic audit ...")
    r = subprocess.run([sys.executable, 'fiber_elimination_audit.py'],
                       capture_output=True, text=True)
    print(r.stdout.rstrip())
    ok = 'PASS_C82_FIBER_ELIMINATION_AUDIT' in r.stdout
    print("      -> %s" % ('PASS' if ok else 'FAIL'))
    return ok


def check_map():
    print("[2/5] Critical-energy sanity + spot reconnaissance ...")
    from convexity_map import cL1_of, scan
    ok = True

    c_eq = cL1_of(0.5)
    print("      c_L1(equal masses) = %.6f (expected 2.000000)" % c_eq)
    ok &= abs(c_eq - 2.0) < 1e-6

    # (mu, delta_c, expected sign of D_true_min): sign +1 convex, -1 non-convex.
    # mu = distant-primary mass; the regularized component has m_r = 1-mu.
    cases = [
        (0.987849414390376, 1e-3, +1, 'MOON component (Earth--Moon)'),
        (0.012150585609624, 1e-4, -1, 'EARTH component, near critical'),
        (0.012150585609624, 2e-2, +1, 'EARTH component, dc=0.02'),
        (0.50, 1e-3, -1, 'equal masses, near critical'),
        (0.20, 0.30, +1, 'unequal masses, larger energy'),
    ]
    for muf, dc, expect, label in cases:
        out = scan(muf, cL1_of(muf) + dc)
        if out is None:
            print("      %-34s : NO ROOT (fail)" % label)
            ok = False
            continue
        _, _, dtrue = out
        got = +1 if dtrue > 0 else -1
        status = 'ok' if got == expect else 'FAIL'
        print("      %-34s : D_true_min=%+.3e  %s" % (label, dtrue, status))
        ok &= (got == expect)

    print("      -> %s" % ('PASS' if ok else 'FAIL'))
    return ok


def check_relevance():
    print("[3/5] Dynamical-relevance observation (on-axis check) ...")
    import numpy as np
    from dynamical_relevance import Dmin, sstar
    from convexity_map import cL1_of
    mf, cc = 0.5, cL1_of(0.5) + 1e-3
    ok = True
    for th, name in [(0.0, 'retrograde axis (z2=0)'), (np.pi/2, 'direct axis (z1=0)')]:
        r0 = sstar(th, mf, cc)
        vals = [Dmin(np.sqrt(fr*r0)*np.cos(th), np.sqrt(fr*r0)*np.sin(th), mf, cc)
                for fr in np.linspace(0.1, 0.999, 40)]
        mn = min(vals)
        print("      min D along %-22s = %+.4e  (positive: %s)" % (name, mn, mn > 0))
        ok &= (mn > 0)
    print("      (figure: run `python dynamical_relevance.py` to regenerate)")
    print("      -> %s" % ('PASS' if ok else 'FAIL'))
    return ok


def check_lobe_structure():
    print("[4/5] Fine structure of the non-convexity lobes ...")
    from lobe_structure import analyze
    ok_fold = True
    ok_star = True
    for mf in (0.20, 0.50, 0.80):
        r = analyze(mf, n_theta=120, n_rad=60)
        band = ("[%.1f, %.1f]" % r['band']) if r['band'] else "none"
        print("      mu=%.2f lobe band=%-14s min calT=%+.3e  max dF/dr=%+.3e"
              % (mf, band, r['minT'], r['max_dFdr']))
        ok_fold &= (r['minT'] is not None and r['minT'] > 0)
        ok_star &= (r['max_dFdr'] < 0)
    print("      simple fold (calT>0 in lobes): %s" % ('ok' if ok_fold else 'FAIL'))
    print("      radial transversality (dF/dr<0): %s" % ('ok' if ok_star else 'FAIL'))
    ok = ok_fold and ok_star
    print("      -> %s" % ('PASS' if ok else 'FAIL'))
    return ok


def check_rotation_term():
    print("[5/5] Rotation-term closed form (elliptic coordinates) ...")
    r = subprocess.run([sys.executable, 'rotation_term_audit.py'],
                       capture_output=True, text=True)
    print("     ", r.stdout.rstrip().replace("\n", "\n      "))
    ok = 'PASS_ROTATION_TERM_AUDIT' in r.stdout
    print("      -> %s" % ('PASS' if ok else 'FAIL'))
    return ok


def main():
    print("=" * 66)
    print("PCR3BP fiber-elimination: reproducibility driver")
    print("=" * 66)
    results = [check_audit(), check_map(), check_relevance(),
               check_lobe_structure(), check_rotation_term()]
    print("=" * 66)
    if all(results):
        print("ALL CHECKS PASSED")
        return 0
    print("SOME CHECKS FAILED")
    return 1


if __name__ == '__main__':
    raise SystemExit(main())

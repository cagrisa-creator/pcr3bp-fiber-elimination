"""Exact independent reconstruction of the scalar convexity criterion.

This is an algebraic audit; it supplies no numerical positivity certificate.
"""
import json
from pathlib import Path
import sympy as sp

x, y, mu, c = sp.symbols('x y mu c', real=True)
v1, v2, eta1, eta2 = sp.symbols('v1 v2 eta1 eta2', real=True)
xi1, xi2, u1, u2 = sp.symbols('xi1 xi2 u1 u2', real=True)
g1, g2, h11, h12, h22 = sp.symbols('g1 g2 h11 h12 h22', real=True)
z = sp.Matrix([x, y])
s = x*x + y*y
a = x*x - y*y
q = sp.Matrix([2*a, -4*x*y])
R = sp.Matrix([[x, -y], [-y, -x]])
Q = q.jacobian(z)
A = sp.Matrix([q[1], mu-q[0]])
B = sp.Matrix([-(2*s+mu)*y, (2*s-mu)*x])
M1 = sp.Matrix([[-4*y, -4*x], [-4*x, -12*y]])
M2 = sp.Matrix([[12*x, 4*y], [4*y, 4*x]])
v = sp.Matrix([v1, v2])
eta = sp.Matrix([eta1, eta2])
xi = sp.Matrix([xi1, xi2])
u = sp.Matrix([u1, u2])
g = sp.Matrix([g1,g2])
HF = sp.Matrix([[h11,h12],[h12,h22]])
J = sp.Matrix([[0,-1],[1,0]])

checks = {}
def assert_zero(name, expr):
    expressions = list(expr) if isinstance(expr, sp.MatrixBase) else [expr]
    residuals = [sp.factor(sp.together(e)) for e in expressions]
    assert all(r == 0 for r in residuals), (name, residuals)
    checks[name] = 'PASS'

assert_zero('LC_jacobian', Q-4*R)
assert_zero('LC_norm', R.T*R-s*sp.eye(2))
assert_zero('completed_square_shift', R*A-B)
assert_zero('M1_is_Hessian_B1', sp.hessian(B[0],z)-M1)
assert_zero('M2_is_Hessian_B2', sp.hessian(B[1],z)-M2)

# w=4v is the canonical LC momentum.  Away from s=0, p=R*v/s.
p = R*v/s
assert_zero('canonical_LC_momentum', Q.T*p-4*v)
assert_zero('regularized_kinetic_square', s*((p+A).dot(p+A))-(v+B).dot(v+B))

dist = sp.sqrt(1-4*a+4*s*s)
F = 1-mu+2*mu*s/dist-2*c*s+mu*mu*s+4*s**3-4*mu*s*a
r = 2*s
Omega = (q-sp.Matrix([mu,0])).dot(q-sp.Matrix([mu,0]))/2+(1-mu)/r+mu/dist
assert_zero('F_is_r_times_Omega_minus_c', F-r*(Omega-c))

# Use arbitrary F second-order data, so no numerical differentiation is involved.
P = (v+B).dot(v+B)
PB = sp.hessian(P, [x,y,v1,v2])/2
PB[:2,:2] -= HF/2
dv = u-B.jacobian(z)*xi
variation = xi.col_join(dv)
qform = (variation.T*PB*variation)[0]
L = -HF/2+eta1*M1+eta2*M2
qform_eta = qform.subs({v1:eta1-B[0],v2:eta2-B[1]}, simultaneous=True)
assert_zero('half_Hessian_quadratic_form', qform_eta-u.dot(u)-(xi.T*L*xi)[0])

f = sp.symbols('f', positive=True)
L11,L12,L22=sp.symbols('L11 L12 L22',real=True)
Lg=sp.Matrix([[L11,L12],[L12,L22]])
C=Lg+g*g.T/(4*f)
T=g.dot(g)+4*f*sp.trace(Lg)
D=(g.T*Lg.adjugate()*g)[0]+4*f*Lg.det()
assert_zero('T_is_4F_trace_C', T-4*f*sp.trace(C))
assert_zero('D_is_4F_det_C', D-4*f*C.det())

theta,t=sp.symbols('theta t',real=True)
e=sp.Matrix([sp.cos(theta),sp.sin(theta)])
uperp=(g.dot(xi))/(2*sp.sqrt(f))*e+t*J*e
assert_zero('interior_tangent_constraint', sp.trigsimp(2*sp.sqrt(f)*e.dot(uperp)-g.dot(xi)))
assert_zero('interior_restricted_Hessian', sp.trigsimp(uperp.dot(uperp)+(xi.T*Lg*xi)[0]-t*t-(xi.T*C*xi)[0]))

assert_zero('boundary_adjugate_identity', (g.T*Lg.adjugate()*g)[0]-((J*g).T*Lg*(J*g))[0])
assert_zero('boundary_tangent_constraint', g.dot(J*g))

# Both coordinate reflections lift to orthogonal symmetries of the ambient
# regularized surface.  They are not quarter-turn rotational symmetries.
for axis, ref in [('z1',sp.diag(-1,1)),('z2',sp.diag(1,-1))]:
    substitutions={x:(ref*z)[0],y:(ref*z)[1]}
    reflected_B=B.subs(substitutions,simultaneous=True)
    assert_zero('reflection_'+axis+'_shift', reflected_B+ref*B)
    reflected_v=-ref*v
    assert_zero('reflection_'+axis+'_kinetic',
                (reflected_v+reflected_B).dot(reflected_v+reflected_B)-P)
    assert_zero('reflection_'+axis+'_F',F.subs(substitutions,simultaneous=True)-F)
    reflected_eta=-ref*eta
    reflected_L=-(ref*HF*ref)/2+reflected_eta[0]*M1.subs(substitutions,simultaneous=True)+reflected_eta[1]*M2.subs(substitutions,simultaneous=True)
    assert_zero('reflection_'+axis+'_L_conjugacy',reflected_L-ref*L*ref)

result={
 'status':'PASS_GEOMETRIC_CRITERION_RECONSTRUCTION',
 'kind':'exact symbolic identities, no interval positivity certificate',
 'sympy_version':sp.__version__,
 'checks':checks,
 'total_checks':len(checks),
 'conventions':{
    'physical_H':'0.5*|p|^2 + q2*p1 - (q1-mu)*p2 - (1-mu)/|q| - mu/|q-e|',
    'LC_map':'q=(2*(z1^2-z2^2),-4*z1*z2)',
    'canonical_momentum':'w=(Dq)^T p',
    'scaled_momentum':'v=w/4; constant invertible linear scaling preserves convexity',
    'regularized_K':'K=|v+B(z)|^2-F(z)=|q|*(H+c)',
    'B':'(-(2*s+mu)*z2,(2*s-mu)*z1)'
 }
}
print(json.dumps(result,indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

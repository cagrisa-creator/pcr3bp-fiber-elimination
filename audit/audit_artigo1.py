#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUDITORIA — Artigo 1 (pcr3bp-fiber-elimination): convenção de μ e lóbulos perdidos
Autossuficiente: só precisa de sympy, numpy, scipy. Não depende do repo.
Rode:  python3 audit_artigo1.py
Cada bloco imprime PASS/FAIL. A tese auditada:
  (T1) Na eq.(2.1) do paper, μ é a massa do primário DISTANTE; o regularizado tem massa 1−μ.
  (T2) O phase map publicado tem entradas falso-convexas (heurística dos piores pontos).
"""
import sympy as sp
import numpy as np
from scipy.optimize import brentq
import json

report = {}

# ═══════════════════════════════════════════════════════════════════
# BLOCO A — DERIVAÇÃO SIMBÓLICA DE F COM MASSAS ETIQUETADAS
# ═══════════════════════════════════════════════════════════════════
# Física: frame rotante, primário REGULARIZADO de massa alpha na ORIGEM,
# primário distante de massa beta=1−alpha em e=(1,0). Baricentro em (beta,0)
# (pois alpha·0 + beta·1 = beta). Potencial efetivo:
#     Omega(q) = ½|q−b|² + alpha/|q| + beta/|q−e|,   b=(beta,0)
# Região de Hill {Omega ≥ c}; função definidora regularizada (multiplicar por |q|):
#     G(q) = |q|·(Omega(q) − c)
# LC: q = 2 z̄²  ⇒  |q|=2s, q1=2a, |q−e|² = 4s²−4a+1.
# TESTE: G ≡ F da eq.(2.1) com beta=μ (⇒ alpha=1−μ)?
print("═"*70)
print("BLOCO A — derivação simbólica de F (identificação das massas)")
print("═"*70)
z1, z2, mu, c = sp.symbols('z1 z2 mu c', real=True)
s = z1**2 + z2**2
a = z1**2 - z2**2
# q = 2 conj(z)^2 = (2(z1²−z2²), −4 z1 z2)
q1, q2 = 2*a, -4*z1*z2
alpha, beta = 1 - mu, mu           # hipótese T1: regularizado = 1−μ
b1 = beta                          # baricentro em (beta, 0)
absq  = sp.sqrt(q1**2 + q2**2)     # = 2s
dist2 = sp.sqrt((q1-1)**2 + q2**2) # distância ao primário distante
Omega = sp.Rational(1,2)*((q1-b1)**2 + q2**2) + alpha/absq + beta/dist2
G = sp.simplify(absq*(Omega - c))
F_paper = 1 - mu + 2*mu*s/sp.sqrt(1 - 4*a + 4*s**2) - 2*c*s + mu**2*s + 4*s**3 - 4*mu*s*a
resid = sp.simplify(sp.radsimp(G - F_paper))
# simplificação com s=|z|²>0: substituir amostras racionais exatas como fallback
ok_A = (resid == 0)
if not ok_A:  # verificação exata por amostragem racional (12 pontos, aritmética exata)
    import itertools, random
    random.seed(7); ok_A = True
    for _ in range(12):
        sub = {z1: sp.Rational(random.randint(1,9), random.randint(10,19)),
               z2: sp.Rational(random.randint(1,9), random.randint(10,19)),
               mu: sp.Rational(random.randint(1,9), 10),
               c:  sp.Rational(random.randint(1,40), 10)}
        if sp.simplify((G - F_paper).subs(sub)) != 0:
            ok_A = False; break
print(f"G(q)=|q|(Omega−c) com regularizado=1−μ, distante=μ  ≡  F eq.(2.1)?  "
      f"{'PASS (identidade exata)' if ok_A else 'FAIL'}")
# contra-prova: a hipótese OPOSTA (regularizado=μ) NÃO deve reproduzir F
alpha2, beta2 = mu, 1 - mu
Omega2 = sp.Rational(1,2)*((q1-beta2)**2 + q2**2) + alpha2/absq + beta2/dist2
G2 = absq*(Omega2 - c)
sub = {z1: sp.Rational(1,3), z2: sp.Rational(1,5), mu: sp.Rational(3,10), c: sp.Rational(9,5)}
ok_A2 = sp.simplify((G2 - F_paper).subs(sub)) != 0
print(f"Contra-prova (regularizado=μ) difere de F?                        "
      f"{'PASS (difere, como esperado)' if ok_A2 else 'FAIL'}")
report['A_derivacao'] = dict(identidade=bool(ok_A), contraprova=bool(ok_A2))

# ═══════════════════════════════════════════════════════════════════
# Setup numérico comum (fórmulas idênticas às do repo, re-declaradas aqui)
# ═══════════════════════════════════════════════════════════════════
Fz  = sp.Matrix([sp.diff(F_paper, z1), sp.diff(F_paper, z2)])
Fzz = sp.hessian(F_paper, (z1, z2))
M1  = sp.Matrix([[-4*z2, -4*z1], [-4*z1, -12*z2]])
M2  = sp.Matrix([[12*z1, 4*z2], [4*z2, 4*z1]])
L0  = -sp.Rational(1,2)*Fzz
adjm = lambda M: sp.Matrix([[M[1,1], -M[0,1]], [-M[1,0], M[0,0]]])
sig  = lambda X, Y: sp.trace(X)*sp.trace(Y) - sp.trace(X*Y)
A0e = (Fz.T*adjm(L0)*Fz)[0] + 4*F_paper*L0.det() + 64*F_paper**2*s
A1e = sp.sqrt(F_paper)*((Fz.T*adjm(M1)*Fz)[0] + 4*F_paper*sig(L0, M1))
B1e = sp.sqrt(F_paper)*((Fz.T*adjm(M2)*Fz)[0] + 4*F_paper*sig(L0, M2))
A2e = -128*F_paper**2*a
B2e = -256*F_paper**2*z1*z2
args = (z1, z2, mu, c)
fA0, fA1, fB1, fA2, fB2 = (sp.lambdify(args, e, 'numpy') for e in (A0e, A1e, B1e, A2e, B2e))
S, Q = sp.symbols('S Q', real=True)
fF = sp.lambdify((S, Q, mu, c),
     1 - mu + 2*mu*S/sp.sqrt(1 - 4*S*Q + 4*S**2) - 2*c*S + mu**2*S + 4*S**3 - 4*mu*S**2*Q, 'numpy')

def cL1_of(muf):
    eps = float(min(muf, 1 - muf))
    dOm = lambda x: x - (1 - eps)/(x + eps)**2 + eps/((1 - eps) - x)**2
    xL1 = brentq(dOm, -eps + 1e-9, (1 - eps) - 1e-9, xtol=1e-15)
    return float(xL1**2/2 + (1 - eps)/(xL1 + eps) + eps/((1 - eps) - xL1))

def collar(muf, cc, n_th):
    th = np.linspace(0, np.pi/2, n_th); qq = np.cos(2*th)
    hi = np.where(np.abs(qq-1) < 0.02, 0.499, 1.0)  # mesma guarda do repo
    Sg = np.linspace(1e-9, 1, 4000)[:, None]*hi[None, :]
    V  = fF(Sg, qq[None, :], muf, cc)
    ch = np.diff(np.sign(V), axis=0) != 0
    i0 = ch.argmax(axis=0); ar = np.arange(n_th)
    lo, hi2 = Sg[i0, ar], Sg[i0+1, ar]
    for _ in range(55):
        mid = (lo+hi2)/2
        sw  = np.sign(fF(mid, qq, muf, cc)) == np.sign(fF(lo, qq, muf, cc))
        lo, hi2 = np.where(sw, mid, lo), np.where(sw, hi2, mid)
    return th, (lo+hi2)/2

def D_full_min(muf, cc, n_th=600, n_s=240, n_phi=512, loc=False):
    """min de D sobre TODOS os pontos base × fibra — SEM heurística de piores pontos."""
    th, r0 = collar(muf, cc, n_th)
    sv = (np.arange(1, n_s+1)/n_s)[:, None]*r0[None, :]
    r  = np.sqrt(sv); X = r*np.cos(th)[None, :]; Y = r*np.sin(th)[None, :]
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
    j = np.unravel_index(np.nanargmin(Dm), Dm.shape)
    out = float(np.nanmin(Dm))
    return (out, float(np.degrees(th[j[1]])), float(sv[j]/r0[j[1]])) if loc else out

def D_heuristica(muf, cc, n_th=600, n_s=400, n_worst=120):
    """Réplica exata da heurística publicada (critical_shell.D_true_min)."""
    th, r0 = collar(muf, cc, n_th)
    sv = (np.arange(1, n_s+1)/n_s)[:, None]*r0[None, :]
    r  = np.sqrt(sv); X = r*np.cos(th)[None, :]; Y = r*np.sin(th)[None, :]
    with np.errstate(invalid='ignore'):
        A0 = fA0(X, Y, muf, cc); A1 = fA1(X, Y, muf, cc); B1 = fB1(X, Y, muf, cc)
        A2 = fA2(X, Y, muf, cc); B2 = fB2(X, Y, muf, cc)
        Db = A0 - np.sqrt(A1**2+B1**2) - np.sqrt(A2**2+B2**2)
    flat = np.argsort(np.nan_to_num(Db, nan=np.inf), axis=None)[:n_worst]
    ii, jj = np.unravel_index(flat, Db.shape)
    PHI = np.linspace(0, 2*np.pi, 1024, endpoint=False)
    Dt = (A0[ii, jj][:, None] + A1[ii, jj][:, None]*np.cos(PHI) + B1[ii, jj][:, None]*np.sin(PHI)
          + A2[ii, jj][:, None]*np.cos(2*PHI) + B2[ii, jj][:, None]*np.sin(2*PHI))
    fracs = sv[ii, jj]/r0[jj]
    return float(np.nanmin(Dt)), fracs

# ═══════════════════════════════════════════════════════════════════
# BLOCO B — TESTE GEOMÉTRICO: tamanho do componente vs raio de Hill
# ═══════════════════════════════════════════════════════════════════
print()
print("═"*70)
print("BLOCO B — tamanho físico do componente (raio de Hill)")
print("═"*70)
ok_B = True
for mf in [0.95, 0.999]:
    _, r0 = collar(mf, cL1_of(mf)+1e-3, 360)
    qmax = float(np.max(2*r0))
    m_reg_T1 = 1 - mf                       # hipótese T1
    rhill = (m_reg_T1/3)**(1/3)             # escala de Hill da massa 1−μ
    consistente = 0.5*rhill < qmax < 1.6*rhill
    ok_B &= consistente
    print(f"  μ={mf}: |q|_max={qmax:.4f}; Hill(1−μ={m_reg_T1})≈{rhill:.4f}; "
          f"Hill(μ={mf})≈{(mf/3)**(1/3):.3f}  → compatível com massa 1−μ? "
          f"{'PASS' if consistente else 'FAIL'}")
report['B_geometria'] = bool(ok_B)

# ═══════════════════════════════════════════════════════════════════
# BLOCO C — LÓBULOS PERDIDOS: varredura completa vs heurística publicada
# ═══════════════════════════════════════════════════════════════════
print()
print("═"*70)
print("BLOCO C — varredura completa (sem heurística) vs heurística publicada")
print("═"*70)
csv_publicado = {  # (mu_formula, dc) -> D_true_min do results/phase_map_data.csv
    (0.05, 1e-4): 'positivo (~+5.6)', (0.02, 1e-4): 'positivo', (0.001, 1e-4): 'positivo'}
casos = [(0.05, 1e-4), (0.02, 1e-4), (0.012151, 1e-4), (0.001, 1e-4)]
res_C = {}
for mf, dc in casos:
    cc = cL1_of(mf) + dc
    d_full, tst, sfr = D_full_min(mf, cc, loc=True)
    d_heu, fracs = D_heuristica(mf, cc)
    res_C[f"mu={mf},dc={dc}"] = dict(full=d_full, heuristica=d_heu,
                                     theta_deg=tst, s_frac=sfr,
                                     piores_pontos_interior=bool(np.max(fracs) < 0.9))
    print(f"  μ={mf:<9g} dc={dc:g}: FULL={d_full:+.4e} @θ={tst:.1f}°,s/s_bd={sfr:.3f} | "
          f"heurística={d_heu:+.4e} | piores-pontos todos interior (s/s_bd max={np.max(fracs):.2f})")
# veredicto: heurística positiva onde full é negativo = falso-convexo
falso_convexo = all(res_C[k]['full'] < 0 < res_C[k]['heuristica'] for k in
                    [f"mu={m},dc={d}" for m, d in casos[:2]])
print(f"  → Entradas falso-convexas demonstradas (μ=0.05 e 0.02, dc=1e-4)?  "
      f"{'PASS' if falso_convexo else 'FAIL'}")
# N-estabilidade do achado principal: CONVERGÊNCIA em 3 malhas
d1 = D_full_min(0.05, cL1_of(0.05)+1e-4, n_th=600,  n_s=240)
d2 = D_full_min(0.05, cL1_of(0.05)+1e-4, n_th=1200, n_s=480)
def D_zoom(muf, cc, th_lo, th_hi, s_lo=0.90, n_th=900, n_s=600, n_phi=1024):
    th = np.linspace(np.radians(th_lo), np.radians(th_hi), n_th)
    qq = np.cos(2*th)
    Sg = np.linspace(1e-9, 0.499, 4000)[:, None]*np.ones(n_th)[None, :]
    V  = fF(Sg, qq[None, :], muf, cc)
    ch = np.diff(np.sign(V), axis=0) != 0
    i0 = ch.argmax(axis=0); ar = np.arange(n_th)
    lo, hi2 = Sg[i0, ar], Sg[i0+1, ar]
    for _ in range(55):
        mid = (lo+hi2)/2
        sw  = np.sign(fF(mid, qq, muf, cc)) == np.sign(fF(lo, qq, muf, cc))
        lo, hi2 = np.where(sw, mid, lo), np.where(sw, hi2, mid)
    r0 = (lo+hi2)/2
    sv = (s_lo + (1-s_lo)*np.arange(1, n_s+1)/n_s)[:, None]*r0[None, :]
    r  = np.sqrt(sv); X = r*np.cos(th)[None, :]; Y = r*np.sin(th)[None, :]
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
d3 = D_zoom(0.05, cL1_of(0.05)+1e-4, 4.0, 8.0)
# critérios: mesmo sinal nas 3; incrementos encolhendo (Cauchy); variação total <5%
ok_N = (d1 < 0 and d2 < 0 and d3 < 0 and abs(d3-d2) < abs(d2-d1)
        and abs(d3-d1) < 0.05*abs(d3))
print(f"  N-convergência (3 malhas): {d1:+.5e} → {d2:+.5e} → {d3:+.5e}  "
      f"(Δ12={abs(d2-d1):.1e}, Δ23={abs(d3-d2):.1e})  {'PASS' if ok_N else 'FAIL'}")
report['C_lobulos'] = dict(casos=res_C, falso_convexo=bool(falso_convexo),
                           N_convergencia=[d1, d2, d3], N_pass=bool(ok_N))

# ═══════════════════════════════════════════════════════════════════
# BLOCO D — CALIBRAÇÃO: reproduzir números PUBLICADOS onde eles valem
# ═══════════════════════════════════════════════════════════════════
print()
print("═"*70)
print("BLOCO D — calibração contra números publicados (onde a heurística acerta)")
print("═"*70)
# D1: critical_shell publicado: μ=0.987849..., dc=1e-4 → D_min=+1.5585e-5
mf = 0.987849414390376
d_cal = D_full_min(mf, cL1_of(mf)+1e-4)
ok_D1 = abs(d_cal - 1.5585e-5) < 2e-7
print(f"  critical_shell publicado +1.5585e-05 | full scan: {d_cal:+.4e}  "
      f"{'PASS' if ok_D1 else 'FAIL'}")
# D2: tabela de lóbulos publicada (c=cL1+1e-3): μ=0.5 → min D=-0.3502; μ=0.35 → -0.8779
ok_D2 = True
for mf, alvo in [(0.5, -0.3502), (0.35, -0.8779)]:
    d_l = D_full_min(mf, cL1_of(mf)+1e-3)
    hit = abs(d_l - alvo) < 0.03*abs(alvo)
    ok_D2 &= hit
    print(f"  lóbulo publicado μ={mf}: {alvo:+.4f} | full scan: {d_l:+.4f}  "
          f"{'PASS' if hit else 'FAIL'}")
report['D_calibracao'] = dict(critical_shell=bool(ok_D1), lobulos=bool(ok_D2))

# ═══════════════════════════════════════════════════════════════════
# VEREDICTO
# ═══════════════════════════════════════════════════════════════════
print()
print("═"*70)
tudo = ok_A and ok_A2 and ok_B and falso_convexo and ok_N and ok_D1 and ok_D2
print("VEREDICTO GERAL:", "AUDITORIA CONFIRMA T1 e T2 (todos os blocos PASS)"
      if tudo else "AUDITORIA INCONCLUSIVA — ver blocos FAIL acima")
print("═"*70)
report['veredicto'] = bool(tudo)
json.dump(report, open('audit_report.json', 'w'), indent=2)
print("relatório salvo em audit_report.json")

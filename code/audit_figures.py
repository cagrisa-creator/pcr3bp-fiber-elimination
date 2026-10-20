"""Reproduce Figures 2--3 and audit the reflection error (floating point).

Requires the original lobe_structure.py beside this script or on PYTHONPATH.
All outputs are numerical samples, never interval certificates.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import time
import warnings

import numpy as np
import scipy
import sympy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import lobe_structure as core

HERE = Path(__file__).resolve().parent
REPORTS = HERE.parent / 'results' / 'figure_audit'
DATA = REPORTS / 'data'
FIGS = HERE.parent / 'paper' / 'figures'
for folder in (DATA, FIGS):
    folder.mkdir(parents=True, exist_ok=True)


def gate_blocks(x, y, mu, c):
    """Vectorized assembly of the same gates as core.gate_values."""
    x, y = np.broadcast_arrays(np.asarray(x), np.asarray(y))
    F, gx, gy, hxx, hxy, hyy = core._blocks(x, y, mu, c)
    s, a = x*x + y*y, x*x - y*y
    l11, l12, l22 = -.5*hxx, -.5*hxy, -.5*hyy
    def gag(m11, m12, m22):
        return m22*gx*gx - 2*m12*gx*gy + m11*gy*gy
    def sig(m11, m12, m22):
        return (l11+l22)*(m11+m22) - (l11*m11+2*l12*m12+l22*m22)
    rt = np.sqrt(np.maximum(F, 0))
    T = gx*gx+gy*gy-2*F*(hxx+hyy)-64*F*rt*np.sqrt(s)
    A0 = gag(l11,l12,l22)+4*F*(l11*l22-l12*l12)+64*F*F*s
    A1 = rt*(gag(-4*y,-4*x,-12*y)+4*F*sig(-4*y,-4*x,-12*y))
    B1 = rt*(gag(12*x,4*y,4*x)+4*F*sig(12*x,4*y,4*x))
    A2, B2 = -128*F*F*a, -256*F*F*x*y
    return F, T, A0, A1, B1, A2, B2


def sampled_min(coeffs, n_phi=192):
    """Finite fiber sampling, retaining the actual minimizing angle."""
    coeffs = np.asarray(coeffs)
    shape = coeffs.shape[1:]
    flat = coeffs.reshape(5,-1)
    values = np.empty(flat.shape[1]); angles = np.empty(flat.shape[1])
    ph = np.arange(n_phi)*(2*np.pi/n_phi)
    harmonics = np.array([np.ones(n_phi),np.cos(ph),np.sin(ph),np.cos(2*ph),np.sin(2*ph)])
    for start in range(0,flat.shape[1],4096):
        rows = flat[:,start:start+4096].T @ harmonics
        idx = np.argmin(rows,axis=1)
        values[start:start+len(rows)] = rows[np.arange(len(rows)),idx]
        angles[start:start+len(rows)] = ph[idx]
    return values.reshape(shape), angles.reshape(shape)


def stationary_min(coeffs):
    """Independent floating-point minimization via all real stationary roots.

    This polynomial-root calculation is not validated root isolation.
    """
    A0,A1,B1,A2,B2 = map(float,coeffs)
    poly = [-B1+2*B2,-2*A1+8*A2,-12*B2,-2*A1-8*A2,B1+2*B2]
    roots = np.roots(np.trim_zeros(poly,'f'))
    candidates = [0.,np.pi]
    candidates.extend(float((2*np.arctan(r.real))%(2*np.pi)) for r in roots if abs(r.imag)<1e-8)
    def value(phi):
        return A0+A1*np.cos(phi)+B1*np.sin(phi)+A2*np.cos(2*phi)+B2*np.sin(2*phi)
    values = np.array([value(p) for p in candidates]); idx=int(np.argmin(values))
    return float(values[idx]),float(candidates[idx])


def write_csv(name, fields, rows):
    with (DATA/name).open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(fields);writer.writerows(rows)


def audit_symmetry():
    mu,c=.5,2.001
    roots={}
    for name,theta in [('x',0),('y',np.pi/2)]:
        roots[name]=float(np.sqrt(core.sstar(theta,mu,c)))
    witnesses=[]
    for sx in (-1,1):
        for sy in (-1,1):
            x,y=sx*.4638427480,sy*.0569527281
            F,T,*co=gate_blocks(x,y,mu,c)
            dm,ph=sampled_min(co,4096)
            stationary,statphi=stationary_min(co)
            _,gx,gy,hxx,hxy,hyy=core._blocks(x,y,mu,c)
            g=np.array([gx,gy]);H=np.array([[hxx,hxy],[hxy,hyy]])
            M1=np.array([[-4*y,-4*x],[-4*x,-12*y]])
            M2=np.array([[12*x,4*y],[4*y,4*x]])
            L=-.5*H+np.sqrt(F)*(np.cos(ph)*M1+np.sin(ph)*M2)
            adj=np.array([[L[1,1],-L[0,1]],[-L[1,0],L[0,0]]])
            direct_D=float(g@adj@g+4*F*np.linalg.det(L))
            direct_T=float(g@g+4*F*np.trace(L))
            eigenvalues=np.linalg.eigvalsh(L+np.outer(g,g)/(4*F))
            assert abs(direct_D-float(dm))<1e-12
            witnesses.append(dict(mu=mu,c=c,x=x,y=y,F=float(F),calT=float(T),
                sampled_Dmin_4096=float(dm),phi_4096=float(ph),
                direct_matrix_D=direct_D,direct_matrix_T=direct_T,
                C_eigenvalue_negative=float(eigenvalues[0]),C_eigenvalue_positive=float(eigenvalues[1]),
                stationary_Dmin=stationary,stationary_phi=statphi))
    # Test reflection laws on nonspecial first-quadrant interior points/angles.
    rng=np.random.default_rng(20261020);max_error=[0.,0.]
    for _ in range(100):
        theta=rng.uniform(.001,np.pi/2-.001);fr=rng.uniform(.01,.99)
        r=np.sqrt(fr*core.sstar(theta,mu,c));x,y=r*np.cos(theta),r*np.sin(theta)
        phi=rng.uniform(0,2*np.pi)
        def direct(xx,yy,pp):
            _,_,A0,A1,B1,A2,B2=gate_blocks(xx,yy,mu,c)
            return float(A0+A1*np.cos(pp)+B1*np.sin(pp)+A2*np.cos(2*pp)+B2*np.sin(2*pp))
        max_error[0]=max(max_error[0],abs(direct(-x,y,phi)-direct(x,y,-phi)))
        max_error[1]=max(max_error[1],abs(direct(x,-y,phi)-direct(x,y,np.pi-phi)))
        assert np.allclose(gate_blocks(x,y,mu,c)[1:],core.gate_values(x,y,mu,c),rtol=1e-12,atol=1e-12)
    write_csv('referee_witnesses.csv',list(witnesses[0]),[list(w.values()) for w in witnesses])
    return dict(axis_radii=roots,witnesses=witnesses,reflection_law_max_abs_error=max_error,
                vectorized_vs_canonical_points=100)


def figure2():
    mu,c=.5,core.cL1_of(.5)+1e-3
    theta=np.linspace(0,2*np.pi,720,endpoint=False)
    fr=np.linspace(.5,1.,60)
    corrected_s=np.array([core.sstar(t,mu,c) for t in theta])
    old_s=np.array([core.sstar(min(abs(t%(np.pi/2)),np.pi/2-1e-6),mu,c) for t in theta])
    diagnostics={}; datasets={}
    for name,radii in [('old_quarter_turn',old_s),('corrected_reflection',corrected_s)]:
        r=np.sqrt(fr[:,None]*radii[None,:]);x=r*np.cos(theta);y=r*np.sin(theta)
        F,T,*co=gate_blocks(x,y,mu,c);dm,ph=sampled_min(co)
        # Original code's sqrt(F) makes materially exterior points NaN; exclude them.
        valid=F>=-1e-12
        negative=(dm<0)&valid
        if name.startswith('old'):
            negative &= F>=0
        byquad=[int(np.sum(negative&(x*sx>0)&(y*sy>0))) for sx,sy in [(1,1),(-1,1),(-1,-1),(1,-1)]]
        diagnostics[name]=dict(negative_point_count=int(negative.sum()),
            negative_by_quadrant_Q1_Q2_Q3_Q4=byquad,
            exterior_sample_count=int((F< -1e-12).sum()),
            sampled_min_D=float(dm[negative].min()) if negative.any() else None)
        datasets[name]=(x,y,F,T,dm,ph,negative)
    # Save every corrected base sample, including sampled minimizer angles.
    x,y,F,T,dm,ph,negative=datasets['corrected_reflection']
    rows=zip(np.broadcast_to(theta,x.shape).flat,np.broadcast_to(fr[:,None],x.shape).flat,
        x.flat,y.flat,F.flat,T.flat,dm.flat,ph.flat)
    write_csv('figure2_base_samples.csv',['theta','s_fraction','z1','z2','F','calT','sampled_Dmin_192','minimizing_phi_192'],rows)
    t_boundary=np.linspace(0,2*np.pi,400)
    ss=np.array([core.sstar(t,mu,c) for t in t_boundary]);xx=np.sqrt(ss)*np.cos(t_boundary);yy=np.sqrt(ss)*np.sin(t_boundary)
    write_csv('figure2_boundary.csv',['theta','s_boundary','z1','z2'],zip(t_boundary,ss,xx,yy))
    fig,ax=plt.subplots(1,2,figsize=(12,5.4))
    ax[0].plot(x[negative],y[negative],'s',color='#d62728',ms=2.5)
    ax[0].plot(xx,yy,'-',color='#888',lw=1,label='collar ($F=0$)')
    ax[0].plot([-.75,.75],[0,0],'-',color='#1f77b4',lw=2,label='axis $z_2=0$')
    ax[0].plot([0,0],[-.75,.75],'-',color='#2ca02c',lw=1.5,label='axis $z_1=0$')
    ax[0].plot([],[],'s',color='#d62728',ms=6,label='sampled $D_{\\min}<0$')
    ax[0].set_aspect('equal');ax[0].set_xlim(-.8,.8);ax[0].set_ylim(-.8,.8)
    ax[0].set_title('Where $D_{\\min}<0$ lives (equal masses, $c=c_{L_1}+10^{-3}$)')
    ax[0].legend(fontsize=8,loc='upper right');ax[0].set_xlabel('$z_1$');ax[0].set_ylabel('$z_2$')
    profile_rows=[];axis_minima={}
    for ang,lbl,col in [(0.,r'on-axis ($\theta=0^\circ$)','#1f77b4'),
        (np.radians(7),r'lobe ($\theta=7^\circ$)','#d62728'),
        (np.radians(45),r'diagonal ($\theta=45^\circ$)','#2ca02c')]:
        frs=np.linspace(.1,.999,60);r=np.sqrt(frs*core.sstar(ang,mu,c))
        Fp,Tp,*co=gate_blocks(r*np.cos(ang),r*np.sin(ang),mu,c);d,p=sampled_min(co)
        ax[1].plot(frs,d,'-',color=col,label=lbl)
        profile_rows.extend(zip(np.repeat(ang,len(frs)),frs,r*np.cos(ang),r*np.sin(ang),Fp,Tp,d,p))
    for th,name in [(0.,'z2=0'),(np.pi/2,'z1=0')]:
        frs=np.linspace(.1,.999,50);r=np.sqrt(frs*core.sstar(th,mu,c));_,_,*co=gate_blocks(r*np.cos(th),r*np.sin(th),mu,c)
        d,_=sampled_min(co);axis_minima[name]=float(d.min())
    write_csv('figure2_ray_profiles.csv',['theta','s_fraction','z1','z2','F','calT','sampled_Dmin_192','minimizing_phi_192'],profile_rows)
    ax[1].axhline(0,color='k',lw=.8,ls='--');ax[1].set_xlabel('radial fraction $s/s^*$ (1.0 = collar)')
    ax[1].set_ylabel('sampled determinant minimum over fiber')
    ax[1].set_title('Sampled determinant profiles along three base rays')
    ax[1].legend(fontsize=9);ax[1].grid(alpha=.3);fig.tight_layout()
    fig.savefig(FIGS/'dynamical_relevance.png',dpi=150);plt.close(fig)
    diagnostics['axis_sampled_minima_192']=axis_minima
    diagnostics['max_wrong_boundary_s_error']=float(np.max(np.abs(old_s-corrected_s)))
    diagnostics['parameters']=dict(mu=mu,c=c,n_theta=720,n_radial=60,n_phi=192,
        radial_fractions=[.5,1.],boundary_points=400,ray_profile_points=60)
    return diagnostics


def analyze_grid(mu,n_theta=200,n_radial=80,n_phi=192):
    c=core.cL1_of(mu)+1e-3;theta=np.linspace(1e-4,np.pi/2-1e-4,n_theta)
    roots=np.array([core.sstar(t,mu,c) for t in theta]);fr=np.linspace(.5,.9999,n_radial)
    r=np.sqrt(fr[:,None]*roots[None,:]);x,y=r*np.cos(theta),r*np.sin(theta)
    F,T,*co=gate_blocks(x,y,mu,c);dm,ph=sampled_min(co,n_phi);negative=dm<0
    collar_derivative=2*np.sqrt(roots)*core.fdFds(roots,np.cos(2*theta),mu,c)
    if not negative.any():
        raise RuntimeError('No sampled negative lobe at '+str(mu))
    it=np.argmin(np.where(negative,T,np.inf));id_=np.argmin(dm)
    it=np.unravel_index(it,T.shape);id_=np.unravel_index(id_,T.shape)
    bands=np.degrees(theta[np.any(negative,axis=0)])
    row=dict(mu=float(mu),m_r=float(1-mu),c=float(c),minT=float(T[it]),minD=float(dm[id_]),
        band_min_deg=float(bands.min()),band_max_deg=float(bands.max()),
        negative_min_s_fraction=float(np.broadcast_to(fr[:,None],negative.shape)[negative].min()),
        negative_max_s_fraction=float(np.broadcast_to(fr[:,None],negative.shape)[negative].max()),
        max_dFdr=float(collar_derivative.max()),
        T_witness_z1=float(x[it]),T_witness_z2=float(y[it]),T_witness_phi=float(ph[it]),
        D_witness_z1=float(x[id_]),D_witness_z2=float(y[id_]),D_witness_phi=float(ph[id_]),
        negative_sample_count=int(negative.sum()))
    return row


def figure3():
    # Original plot used thirteen masses; the five canonical log masses are a subset.
    masses=np.linspace(.2,.8,13);rows=[]
    for mu in masses:
        row=analyze_grid(mu);rows.append(row);print('lobe',json.dumps(row),flush=True)
    write_csv('figure3_mass_series.csv',list(rows[0]),[list(r.values()) for r in rows])
    mu=.5;c=core.cL1_of(mu)+1e-3;theta=np.radians(7)
    fr=np.linspace(.1,.9999,160);r=np.sqrt(fr*core.sstar(theta,mu,c));x,y=r*np.cos(theta),r*np.sin(theta)
    F,T,*co=gate_blocks(x,y,mu,c);d,ph=sampled_min(co)
    write_csv('figure3_lobe_ray.csv',['theta','s_fraction','z1','z2','F','calT','sampled_Dmin_192','minimizing_phi_192'],
        zip(np.repeat(theta,len(fr)),fr,x,y,F,T,d,ph))
    fig,ax=plt.subplots(1,2,figsize=(12,5.4))
    ax[0].plot(fr,d,color='red',label=r'$D_{\min}$ (sampled fiber minimum)')
    ax[0].plot(fr,T,color='#1f77b4',label=r'$\mathcal{T}$ (tangential gate)')
    ax[0].axhline(0,color='k',lw=.8,ls='--');ax[0].set_title(r'Lobe ray ($\theta=7^\circ$): one negative direction')
    ax[0].set_xlabel('radial fraction $s/s^*$ (1.0 = collar)');ax[0].set_ylabel('sampled gate value')
    ax[0].legend(fontsize=9);ax[0].grid(alpha=.3)
    ax[1].plot(masses,[r['minT'] for r in rows],'-o',color='#1f77b4',label=r'min $\mathcal{T}$ at sampled lobe points')
    ax[1].plot(masses,[r['minD'] for r in rows],'-s',color='#d62728',label=r'min $D_{\min}$ at sampled lobe points')
    ax[1].axhline(0,color='k',lw=.8,ls='--');ax[1].set_xlabel(r'distant-primary mass fraction $\mu$')
    ax[1].set_ylabel('sampled gate minimum over the lobes')
    ax[1].set_title(r'Across the band: $\mathcal{T}>0$ while $D_{\min}<0$')
    ax[1].legend(fontsize=8);ax[1].grid(alpha=.3);fig.tight_layout()
    fig.savefig(FIGS/'lobe_structure.png',dpi=150);plt.close(fig)
    return dict(parameters=dict(delta_c=.001,n_theta=200,n_radial=80,n_phi=192,
        theta_interval=[.0001,float(np.pi/2-.0001)],radial_fraction=[.5,.9999],
        mass_values=masses.tolist(),lobe_ray_theta_deg=7,lobe_ray_n_radial=160),rows=rows)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--skip-figure3',action='store_true');args=parser.parse_args()
    started=time.perf_counter()
    report=dict(status='FLOATING_POINT_AUDIT_NOT_INTERVAL_CERTIFICATION',
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
            sympy=sympy.__version__,matplotlib=matplotlib.__version__),
        symmetry=audit_symmetry(),figure2=figure2())
    if not args.skip_figure3:report['figure3']=figure3()
    report['elapsed_seconds']=time.perf_counter()-started
    with (REPORTS/'figure_audit_report.json').open('w',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""Reproduce the finite PCR3BP samples, exporting locations and actual gates.

All numerical outputs are float64 reconnaissance. No positivity certificate is
constructed. Historical scans preserve the worst-bound union outer-three-row
selection, angular samples, root search and NaN omission of the released code.
The boundary-complete comparison instead imposes F=0 on the numerically located
boundary and evaluates its nonsingular gate formulas there; this remains floating
point and does not enclose root errors.

Usage (installed numpy/scipy/sympy, from any directory):
 python audit_parameters.py --source ../code --reference ../results/phase_map_data.csv \
     --output ../results/parameter_audit --suite all
"""
import argparse, csv, hashlib, json, math, sys, time, warnings
from pathlib import Path
sys.dont_write_bytecode = True
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',type=Path,required=True)
p.add_argument('--reference',type=Path)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--suite',choices=['phase','critical','spot','witnesses','corrected','threshold','all'],default='all')
p.add_argument('--corrected-source',type=Path)
p.add_argument('--deps',type=Path)
p.add_argument('--export-only',action='store_true',help='Re-export existing JSON results without recomputing numerical scans.')
args=p.parse_args()
if args.deps: sys.path.insert(0,str(args.deps.resolve()))
sys.path.insert(0,str(args.source.resolve()))
import numpy as np
import scipy, sympy
import convexity_map as cm
import lobe_structure as blocks
warnings.filterwarnings('ignore',category=RuntimeWarning)
args.output.mkdir(parents=True,exist_ok=True)

def safe(v):
    if isinstance(v,dict): return {str(k):safe(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [safe(x) for x in v]
    if isinstance(v,np.ndarray): return safe(v.tolist())
    if isinstance(v,(np.floating,float)): return float(v) if math.isfinite(float(v)) else None
    if isinstance(v,np.integer): return int(v)
    if isinstance(v,np.bool_): return bool(v)
    return v

def write_json(name,data):
    (args.output/name).write_text(json.dumps(safe(data),indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

def write_csv(name,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (args.output/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(safe(rows))

def actual_at(x,y,mf,cc,phi,force_boundary=False):
    """Direct evaluation from original equations (2)-(5), not lower bound."""
    F,gx,gy,hxx,hxy,hyy=map(float,blocks._blocks(x,y,mf,cc))
    rawF=F
    if force_boundary: F=0.0
    if F<0:
        return dict(x=x,y=y,phi=phi,F=rawF,valid_in_float_base=False)
    s=x*x+y*y; rt=math.sqrt(F)
    l11=-.5*hxx+rt*(-4*y*math.cos(phi)+12*x*math.sin(phi))
    l12=-.5*hxy+rt*(-4*x*math.cos(phi)+4*y*math.sin(phi))
    l22=-.5*hyy+rt*(-12*y*math.cos(phi)+4*x*math.sin(phi))
    T=gx*gx+gy*gy+4*F*(l11+l22)
    D=l22*gx*gx-2*l12*gx*gy+l11*gy*gy+4*F*(l11*l22-l12*l12)
    calT=gx*gx+gy*gy-2*F*(hxx+hyy)-64*F*rt*math.sqrt(s)
    return dict(x=x,y=y,phi=phi,F=rawF,F_used=F,valid_in_float_base=rawF>=0,
                boundary_formula=force_boundary,calT=calT,T_actual=T,D_actual=D)

def grid(mf,cc,nth,ns,npre,nbisect):
    th=np.linspace(0,np.pi/2,nth);q=np.cos(2*th)
    guard=np.where(np.abs(q-1)<.02,.499,1.)
    sg=np.linspace(1e-9,1,npre)[:,None]*guard[None,:]
    values=cm.fF(sg,q[None,:],mf,cc)
    change=np.diff(np.sign(values),axis=0)!=0
    if not change.any(axis=0).all():
        raise ValueError(f'Unbracketed boundary mu={mf} c={cc}')
    i=change.argmax(axis=0);ar=np.arange(nth)
    lo,hi=sg[i,ar],sg[i+1,ar]
    for _ in range(nbisect):
        mid=(lo+hi)/2
        same=np.sign(cm.fF(mid,q,mf,cc))==np.sign(cm.fF(lo,q,mf,cc))
        lo,hi=np.where(same,mid,lo),np.where(same,hi,mid)
    s0=(lo+hi)/2
    sv=(np.arange(1,ns+1)/ns)[:,None]*s0[None,:]
    radius=np.sqrt(sv)
    X=radius*np.cos(th)[None,:];Y=radius*np.sin(th)[None,:]
    return th,s0,X,Y

def block_coeffs(X,Y,mf,cc,boundary_complete=False):
    """Algebraically identical cheap-block evaluation, vectorized.

    Boundary-complete mode uses exact F=0 formulas at the last sampled row;
    no assertion of rigorously enclosing that boundary is made.
    """
    F,gx,gy,hxx,hxy,hyy=blocks._blocks(X,Y,mf,cc)
    F=np.array(F,copy=True)
    if boundary_complete: F[-1,:]=0.
    rt=np.sqrt(F)
    s=X*X+Y*Y;a=X*X-Y*Y
    l11,l12,l22=-.5*hxx,-.5*hxy,-.5*hyy
    gAg=lambda a11,a12,a22:a22*gx*gx-2*a12*gx*gy+a11*gy*gy
    sig=lambda a11,a12,a22:l11*a22+l22*a11-2*l12*a12
    calT=gx*gx+gy*gy-2*F*(hxx+hyy)-64*F*rt*np.sqrt(s)
    A0=gAg(l11,l12,l22)+4*F*(l11*l22-l12*l12)+64*F*F*s
    A1=rt*(gAg(-4*Y,-4*X,-12*Y)+4*F*sig(-4*Y,-4*X,-12*Y))
    B1=rt*(gAg(12*X,4*Y,4*X)+4*F*sig(12*X,4*Y,4*X))
    A2=-128*F*F*a;B2=-256*F*F*X*Y
    return calT,A0,A1,B1,A2,B2

def scan(mf,cc,nth=192,ns=80,nphi=256,npre=2000,nbisect=45,worst=80,
         coeff_source='cache',boundary_complete=False,full=False):
    th,s0,X,Y=grid(mf,cc,nth,ns,npre,nbisect)
    rawF=blocks.fF(X,Y,mf,cc)
    if coeff_source=='cache':
        T=cm.fT(X,Y,mf,cc)
        A0,A1,B1,A2,B2=(f(X,Y,mf,cc) for f in (cm.fA0,cm.fA1,cm.fB1,cm.fA2,cm.fB2))
        if boundary_complete:
            collar=block_coeffs(X[-1:,:],Y[-1:,:],mf,cc,True)
            T[-1]=collar[0][-1]
            for target,value in zip([A0,A1,B1,A2,B2],collar[1:]):target[-1]=value[-1]
    else:
        T,A0,A1,B1,A2,B2=block_coeffs(X,Y,mf,cc,boundary_complete)
    Db=A0-np.hypot(A1,B1)-np.hypot(A2,B2)
    # Keep original sqrt form in historical cache mode, for bitwise operation order.
    if coeff_source=='cache':
        Db=A0-np.sqrt(A1**2+B1**2)-np.sqrt(A2**2+B2**2)
    if full:
        ii,jj=np.indices(Db.shape);ii=ii.ravel();jj=jj.ravel()
    else:
        flat=np.argsort(np.nan_to_num(Db,nan=np.inf),axis=None)[:worst]
        ii,jj=np.unravel_index(flat,Db.shape)
        rows=np.arange(max(0,ns-3),ns)
        ii=np.concatenate([ii,np.repeat(rows,nth)])
        jj=np.concatenate([jj,np.tile(np.arange(nth),len(rows))])
    phi=np.linspace(0,2*np.pi,nphi,endpoint=False)
    c1,s1,c2,s2=np.cos(phi),np.sin(phi),np.cos(2*phi),np.sin(2*phi)
    best=np.inf; best_candidate=best_phase=None
    # Batching preserves each polynomial's arithmetic order while bounding RAM.
    for start in range(0,len(ii),1024):
        ir,jr=ii[start:start+1024],jj[start:start+1024]
        vals=A0[ir,jr,None]+A1[ir,jr,None]*c1+B1[ir,jr,None]*s1+A2[ir,jr,None]*c2+B2[ir,jr,None]*s2
        if np.isfinite(vals).any():
            loc=np.unravel_index(np.nanargmin(vals),vals.shape)
            if vals[loc]<best:
                best=float(vals[loc]);best_candidate=start+loc[0];best_phase=loc[1]
    row,col=int(ii[best_candidate]),int(jj[best_candidate])
    witness=actual_at(float(X[row,col]),float(Y[row,col]),mf,cc,float(phi[best_phase]),boundary_complete and row==ns-1)
    witness.update(D_harmonic=best,theta=float(th[col]),s_fraction=(row+1)/ns,
                   direct_minus_harmonic=witness.get('D_actual',float('nan'))-best,
                   row=row,column=col,phi_index=int(best_phase))
    ti,tj=np.unravel_index(np.nanargmin(T),T.shape)
    tphi=math.atan2(-float(X[ti,tj]),float(Y[ti,tj]))%(2*math.pi)
    tw=actual_at(float(X[ti,tj]),float(Y[ti,tj]),mf,cc,tphi,boundary_complete and ti==ns-1)
    tw.update(calT_sampled=float(T[ti,tj]),theta=float(th[tj]),s_fraction=(ti+1)/ns)
    collarT=block_coeffs(X[-1:,:],Y[-1:,:],mf,cc,True)
    # At F=0, D equals A0 and is independent of the fiber angle.
    collar_min=float(np.min(collarT[1]))
    ret=dict(mu_formula=mf,m_r=1-mf,c=cc,c_L1=cm.cL1_of(mf),delta_c=cc-cm.cL1_of(mf),
             nth=nth,ns=ns,nphi=nphi,npre=npre,nbisect=nbisect,worst=worst,
             method='exhaustive finite base x fiber' if full else 'worst amplitude bound union outer three rows',
             coefficient_source=coeff_source,boundary_complete=boundary_complete,
             calT_min=float(np.nanmin(T)),calD_bound_min=float(np.nanmin(Db)),D_sampled_min=best,
             collar_D_F0_min=collar_min,witness=witness,tangential_minimum_sample=tw,
             base_samples=int(X.size),fiber_base_entries=int(len(ii)),
             unique_fiber_base_entries=int(len(set(zip(ii.tolist(),jj.tolist())))),
             F_boundary_min=float(np.min(rawF[-1])),F_boundary_max=float(np.max(rawF[-1])),
             F_negative_boundary_samples=int(np.count_nonzero(rawF[-1]<0)),
             F_negative_interior_samples=int(np.count_nonzero(rawF[:-1]<0)),
             invalid_T_samples=int(np.count_nonzero(~np.isfinite(T))),
             invalid_Dbound_samples=int(np.count_nonzero(~np.isfinite(Db))),
             invalid_boundary_Dbound_samples=int(np.count_nonzero(~np.isfinite(Db[-1]))),
             boundary_root_s_min=float(np.min(s0)),boundary_root_s_max=float(np.max(s0)))
    return ret

def flattened(r,label):
    return dict(label=label,**{k:v for k,v in r.items() if not isinstance(v,dict) and k!='label'},
                **{'witness_'+k:v for k,v in r['witness'].items()},
                **{'Tmin_'+k:v for k,v in r['tangential_minimum_sample'].items()})

def phase_suite():
    start=time.monotonic()
    mus=[.001,.005,.012150585609624,.02,.05,.10,.20,.35,.50,.65,.80,.90,.941,.97,.987849414390376,.999]
    offsets=[1e-4,1e-3,.01,.03,.1,.3,1.]
    ref=[]
    if args.reference:
        with args.reference.open(newline='',encoding='utf-8') as f:ref=list(csv.DictReader(f))
    results=[];comparisons=[];original_checks=[];boundary=[]
    for mf in mus:
        cl=cm.cL1_of(mf)
        for dc in offsets:
            r=scan(mf,cl+dc)
            r['nominal_delta_c']=dc
            original=cm.scan(mf,cl+dc,n_th=192,n_s=80)
            original_checks.append(dict(mu_formula=mf,delta_c=dc,max_abs_wrapper_vs_original=max(abs(a-b) for a,b in zip(original,[r['calT_min'],r['calD_bound_min'],r['D_sampled_min']]))))
            results.append(r)
            rb=scan(mf,cl+dc,coeff_source='blocks',boundary_complete=True)
            rb['nominal_delta_c']=dc;boundary.append(rb)
            if ref:
                old=next(row for row in ref if abs(float(row['mu_formula'])-mf)<1e-14 and abs(float(row['delta_c'])-dc)<1e-14)
                co=dict(mu_formula=mf,m_r=1-mf,delta_c=dc)
                for newkey,oldkey in [('calT_min','calT_min'),('calD_bound_min','calD_bound_min'),('D_sampled_min','D_true_min')]:
                    a,b=r[newkey],float(old[oldkey]);co[newkey+'_old']=b;co[newkey+'_rerun']=a;co[newkey+'_difference']=a-b
                co['D_sign_changed']=np.sign(r['D_sampled_min'])!=np.sign(float(old['D_true_min']))
                co['D_boundary_complete']=rb['D_sampled_min'];co['D_boundary_complete_minus_historical']=rb['D_sampled_min']-r['D_sampled_min']
                co['boundary_complete_sign_changed']=np.sign(rb['D_sampled_min'])!=np.sign(r['D_sampled_min'])
                comparisons.append(co)
        print('PHASE done mu',mf,flush=True)
        write_json('phase_details.json',results);write_json('phase_boundary_complete_details.json',boundary)
    write_csv('phase_map_rerun.csv',[dict(mu_formula=r['mu_formula'],m_r=r['m_r'],c_L1=r['c_L1'],delta_c=r['nominal_delta_c'],calT_min=r['calT_min'],calD_bound_min=r['calD_bound_min'],D_true_min=r['D_sampled_min']) for r in results])
    write_csv('phase_map_boundary_complete.csv',[dict(mu_formula=r['mu_formula'],m_r=r['m_r'],c_L1=r['c_L1'],delta_c=r['nominal_delta_c'],calT_min=r['calT_min'],calD_bound_min=r['calD_bound_min'],D_true_min=r['D_sampled_min']) for r in boundary])
    write_csv('phase_witnesses.csv',[flattened(r,f"phase mu={r['mu_formula']} dc={r['nominal_delta_c']}") for r in results])
    write_csv('phase_comparison.csv',comparisons)
    write_csv('phase_wrapper_original_checks.csv',original_checks)
    summary=dict(elapsed_seconds=time.monotonic()-start,pairs=len(results),numerical_type='float64; no interval certificates',
                 historical_sign_changes=[r for r in comparisons if r['D_sign_changed']],
                 boundary_complete_sign_changes=[r for r in comparisons if r['boundary_complete_sign_changed']],
                 max_wrapper_original_difference=max(r['max_abs_wrapper_vs_original'] for r in original_checks),
                 max_csv_D_abs_difference=max(abs(r['D_sampled_min_difference']) for r in comparisons) if comparisons else None,
                 max_csv_T_abs_difference=max(abs(r['calT_min_difference']) for r in comparisons) if comparisons else None,
                 max_csv_bound_abs_difference=max(abs(r['calD_bound_min_difference']) for r in comparisons) if comparisons else None,
                 historical_boundary_samples_omitted=sum(r['invalid_boundary_Dbound_samples'] for r in results),
                 negatives_at_dc_01=[r['m_r'] for r in results if r['nominal_delta_c']==.01 and r['D_sampled_min']<0])
    write_json('phase_summary.json',summary);print(json.dumps(safe(summary)),flush=True)

def critical_suite():
    start=time.monotonic();mf=.987849414390376;cl=cm.cL1_of(mf);results=[]
    common=dict(nth=600,ns=400,nphi=1024,npre=4000,nbisect=55,worst=120,coeff_source='blocks')
    for dc in [1e-3,3e-4,1e-4,3e-5,1e-5,5e-6,2e-6,1e-6,5e-7,1e-7]:
        r=scan(mf,cl+dc,**common);r['nominal_delta_c']=dc;r['label']='track';results.append(r)
        print('LUNAR track',dc,r['D_sampled_min'],flush=True)
        write_json('lunar_details.json',results)
    lo,hi=1e-6,2e-5;bis=[]
    for k in range(30):
        dc=(lo+hi)/2;r=scan(mf,cl+dc,**common)
        bis.append(dict(iteration=k,lo_before=lo,hi_before=hi,delta_c=dc,D_sampled_min=r['D_sampled_min'],witness=r['witness']))
        if r['D_sampled_min']>0:hi=dc
        else:lo=dc
        if k%5==0:print('LUNAR bisect',k,lo,hi,flush=True)
    stability=[]
    for ns in [60,150,400]:
        for nphi in [512,1024]:
            opts=dict(common,ns=ns,nphi=nphi)
            r=scan(mf,cl+3e-6,**opts);r['nominal_delta_c']=3e-6;r['label']='radial_fiber_stability';stability.append(r)
            print('LUNAR stability',ns,nphi,r['D_sampled_min'],flush=True)
    # Independent formulas at each exported minimizer are recorded by actual_at.
    write_json('lunar_details.json',results);write_json('lunar_bisection.json',bis);write_json('lunar_stability.json',stability)
    write_csv('lunar_track_witnesses.csv',[flattened(r,str(r['nominal_delta_c'])) for r in results])
    write_csv('lunar_stability_witnesses.csv',[flattened(r,f"ns={r['ns']} nphi={r['nphi']}") for r in stability])
    summary=dict(elapsed_seconds=time.monotonic()-start,c_L1=cl,nominal_sampled_shell_endpoint=(lo+hi)/2,
                 final_float_bisection_bracket=[lo,hi],bracket_is_interval_certificate=False,
                 no_limit_as_dc_tends_zero_proved=True,
                 max_track_direct_harmonic_abs_difference=max(abs(r['witness']['direct_minus_harmonic']) for r in results),
                 boundary_nan_omission_in_historical_method=True)
    write_json('lunar_summary.json',summary);print(json.dumps(safe(summary)),flush=True)

def reexport_saved():
    if args.suite in ['critical','all']:
        results=json.loads((args.output/'lunar_details.json').read_text())
        stability=json.loads((args.output/'lunar_stability.json').read_text())
        bis=json.loads((args.output/'lunar_bisection.json').read_text())
        last=bis[-1];lo=last['lo_before'];hi=last['hi_before']
        if last['D_sampled_min']>0:hi=last['delta_c']
        else:lo=last['delta_c']
        write_csv('lunar_track_witnesses.csv',[flattened(r,str(r['nominal_delta_c'])) for r in results])
        write_csv('lunar_stability_witnesses.csv',[flattened(r,f"ns={r['ns']} nphi={r['nphi']}") for r in stability])
        summary=dict(c_L1=results[0]['c_L1'],nominal_sampled_shell_endpoint=(lo+hi)/2,
                     final_float_bisection_bracket=[lo,hi],bracket_is_interval_certificate=False,
                     no_limit_as_dc_tends_zero_proved=True,
                     max_track_direct_harmonic_abs_difference=max(abs(r['witness']['direct_minus_harmonic']) for r in results),
                     boundary_nan_omission_in_historical_method=True)
        write_json('lunar_summary.json',summary);print(json.dumps(safe(summary)),flush=True)
    if args.suite in ['spot','all']:
        results=json.loads((args.output/'spot_details.json').read_text());s4=results[:4]
        write_csv('spot_witnesses.csv',[flattened(r,r['label']) for r in results])
        summary=dict(section4_grid='fresh explicit 400 angular rays x 160 radial rows; all 512 fibers for D, exact reduced T and boundary formula for collar D; F=0 enforced only on numerical boundary row; no interval arithmetic',
                     section4_minimum_sampled_T=min(r['calT_min'] for r in s4),section4_minimum_sampled_collar_D=min(r['collar_D_F0_min'] for r in s4),
                     section4_original_grid_artifact_found=False,
                     max_direct_harmonic_abs_difference=max(abs(r['witness']['direct_minus_harmonic']) for r in results))
        write_json('spot_summary.json',summary);print(json.dumps(safe(summary)),flush=True)

def spot_suite():
    start=time.monotonic();results=[]
    # No dedicated Section-4 corner-grid artifact was recovered. Declare this
    # fresh grid explicitly; it must not be attributed to the old manuscript.
    for mf in [.987839414390376,.987859414390376]:
        for cc in [1.5995,1.6005]:
            r=scan(mf,cc,nth=400,ns=160,nphi=512,npre=4000,nbisect=55,worst=120,coeff_source='blocks',boundary_complete=True,full=True)
            r['label']='section4_fresh_declared_grid';results.append(r)
            print('SECTION4',mf,cc,r['calT_min'],r['collar_D_F0_min'],flush=True)
    s4=results.copy()
    # Dedicated quoted values: original validation-artifact full finite scan.
    for label,mf,dc in [('validation_mu_005',.05,1e-4),('validation_mu_09',.9,1e-3),('validation_mu_097',.97,1e-4),('earth_nearcritical',.012150585609624,1e-4),('earth_dc_01',.012150585609624,.01),('earth_dc_011',.012150585609624,.011),('earth_dc_012',.012150585609624,.012),('light_mr_04',.96,1e-4),('light_mr_05',.95,1e-4),('moon_dc_1e4',.987849414390376,1e-4)]:
        r=scan(mf,cm.cL1_of(mf)+dc,nth=400,ns=160,nphi=512,npre=4000,nbisect=55,worst=120,coeff_source='cache',full=True)
        r['label']=label;r['nominal_delta_c']=dc;results.append(r)
        print('DEDICATED',label,r['D_sampled_min'],flush=True)
    write_json('spot_details.json',results);write_csv('spot_witnesses.csv',[flattened(r,r['label']) for r in results])
    summary=dict(elapsed_seconds=time.monotonic()-start,section4_grid='fresh explicit 400 angular rays x 160 radial rows; all 512 fibers for D, exact reduced T and boundary formula for collar D; F=0 enforced only on numerical boundary row; no interval arithmetic',
                 section4_minimum_sampled_T=min(r['calT_min'] for r in s4),section4_minimum_sampled_collar_D=min(r['collar_D_F0_min'] for r in s4),
                 section4_original_grid_artifact_found=False,
                 max_direct_harmonic_abs_difference=max(abs(r['witness']['direct_minus_harmonic']) for r in results))
    write_json('spot_summary.json',summary);print(json.dumps(safe(summary)),flush=True)

def witnesses_suite():
    """Independent 80-digit checks of actual D at negative sampled witnesses.

    This is high-precision floating point, without directed rounding or interval
    enclosure. A boundary minimizer may be moved a documented 1e-8 relative
    distance inward to provide an unambiguously interior numerical witness.
    """
    import mpmath as mp
    mp.mp.dps=80
    expr=blocks.F
    xyz=(blocks.z1,blocks.z2,blocks.mu,blocks.c)
    functions=[sympy.lambdify(xyz,e,'mpmath') for e in [expr,sympy.diff(expr,blocks.z1),sympy.diff(expr,blocks.z2),sympy.diff(expr,blocks.z1,2),sympy.diff(expr,blocks.z1,blocks.z2),sympy.diff(expr,blocks.z2,2)]]
    rows=[]
    for name in ['phase_details.json','phase_boundary_complete_details.json','phase_revised_details.json','lunar_details.json','lunar_revised_details.json','lunar_stability.json','spot_details.json','light_threshold_details.json','light_threshold_bisection_details.json']:
        if not (args.output/name).exists():continue
        records=json.loads((args.output/name).read_text())
        for index,r in enumerate(records):
            if r['D_sampled_min']>=0:continue
            w=r['witness'];mf,cc=map(lambda x:mp.mpf(repr(x)),[r['mu_formula'],r['c']])
            x0,y0=w['x'],w['y'];x,y=mp.mpf(repr(x0)),mp.mpf(repr(y0));phi=mp.mpf(repr(w['phi']))
            vals=[f(x,y,mf,cc) for f in functions]
            shifted=False
            if vals[0]<=mp.mpf('1e-20') or w.get('boundary_formula'):
                x0=float(x0*(1.-1e-8));y0=float(y0*(1.-1e-8))
                x,y=mp.mpf(repr(x0)),mp.mpf(repr(y0));vals=[f(x,y,mf,cc) for f in functions];shifted=True
            F,gx,gy,hxx,hxy,hyy=vals
            if F<0:raise AssertionError((name,index,'negative F after inward shift',str(F)))
            rt=mp.sqrt(F);co,si=mp.cos(phi),mp.sin(phi)
            l11=-hxx/2+rt*(-4*y*co+12*x*si)
            l12=-hxy/2+rt*(-4*x*co+4*y*si)
            l22=-hyy/2+rt*(-12*y*co+4*x*si)
            D=l22*gx*gx-2*l12*gx*gy+l11*gy*gy+4*F*(l11*l22-l12*l12)
            T=gx*gx+gy*gy+4*F*(l11+l22)
            calT=gx*gx+gy*gy-2*F*(hxx+hyy)-64*F*rt*mp.sqrt(x*x+y*y)
            row=dict(source=name,index=index,label=r.get('label','phase'),mu_formula=repr(r['mu_formula']),m_r=repr(r['m_r']),c=repr(r['c']),delta_c=repr(r.get('nominal_delta_c',r['delta_c'])),x=repr(x0),y=repr(y0),phi=repr(w['phi']),inward_relative_radius_shift=1e-8 if shifted else 0.,F_80digits=mp.nstr(F,60),D_actual_80digits=mp.nstr(D,60),T_actual_80digits=mp.nstr(T,60),calT_80digits=mp.nstr(calT,60),negative_D_observed=bool(D<0),calT_positive_observed=bool(calT>0),precision_decimal_digits=80,interval_certificate=False)
            if D>=0:raise AssertionError(row)
            rows.append(row)
    write_json('negative_witnesses_80digits.json',rows);write_csv('negative_witnesses_80digits.csv',rows)
    write_json('negative_witnesses_summary.json',dict(count=len(rows),all_actual_D_negative=all(r['negative_D_observed'] for r in rows),all_reduced_T_positive=all(r['calT_positive_observed'] for r in rows),precision_decimal_digits=80,interval_certificate=False,inward_shifts=sum(r['inward_relative_radius_shift']>0 for r in rows)))
    print('80-DIGIT NEGATIVE WITNESSES',len(rows),'all actual D<0; no interval certificate',flush=True)

def corrected_suite():
    import importlib.util
    corrected=args.corrected_source or Path(__file__).parent/'corrected_code'
    spec=importlib.util.spec_from_file_location('convexity_map_revised',corrected/'convexity_map.py')
    revised=importlib.util.module_from_spec(spec);spec.loader.exec_module(revised)
    oldmodule=sys.modules.get('convexity_map');sys.modules['convexity_map']=revised
    try:
        spec=importlib.util.spec_from_file_location('critical_shell_revised',corrected/'critical_shell.py')
        lunar_revised=importlib.util.module_from_spec(spec);spec.loader.exec_module(lunar_revised)
    finally:
        sys.modules['convexity_map']=oldmodule
    history=json.loads((args.output/'phase_details.json').read_text())
    frozen=[]
    if (args.output/'phase_comparison.csv').exists():
        with (args.output/'phase_comparison.csv').open(newline='',encoding='utf-8') as f:frozen=list(csv.DictReader(f))
    revised_data=[];comparison=[]
    for old in history:
        mf,cc=old['mu_formula'],old['c'];dc=old['nominal_delta_c']
        r=scan(mf,cc,coeff_source='cache',boundary_complete=True);r['nominal_delta_c']=dc
        actual=revised.scan(mf,cc,n_th=192,n_s=80)
        discrepancy=max(abs(a-b) for a,b in zip(actual,[r['calT_min'],r['calD_bound_min'],r['D_sampled_min']]))
        if discrepancy>1e-12:raise AssertionError(('corrected wrapper discrepancy',mf,dc,discrepancy))
        revised_data.append(r)
        row=dict(mu_formula=mf,m_r=1-mf,delta_c=dc,old_boundary_invalid_Dbound=old['invalid_boundary_Dbound_samples'],revised_boundary_invalid_Dbound=r['invalid_boundary_Dbound_samples'],revised_interior_invalid_Dbound=r['invalid_Dbound_samples'],wrapper_vs_actual_corrected_max_abs=discrepancy)
        for key in ['calT_min','calD_bound_min','D_sampled_min']:
            row['old_'+key]=old[key];row['revised_'+key]=r[key];row['difference_'+key]=r[key]-old[key]
        if frozen:
            oldcsv=next(v for v in frozen if abs(float(v['mu_formula'])-mf)<1e-14 and abs(float(v['delta_c'])-dc)<1e-14)
            for key in ['calT_min','calD_bound_min','D_sampled_min']:
                row['frozen_csv_'+key]=float(oldcsv[key+'_old'])
                row['revised_minus_frozen_csv_'+key]=r[key]-float(oldcsv[key+'_old'])
        row['D_sign_changed']=np.sign(old['D_sampled_min'])!=np.sign(r['D_sampled_min']);comparison.append(row)
    write_json('phase_revised_details.json',revised_data);write_csv('phase_revised_comparison.csv',comparison)
    write_csv('phase_map_data_revised.csv',[dict(mu_formula=r['mu_formula'],m_r=r['m_r'],c_L1=r['c_L1'],delta_c=r['nominal_delta_c'],calT_min=r['calT_min'],calD_bound_min=r['calD_bound_min'],D_true_min=r['D_sampled_min']) for r in revised_data])
    write_csv('phase_revised_witnesses.csv',[flattened(r,f"mu={r['mu_formula']} dc={r['nominal_delta_c']}") for r in revised_data])
    track=json.loads((args.output/'lunar_details.json').read_text());newtrack=[]
    common=dict(nth=600,ns=400,nphi=1024,npre=4000,nbisect=55,worst=120,coeff_source='cache',boundary_complete=True)
    for old in track:
        r=scan(old['mu_formula'],old['c'],**common);r['nominal_delta_c']=old['nominal_delta_c']
        actual=lunar_revised.D_true_min(old['mu_formula'],old['c'])
        r['actual_corrected_critical_shell_value']=actual
        r['actual_corrected_minus_wrapper']=actual-r['D_sampled_min']
        r['revised_minus_historical']=actual-old['D_sampled_min'];newtrack.append(r)
    mf=.987849414390376;cl=cm.cL1_of(mf);lo,hi=1e-6,2e-5;bis=[]
    for k in range(30):
        dc=(lo+hi)/2;value=lunar_revised.D_true_min(mf,cl+dc)
        bis.append(dict(iteration=k,delta_c=dc,D_sampled_min=value,lo_before=lo,hi_before=hi))
        if value>0:hi=dc
        else:lo=dc
    write_json('lunar_revised_details.json',newtrack);write_csv('lunar_revised_track_witnesses.csv',[flattened(r,str(r['nominal_delta_c'])) for r in newtrack]);write_json('lunar_revised_bisection.json',bis)
    historical_endpoint=json.loads((args.output/'lunar_summary.json').read_text())['nominal_sampled_shell_endpoint']
    summary=dict(pairs=len(revised_data),sign_changes=sum(r['D_sign_changed'] for r in comparison),boundary_NaNs_before=sum(r['old_boundary_invalid_Dbound'] for r in comparison),boundary_NaNs_after=sum(r['revised_boundary_invalid_Dbound'] for r in comparison),max_D_abs_change=max(abs(r['difference_D_sampled_min']) for r in comparison),max_T_abs_change=max(abs(r['difference_calT_min']) for r in comparison),max_Dbound_abs_change=max(abs(r['difference_calD_bound_min']) for r in comparison),max_lunar_track_abs_change=max(abs(r['revised_minus_historical']) for r in newtrack),lunar_endpoint_revised=(lo+hi)/2,lunar_endpoint_historical=historical_endpoint,lunar_endpoint_change=(lo+hi)/2-historical_endpoint,no_interval_certificate=True)
    write_json('corrected_summary.json',summary);print('CORRECTED',json.dumps(safe(summary)),flush=True)

def threshold_suite():
    start=time.monotonic();dc=1e-4
    common=dict(nth=600,ns=400,nphi=1024,npre=4000,nbisect=55,worst=120,coeff_source='cache',boundary_complete=True,full=True)
    cases=[]
    for mr in [.04,.045,.05]:
        mf=1-mr;r=scan(mf,cm.cL1_of(mf)+dc,**common)
        r['label']='light_boundary_dedicated';r['nominal_delta_c']=dc;cases.append(r)
        print('LIGHT THRESHOLD',mr,r['D_sampled_min'],flush=True)
        write_json('light_threshold_details.json',cases)
    lo,hi=.04,.05;bis=[]
    for k in range(12):
        mr=(lo+hi)/2;mf=1-mr;r=scan(mf,cm.cL1_of(mf)+dc,**common)
        r.update(label='light_boundary_mass_bisection',iteration=k,lo_before=lo,hi_before=hi,nominal_delta_c=dc)
        bis.append(r)
        # Empirical sign-transition search only, not a proof of monotonicity.
        if r['D_sampled_min']<0:hi=mr
        else:lo=mr
        print('LIGHT BISECT',k,mr,r['D_sampled_min'],flush=True)
        write_json('light_threshold_bisection_details.json',bis)
    write_csv('light_threshold_witnesses.csv',[flattened(r,str(r['m_r'])) for r in cases+bis])
    summary=dict(elapsed_seconds=time.monotonic()-start,nth=600,ns=400,nphi=1024,all_finite_base_fiber_samples=True,
                 delta_c=dc,mass_fraction_of_regularized_primary=True,empirical_transition_bracket=[lo,hi],
                 bracket_not_interval_certificate=True,continuous_mass_monotonicity_not_proved=True,
                 mr_005_D_sampled_min=cases[2]['D_sampled_min'])
    write_json('light_threshold_summary.json',summary);print('LIGHT SUMMARY',json.dumps(safe(summary)),flush=True)

if __name__=='__main__':
    if args.export_only:
        reexport_saved()
        raise SystemExit(0)
    manifest=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,sympy=sympy.__version__,
                  source_files={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in args.source.glob('*.py')},
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  numerical_scope='finite float64 sampling only; no interval arithmetic or domain-wide numerical certificate')
    write_json('environment_'+args.suite+'.json',manifest)
    if args.suite in ['phase','all']:phase_suite()
    if args.suite in ['critical','all']:critical_suite()
    if args.suite in ['spot','all']:spot_suite()
    if args.suite in ['corrected','all']:corrected_suite()
    if args.suite in ['threshold','all']:threshold_suite()
    if args.suite in ['witnesses','all']:witnesses_suite()

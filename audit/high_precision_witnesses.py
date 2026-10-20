"""Independent 80-digit evaluation of F, its derivatives, T and D from L.

No expression cache or SymPy formula is imported. Numerical differentiation
uses mpmath's high-precision derivative algorithm. This is not an interval
calculation or a rigorous rounding-error enclosure.
"""
from pathlib import Path
import csv
import json
import time
import mpmath as mp

mp.mp.dps=80
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'results/figure_audit/data/referee_witnesses.csv'
started=time.perf_counter();rows=[]

def Ffun(x,y,mu,c):
    s=x*x+y*y;a=x*x-y*y
    return 1-mu+2*mu*s/mp.sqrt(1-4*a+4*s*s)-2*c*s+mu*mu*s+4*s*s*s-4*mu*s*a

for row in csv.DictReader(SOURCE.open()):
    x,y,mu,c,phi=[mp.mpf(row[k]) for k in ('x','y','mu','c','phi_4096')]
    F=Ffun(x,y,mu,c)
    gx=mp.diff(lambda xx:Ffun(xx,y,mu,c),x)
    gy=mp.diff(lambda yy:Ffun(x,yy,mu,c),y)
    hxx=mp.diff(lambda xx:Ffun(xx,y,mu,c),x,2)
    hyy=mp.diff(lambda yy:Ffun(x,yy,mu,c),y,2)
    hxy=mp.diff(lambda xx:mp.diff(lambda yy:Ffun(xx,yy,mu,c),y),x)
    M1=mp.matrix([[-4*y,-4*x],[-4*x,-12*y]])
    M2=mp.matrix([[12*x,4*y],[4*y,4*x]])
    L=-mp.matrix([[hxx,hxy],[hxy,hyy]])/2+mp.sqrt(F)*(mp.cos(phi)*M1+mp.sin(phi)*M2)
    adj=mp.matrix([[L[1,1],-L[0,1]],[-L[1,0],L[0,0]]])
    g=mp.matrix([gx,gy])
    D=(g.T*adj*g)[0]+4*F*mp.det(L)
    T=gx*gx+gy*gy+4*F*(L[0,0]+L[1,1])
    calT=gx*gx+gy*gy-2*F*(hxx+hyy)-64*F*mp.sqrt(F)*mp.sqrt(x*x+y*y)
    C=L+g*g.T/(4*F)
    ce=mp.eigsy(C,eigvals_only=True)
    values=dict(F=F,calT=calT,T=T,D=D,C_negative=ce[0],C_positive=ce[1])
    comparison=dict(F=mp.mpf(row['F']),calT=mp.mpf(row['calT']),
                    T=mp.mpf(row['direct_matrix_T']),D=mp.mpf(row['direct_matrix_D']))
    differences={k:abs(values[k]-v) for k,v in comparison.items()}
    assert F>0 and calT>0 and T>0 and D<0 and ce[0]<0 and ce[1]>0
    rows.append(dict(inputs={k:row[k] for k in ('x','y','mu','c','phi_4096')},
        high_precision={k:mp.nstr(v,80) for k,v in values.items()},
        absolute_float_differences={k:mp.nstr(v,30) for k,v in differences.items()}))
report=dict(status='80_DIGIT_FLOATING_POINT_CHECK_NOT_INTERVAL_CERTIFICATE',precision_digits=80,
    derivative_source='F formula evaluated directly; mpmath differentiation; no SymPy or expression cache',
    mpmath_version=mp.__version__,elapsed_seconds=time.perf_counter()-started,rows=rows,
    maxima_absolute_float_differences={k:mp.nstr(max(mp.mpf(r['absolute_float_differences'][k]) for r in rows),30)
                                    for k in ('F','calT','T','D')})
(HERE.parent/'results/symbolic_audit/high_precision_witnesses.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))

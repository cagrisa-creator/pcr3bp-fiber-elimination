"""Optional resolution check for Figure 3 signs, not an interval certificate."""
import json
import time
from pathlib import Path
from audit_figures import analyze_grid

started=time.perf_counter()
rows=[]
for mu in (.2,.5,.8):
    base=analyze_grid(mu,200,80,192)
    refined=analyze_grid(mu,400,160,768)
    assert refined['minT']>0 and refined['minD']<0 and refined['max_dFdr']<0
    rows.append(dict(mu=mu,baseline=base,refined=refined))
report=dict(status='FLOATING_POINT_SIGN_REFINEMENT_NOT_CERTIFICATION',
    baseline=dict(n_theta=200,n_radial=80,n_phi=192),
    refined=dict(n_theta=400,n_radial=160,n_phi=768),rows=rows,
    elapsed_seconds=time.perf_counter()-started)
(Path(__file__).resolve().parent.parent/'results'/'figure_audit'/'lobe_refinement_report.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))

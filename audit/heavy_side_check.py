"""Reproduce the dedicated heavy-side sample quoted in Section 5.2.

The 600 x 240 x 512 grid matches the preserved audit_report_2026-07-08.json.
The collar is evaluated with its nonsingular F=0 formula. No interval claim.
"""
from pathlib import Path
import importlib.util,sys,json,time
here=Path(__file__).resolve().parent
root=here.parent
if (root/'frozen_inputs').exists():
 source=root/'frozen_inputs/code'
 out=root/'results/parameter_audit'
 helper=here/'audit_parameters.py'
else:
 source=root/'canonical/code'
 out=here/'results'
 helper=here/'audit_parameters.py'
sys.argv=[str(helper),'--source',str(source),'--output',str(out),'--suite','spot']
spec=importlib.util.spec_from_file_location('parameter_functions',helper)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
started=time.perf_counter()
r=audit.scan(.001,audit.cm.cL1_of(.001)+1e-4,nth=600,ns=240,nphi=512,npre=4000,nbisect=55,coeff_source='cache',boundary_complete=True,full=True)
r.update(label='dedicated_heavy_side_check',elapsed_seconds=time.perf_counter()-started,historical_reference=-.5783352446802084,interval_certificate=False)
audit.write_json('heavy_side_check.json',r)
audit.write_csv('heavy_side_witness.csv',[audit.flattened(r,'dedicated_heavy_side')])
print(json.dumps(audit.safe(r),indent=2))

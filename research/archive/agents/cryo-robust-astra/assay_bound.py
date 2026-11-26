"""Post hoc, sharp finite-sample constraints from Gao et al.'s triplicate assay.
No draws or distributional model. Independently parametrize all possible triplets."""
import json
import numpy as np
from investigate import coefficients,designs
mean,sd=1.47,.36
# Centered triplet lies in sum(x)=0 plane on circle of radius sqrt(2)*sample SD.
angles=np.arange(60000)*2*np.pi/60000
triplets=mean+(2*sd/np.sqrt(3))*np.cos(angles[:,None]+np.array([0,2*np.pi/3,4*np.pi/3]))
cap=float(coefficients(np.array(-150.),False)[0])
upper=mean-sd/np.sqrt(3)
rounding_upper=(mean+.005)-(sd-.005)/np.sqrt(3)
report=dict(assumed_sample_SD=sd,mean=mean,sharp_upper_on_smallest=upper,rounded_upper=rounding_upper,
            circle_max_of_minima=triplets.min(axis=1).max(),capacity_at_start=cap,tests={})
assert abs(upper-report['circle_max_of_minima'])<1e-10
for mode,sar in [('clamped',691.),('extrapolated',1013.)]:
    for design in ['caseII','margin','equal','original']:
        factor=60000*sar*designs[design][3]/cap
        threshold=50/factor
        counts=(triplets<threshold).sum(axis=1)
        report['tests'][mode+'_'+design]=dict(required_fe=threshold,min_failures_on_circle=int(counts.min()),
            max_failures_on_circle=int(counts.max()),max_possible_slowest_start_rate=upper*factor,
            rounding_conservative_slowest_start_rate=rounding_upper*factor)
print(json.dumps(report,indent=2))

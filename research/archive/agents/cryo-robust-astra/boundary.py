"""Post hoc bracket around the upper measured-scale scenario, constant clamped duty.
Also revisit the fixed-mean split at the next two spacings."""
import numpy as np
import pandas as pd
import investigate as m
rows=[]
for mm in [.5,.4,.35]:
    # Uniform muscle field, linear spacing in delivered iron. No distribution.
    p,_,_,labels,_=m.mesh(mm)
    for fe in [1.81,1.815,1.82,1.825]:
        _,faces,labels,iron,dt,duty=m.setup(mm,'equal','high',1,70)
        iron[labels==2]=fe
        out=np.zeros(15);end=np.zeros(len(iron))
        m.library(len(iron),len(faces),faces,labels,iron,p.spacing_m,dt,len(duty),duty,1,0.,
                  int(np.ceil(len(duty)/120)),out,end)
        rows.append(dict(mm=mm,fe=fe,min_rate=out[0],cold=out[8],hot=out[9],passes=out[0]>=50 and out[8]>=-35 and out[9]<=0))
    pd.DataFrame(rows).to_csv('outputs/loading_boundary.csv',index=False)
rows=[]
for mm in [.4,.35]:
    for scenario in ['deepLow','deepHigh','nominal']:
        rows.append(m.run(mm,'margin',scenario,'clamped'))
    pd.DataFrame(rows).to_csv('outputs/split_refinement.csv',index=False)

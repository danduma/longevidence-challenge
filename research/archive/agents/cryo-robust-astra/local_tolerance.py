"""Post hoc robustness radius around the high-loading, constant-duty counterexample.
Hold mean 1.83 fixed. Use the preregistered complementary spatial maps."""
import numpy as np
import pandas as pd
import investigate as m
rows=[]
for mm in [.5,.4,.35]:
    p,_,_,_,offsets=m.mesh(mm)
    for amp in [.005,.010,.015]:
        for sign in [1,-1]:
            p,faces,labels,iron,dt,duty=m.setup(mm,'equal','high',1,70)
            iron+=sign*amp*offsets/.36
            assert abs(iron[labels==2].mean()-1.83)<1e-12
            out=np.zeros(15);end=np.zeros(len(iron))
            m.library(len(iron),len(faces),faces,labels,iron,p.spacing_m,dt,len(duty),duty,1,0.,
                      int(np.ceil(len(duty)/120)),out,end)
            rows.append(dict(mm=mm,amplitude=amp,placement='deepLow' if sign==1 else 'deepHigh',
                             mean_fe=iron[labels==2].mean(),min_rate=out[0],window_rate=out[4],cold=out[8],hot=out[9],
                             passes=out[0]>=50 and out[8]>=-35 and out[9]<=0))
    pd.DataFrame(rows).to_csv('outputs/local_tolerance.csv',index=False)

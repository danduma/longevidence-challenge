"""Independent dense-grid method of lines, SSP RK2, adiabatic; data checks."""
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from investigate import coefficients,setup,run,Path,rat_heart_cryobag,designs

def independent(mm,design,scenario,mode):
    p,_,_,iron,dt,_=setup(mm,design,scenario,1,70)
    # Exact final time and continuous duty at physical time; differs from pinned stepping.
    steps=int(np.ceil(70/dt))
    dt=70/steps
    mask=p.active_mask
    spec=p.specimen_mask
    wall=p.labels==0
    fe=np.zeros(p.shape);fe[mask]=iron
    T=np.full(p.shape,-150.)
    links=[]
    for axis in range(3):
        a=[slice(None)]*3;b=a.copy()
        a[axis]=slice(None,-1);b[axis]=slice(1,None)
        a,b=tuple(a),tuple(b)
        links.append((a,b,mask[a]&mask[b]))
    def rhs(T,time):
        c,k=coefficients(T,False)
        c[wall],k[wall]=coefficients(T[wall],True)
        power=np.interp(time,[0,35,70],designs[design][3:])
        effective=np.minimum(T,0) if mode=='extrapolated' else np.clip(T,-80,0)
        q=1000*fe*(691-4.6*(effective+80))*power
        for a,b,link in links:
            flux=2/(1/k[a]+1/k[b])*(T[b]-T[a])*link/p.spacing_m**2
            q[a]+=flux
            q[b]-=flux
        return q/c
    minimum=np.inf; window=np.inf
    for step in range(steps):
        r=rhs(T,step*dt)
        minimum=min(minimum,r[spec].min()*60)
        window=min(window,np.min(60*r,where=spec&(T>=-100)&(T<=-40),initial=np.inf))
        T=T+dt*(r+rhs(T+dt*r,(step+1)*dt))*.5
    return np.array([minimum,window,T[spec].min(),T[spec].max()])

def raw():
    src=Path('sources/liter.xlsx')
    sheet=load_workbook(src,read_only=True,data_only=True)['Fig.7']
    rows=list(sheet.values)
    result=[]
    for offset,vol in [(1,1),(5,2)]:
        for k,probe in enumerate(['left','right','center'],1):
            a=np.array([(r[offset],r[offset+k]) for r in rows[5:] if r[offset] is not None and r[offset+k] is not None])
            time,T=a.T
            def cross(temp):
                idx=np.flatnonzero((T[:-1]<temp)&(T[1:]>=temp))[0]
                return time[idx]+(temp-T[idx])*np.diff(time)[idx]/np.diff(T)[idx]
            # Fig. 7b representative traces may stop before zero. Keep censoring.
            t1=cross(-100)
            inside=(T[:-1]>=-100)&(T[1:]<=0)&(time[:-1]>=t1)
            if (T>=0).any():inside &= time[:-1]<cross(0)
            rates=np.diff(T)[inside]/np.diff(time)[inside]*60
            result.append(dict(volume_L=vol,probe=probe,dt_s=np.median(np.diff(time)),n_intervals=len(rates),
                               raw_min_C_min=rates.min(),zero_intervals=(rates==0).sum(),
                               transit_to_zero=6000/(cross(0)-t1) if max(T)>=0 else np.nan,
                               rate_minus100_to_minus40=3600/(cross(-40)-t1),maximum_recorded_C=T.max()))
    pd.DataFrame(result).to_csv('outputs/raw_measurement_check.csv',index=False)
if __name__=='__main__':
    raw()
    records=[]
    for d,s,mode in [('margin','deepLow','extrapolated'),('equal','high','clamped'),('margin','high','clamped')]:
        mm=.75
        secondary=independent(mm,d,s,mode)
        primary=run(mm,d,s,mode)
        for field,val in zip(['min_rate','window_rate','cold','hot'],secondary):
            records.append(dict(design=d,scenario=s,sar=mode,mm=mm,metric=field,rk2=val,euler=primary[field],difference=val-primary[field]))
    f=pd.DataFrame(records)
    f.to_csv('outputs/independent.csv',index=False)
    print(f.to_string(index=False))
    assert f.difference.abs().max()<.3

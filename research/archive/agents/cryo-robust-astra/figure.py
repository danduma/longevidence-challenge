"""Plot observed grid refinements and the two rate criteria."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR',str(Path('.mplcache').resolve()))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
sweep=pd.read_csv('outputs/sweep.csv')
ref=pd.read_csv('outputs/refinement.csv')
fig,axs=plt.subplots(1,2,figsize=(10,4.3),layout='constrained')
for name,label,color in [('margin','Published margin schedule','#286585'),('equal','Constant, equal integrated duty','#b45323')]:
    a=pd.concat([sweep,ref]).query("sar=='clamped' and scenario=='high' and design==@name and dt_divisor==1").sort_values('mm',ascending=False)
    axs[0].plot(a.mm,a.min_rate,'o-',label=label,color=color,markersize=4)
    # Only computed scores, with no inferred uncertainty between individual cases.
    b=sweep.query("sar=='extrapolated' and mm==.5 and dt_divisor==1 and design==@name").set_index('scenario')
    for i,k in enumerate(['nominal','deepLow','deepHigh']):
        x=i+(.08 if name=='equal' else -.08)
        axs[1].plot([x,x],[b.loc[k,'min_rate'],b.loc[k,'window_rate']],color=color,lw=1)
        axs[1].scatter(x,b.loc[k,'min_rate'],color=color)
        axs[1].scatter(x,b.loc[k,'window_rate'],facecolor='white',edgecolor=color)
axs[0].set(xlabel='Voxel spacing (mm)',ylabel='Minimum warming rate (°C/min)',title='Clamped SAR · muscle at 1.83 mg/mL',xlim=(1.47,.23),ylim=(42,53))
axs[0].legend(loc='upper left',frameon=False,fontsize=8)
axs[1].set(xticks=[0,1,2],xticklabels=['Uniform','Deep half low','Deep half high'],ylabel='Minimum warming rate (°C/min)',title='Extrapolated SAR · same mean muscle iron',ylim=(35,82))
axs[1].text(.04,.97,'○  −100 to −40 °C     ●  whole run',transform=axs[1].transAxes,va='top',fontsize=9)
for a in axs:
    a.axhline(50,color='0.35',ls='--',lw=1)
    a.spines[['right','top']].set_visible(False)
fig.savefig('outputs/thermal_decisions.png',dpi=190)

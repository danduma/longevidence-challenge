from pathlib import Path
import pandas as pd,numpy as np
import os
os.environ['MPLCONFIGDIR']=str(Path(__file__).resolve().parents[1]/'.cache/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
a=pd.read_csv(R/'results/library2025_calibration_pairs.csv');a=a.sort_values('effective_factor')
labels=[s if len(s)<43 else s[:40]+'…' for s in a.species]
fig,ax=plt.subplots(figsize=(10,6.5));yy=np.arange(len(a))
ax.axvspan(1.7,2.8,color='#d1dce8',alpha=.8,label='Published CI for common factor (not prediction interval)')
ax.axvline(2.1292109150998813,color='#31495d',ls='--',label='Deposited factor 2.1292')
ax.scatter(a.effective_factor,yy,color='#0e6d78',s=45,zorder=3)
ax.set_yticks(yy,labels);ax.set_xlabel('Effective factor: authors’ temperature-standardized literature rate / raw coupled rate')
ax.set_title('2025 calibration overlaps: all original literature temperatures are 25°C')
ax.set_xlim(.8,5.2);ax.grid(axis='x',alpha=.2);ax.legend(loc='lower right',fontsize=8);fig.tight_layout()
for ext in ['pdf','png']:fig.savefig(R/'results'/('calibration2025.'+ext),dpi=180)
plt.close(fig)
r=pd.read_csv(R/'results/control_rectangles.csv');fig,ax=plt.subplots(figsize=(6.5,4))
for trait,color in [('S','#0e6d78'),('vC','#bd6538')]:
 z=np.sort(r[(r['mode']=='condition')&(r.trait==trait)].factor);ax.step(z,np.arange(1,len(z)+1)/len(z),where='post',label=f'{trait}, {len(z)} dependent rectangles',color=color)
ax.axvline(1.1,color='gray',ls=':',lw=1);ax.axvline(1.25,color='gray',ls=':',lw=1)
ax.set_xlabel('Factor discrepancy between within-source enzyme ratios');ax.set_ylabel('Empirical cumulative fraction');ax.set_ylim(0,1.02);ax.set_title('Historical control transport: source overlap is incomplete');ax.legend(fontsize=9);ax.grid(alpha=.2);fig.tight_layout()
for ext in ['pdf','png']:fig.savefig(R/'results'/('control_transport.'+ext),dpi=180)
plt.close(fig)

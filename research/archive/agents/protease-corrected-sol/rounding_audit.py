import numpy as np,pandas as pd
from reconstruct import O

def rround(x,digits):
 x=np.asarray(x,dtype=float);sgn=np.sign(x);z=np.abs(x);p=10.**digits;t=z*p;i=np.floor(t);down=i/p;up=np.ceil(t)/p
 return sgn*np.where((up-z<z-down)|((i%2==1)&(up-z==z-down)),up,down)
if __name__=='__main__':
 d=pd.read_csv(O/'DMS_script05_singles.tsv.gz',sep='\t');d['key']=d.TEVp+'_'+d.TEVs
 # Read all source rows once, source discrepancies are insufficient as reference.
 import openpyxl
 w=openpyxl.load_workbook('sources/41467_2025_60622_MOESM4_ESM.xlsx',read_only=True,data_only=True)
 a=pd.DataFrame(w['Fig. S8e'].values);a.columns=a.iloc[0];a=a.iloc[1:];a.columns=['key','target','category']
 d=d.merge(a,on='key');y=d.target.astype(float).to_numpy()
 n=d[[f'N{i}' for i in range(1,6)]].to_numpy();f=d[[f'F{i}' for i in range(1,6)]].to_numpy();q=f/(n+f)
 repl=d.source.str.startswith('repl').to_numpy();rounds=[('numpy',np.round),('R',rround)]
 for label,fun in rounds:
  t=q.copy();t[repl]=fun(t[repl],4)
  for scale in [1,60]:
   z=(scale*t[:,:-1]+scale*(t[:,1:]-t[:,:-1])*.5)
   z=z/(4*scale) if scale==60 else z
   for method in ['numpy','sequential','fsum']:
    if method=='numpy':v=z.sum(axis=1)
    elif method=='sequential':v=((z[:,0]+z[:,1])+z[:,2])+z[:,3]
    else:
     import math
     v=np.array([math.fsum(row) for row in z])
    if scale==1:v=v/4
    # original rows don't have minute scaling
    raw=(q[:,:-1]+(q[:,1:]-q[:,:-1])*.5)
    if method=='fsum':raw=np.array([math.fsum(row) for row in raw])/4
    else:raw=raw.sum(axis=1)/4
    v[~repl]=raw[~repl]
    res=fun(v,4)
    print(label,scale,method,int((abs(res-y)>1e-12).sum()),flush=True)

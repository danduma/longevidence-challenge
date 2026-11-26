"""Post-sweep refinement of the small positive high-loading, clamped-SAR margin."""
from investigate import run, coefficients, designs
import pandas as pd
import numpy as np

rows=[]
for mm in [.625,.4,.35,.3]:
    for d in ['equal','margin']:
        rows.append(run(mm,d,'high','clamped'))
    pd.DataFrame(rows).to_csv('outputs/refinement.csv',index=False)
    print(mm,[(r['min_rate'],r['cold'],r['hot']) for r in rows[-2:]],flush=True)
rows.append(run(.4,'equal','high','clamped',half=2))
pd.DataFrame(rows).to_csv('outputs/refinement.csv',index=False)
for mode,sar in [('clamped',691.),('extrapolated',1013.)]:
    c=coefficients(np.array(-150.),False)[0]
    for name in ['caseII','margin','equal','original']:
        print('necessary initial Fe',mode,name,50*c/(60000*sar*designs[name][3]))

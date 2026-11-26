"""Explore the unresolved tissue-placement comparison below 0.35 mm."""
from investigate import run,pd
rows=[]
for mm in [.3,.25]:
    for scenario in ['deepLow','deepHigh']:
        rows.append(run(mm,'margin',scenario,'clamped'))
        pd.DataFrame(rows).to_csv('outputs/split_fine.csv',index=False)
        print(mm,scenario,rows[-1]['cold'],flush=True)

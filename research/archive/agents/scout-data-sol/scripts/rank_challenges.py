"""Transparent subjective prioritization for bounded public-data computational work."""
import pandas as pd,json
rows=[
(1,'Origins of life',3.5,2.5,5,1.5),
(2,'Cryopreservation',3.5,3.5,5,2),
(3,'The reverse translatase',3,1.5,4.5,2),
(4,'Improve Rubisco',5,4.5,5,4.5),
(5,'The quadruplet cell',4,3.5,4,3.5),
(6,'Somatic limb regeneration',5,3.5,5,2),
(7,'Bacterial production of gene therapies',3,2.5,4.5,2.5),
(8,'Programmable proteases',5,4,5,4),
(9,'Cell-penetrating protein binders',3.5,3,5,3),
(10,'Protein amplification chain reaction',2,1,4.5,1.5),
(11,'5′ polymerases',2.5,3,4,3),
(12,'New nitrogenases',4,2.5,5,2)]
d=pd.DataFrame(rows,columns=['challenge_id','challenge_name','data_availability','computational_tractability','impact','validation_path'])
criteria=d.columns[2:]
weights={'primary':[.35,.30,.20,.15],'equal':[.25]*4,'data_first':[.50,.30,.10,.10],'validation_emphasis':[.20,.20,.15,.45]}
for n,w in weights.items():d[n+'_score']=d[criteria].to_numpy()@w
# Ties broken by validation path, then tractability, then stable challenge id.
for n in weights:
 order=d.sort_values([n+'_score','validation_path','computational_tractability','challenge_id'],ascending=[False,False,False,True]).index
 d[n+'_rank']=pd.Series(range(1,len(d)+1),index=order)
d=d.sort_values('primary_rank');d.to_csv('results/challenge_ranking.csv',index=False)
open('results/ranking_method.json','w').write(json.dumps({'scales':'1 (poor) to 5 (strong), half-points allowed; subjective evidence-informed estimates, not success probabilities','weights':weights,'tie_break':'validation path, computational tractability, challenge id','scope':'meaningful public-data computational subproblems, not likelihood of solving whole experimental challenge'},indent=2))
print(d[['primary_rank','challenge_name','primary_score','equal_rank','data_first_rank','validation_emphasis_rank']].to_string(index=False))

import unittest
import numpy as np
import pandas as pd
from reconstruct import measures, correct

class TransformTests(unittest.TestCase):
 def test_auc_uses_unrounded_fraction_and_counts(self):
  d=pd.DataFrame({**{f'N{i}':[2] for i in range(1,6)},**{f'F{i}':[1] for i in range(1,6)}})
  x=measures(d)
  self.assertEqual(x.AUC.iloc[0],.3333)
  self.assertEqual(x['min'].iloc[0],3)
 def test_matches_deposited_rounding_boundary(self):
  f=[0,161,990,1922,2338]
  d=pd.DataFrame({**{f'N{i+1}':[10000-v] for i,v in enumerate(f)},**{f'F{i+1}':[v] for i,v in enumerate(f)}})
  self.assertEqual(measures(d,fraction_round=4).AUC.iloc[0],.1061)
 def test_outlier_and_duplicate_are_distinct(self):
  rows=[]
  for enzyme,flips in [('outlier',[90,10,10]),('duplicate',[90,10]),('stable',[10,10,10])]:
   for j,f in enumerate(flips):
    rows.append(dict(TEVp=enzyme,TEVs='ENLYFQS',BC_TEVp=enzyme+str(j),BC_TEVs='a',n_mut_TEVs=0,**{f'F{i}':f for i in range(1,6)},**{f'N{i}':100-f for i in range(1,6)}))
  b=measures(pd.DataFrame(rows));res,ledger,bc,meta=correct(b,.15,'none','test',save=False)
  self.assertNotIn('duplicate',res.TEVp.tolist())
  self.assertAlmostEqual(res.set_index('TEVp').loc['outlier','AUC'],.1)
  self.assertEqual(meta['outlier_barcodes'],1)
  self.assertEqual(meta['diverging_pairs'],1)
if __name__=='__main__':unittest.main()

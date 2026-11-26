#!/usr/bin/env python3
"""Independent audit: exact fractions and cubic brackets, no numpy/scipy.

Numbers are independently transcribed from Patel 2021 Table 1 image, not read
from shared CSV or implementation. Checks are conditional algebra, not biology.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from pathlib import Path
import csv, hashlib, json, platform
BASE=Path(__file__).resolve().parent
SHARED=Path('/path/to/longevidence-checkout/research/millennium')
# template incoming kpol +/- KD +/- printed efficiency; C,G is lower-censored.
RAW='''A G 2.84 .08 23 2 .12
A A .032 .001 70 15 .00046
A U .77 .07 140 50 .0055
A C .009 .001 200 130 .000045
C G 10 - 1 - 10
C A .065 .005 160 40 .00041
C U .14 .02 510 160 .00027
C C .10 .01 3000 500 .000033
G G 3.56 .39 22 8 .16
G A .16 .01 46 13 .0035
G U .11 .01 452 170 .00024
G C 1.38 .07 46 9 .030
U G 4.88 .44 12 4 .41
U A .46 .03 20 6 .023
U U .25 .02 362 72 .00069
U C .07 .008 628 258 .00011'''
DATA={}
for line in RAW.splitlines():
 t,n,*v=line.split(); DATA[t,n]=[F(s) if s!='-' else None for s in v]
NT='ACGU'; COMP=dict(zip(NT,'UGCA'))
def dec(q,digits=45):
 with localcontext() as c:
  c.prec=digits
  return str(Decimal(q.numerator)/Decimal(q.denominator))
def matrix(eff):
 return [[F(0) if i==j else eff[COMP[i],j]/eff[COMP[i],i] for j in NT] for i in NT]
def principal(b): return [[b[NT.index(i)][NT.index(j)] for j in 'ACU'] for i in 'ACU']
def cw(b):
 """Quantized rational iteration proposes vector; exact residual certifies it.
 Every positive vector yields valid CW enclosure even if iteration fails.
 """
 scale=10**45; x=[F(1)]*len(b)
 for _ in range(1800):
  y=[sum(v*z for v,z in zip(row,x)) for row in b]
  m=max(y)
  x=[F(max(1,(v*scale/m).__floor__()),scale) for v in y]
 quot=[sum(v*z for v,z in zip(row,x))/xx for row,xx in zip(b,x)]
 lo,hi=min(quot),max(quot)
 assert hi-lo<F(1,10**35), (float(lo),float(hi))
 return {'rho_lower_exact':str(lo),'rho_upper_exact':str(hi),
         'probability_lower_exact':str(1/(1+hi)), 'probability_upper_exact':str(1/(1+lo)),
         'probability_decimal':dec(1/(1+lo)),
         'enclosure_width':dec(1/(1+lo)-1/(1+hi)),
         'positive_vector_exact':[str(v) for v in x],
         'normalized_mixture':[dec(v/sum(x),18) for v in x]}
def cubic(b):
 a,b1=b[0][1:]; c,d=b[1][0],b[1][2]; e,f=b[2][:2]
 s=a*c+b1*e+d*f; t=a*d*e+b1*c*f
 # p(r)=r^3-s*r-t has exactly one positive root, by Descartes.
 p=lambda r:r**3-s*r-t
 lo,hi=F(0),F(1)
 while p(hi)<=0: hi*=2
 for _ in range(180):
  mid=(lo+hi)/2
  if p(mid)<0: lo=mid
  else: hi=mid
 assert p(lo)<0<p(hi)
 return {'coefficient_s_exact':str(s),'coefficient_t_exact':str(t),
         'rho_lower_exact':str(lo),'rho_upper_exact':str(hi),
         'probability_lower_exact':str(1/(1+hi)),'probability_upper_exact':str(1/(1+lo)),
         'probability_decimal':dec(1/(1+lo)),
         'p_at_1_over_99_exact':str(p(F(1,99))),
         'reject_99_percent_exact':p(F(1,99))<0,
         'matrix_exact':[[str(v) for v in row] for row in b]}
def envelope(m):
 efflo={}; effhi={}
 for key,(k,ke,kd,kde,_) in DATA.items():
  if ke is None: continue
  efflo[key]=max(F(0),k-m*ke)/(kd+m*kde)
  # Only the three correct denominators matter; some mismatch KD boxes cross 0.
  if key[1]==COMP[key[0]]:
   assert kd-m*kde>0
   effhi[key]=(k+m*ke)/(kd-m*kde)
 b=[[F(0) if i==j else efflo[COMP[i],j]/effhi[COMP[i],i] for j in 'ACU'] for i in 'ACU']
 return cubic(b)
def main():
 csvpath=SHARED/'agents/scout-mechanisms-sol/data/thg1_table1_2021.csv'
 rows=list(csv.DictReader(csvpath.open()))
 for row in rows:
  expected=DATA[row['template'],row['ntp']]
  actual=[F(row[k]) if row[k] else None for k in ['kpol_per_min','kpol_sd','KD_uM','KD_sd','eff_uM_inv_min_inv']]
  assert actual==expected, (row,expected)
 e={key:v[-1] for key,v in DATA.items()}
 b=matrix(e); sub=principal(b)
 # The requested elementary certificate uses no spectral theorem.
 # U73 (correct A): A/C >=99*(.00011/.023); G73 (correct C): C/A >=99*(.0035/.030).
 cycle=(e['U','C']/e['U','A'])*(e['G','A']/e['G','C'])
 result={'python':platform.python_version(),'method':'exact Fraction Collatz-Wielandt residual + independent cubic rational bisection',
 'independent_table_matches_shared_csv':True,'shared_csv_sha256':hashlib.sha256(csvpath.read_bytes()).hexdigest(),
 'incoming_order':NT,'template_row_order':'UGCA','full_B_exact':[[str(v) for v in row] for row in b],
 'principal':cubic(sub),'principal_CW':cw(sub),
 'cycle':{'ratio_product_exact':str(cycle),'99_required_product_exact':str(99**2*cycle),
          '99_required_product_decimal':dec(99**2*cycle),'reject_99_exact':99**2*cycle>1,
          'A_over_C_required_exact':str(99*e['U','C']/e['U','A']),
          'C_over_A_required_exact':str(99*e['G','A']/e['G','C'])},
 'censoring_cases':{},'rectangles':{},
 'published_competition_comparison':{'source_doi':'10.1093/nar/gkq976',
 'A73_G_over_U_prediction_exact':str(e['A','G']/e['A','U']),
 'A73_G_over_U_prediction_decimal':dec(e['A','G']/e['A','U']),
 'reported_approximate_ratio':5,
 'interpretation':'Extrapolation diagnostic only: 2011 monophosphate/activation and high-NTP conditions differ from the dilute preactivated model.'}}
 # Independently derived rational cubic and CW brackets must intersect.
 pc=result['principal']; pw=result['principal_CW']
 assert max(F(pc['rho_lower_exact']),F(pw['rho_lower_exact'])) <= min(F(pc['rho_upper_exact']),F(pw['rho_upper_exact']))
 for v in [10,100,10000]:
  ee=e.copy(); ee['C','G']=F(v); result['censoring_cases'][str(v)]=cw(matrix(ee))
 for m in [0,1,2]: result['rectangles'][str(m)]=envelope(F(m))
 # Quantify how close the deterministic radius is to the 99% crossing.
 lo,hi=F(2),F(21,10)
 assert F(envelope(lo)['p_at_1_over_99_exact']) < 0 < F(envelope(hi)['p_at_1_over_99_exact'])
 # The leading sign of p(1/99) determines which side of the target the bound lies.
 for _ in range(60):
  mid=(lo+hi)/2
  if F(envelope(mid)['p_at_1_over_99_exact'])<0: lo=mid
  else: hi=mid
 result['rectangle_multiplier_99_crossing']={'lower':dec(lo,20),'upper':dec(hi,20)}
 shared=json.loads((SHARED/'polymerase/results/thg1_independent_verification.json').read_text())
 assert abs(float(result['principal']['probability_decimal'])-shared['principal_submatrix_probability_upper_bound_for_any_censored_correct_G_rate'])<1e-14
 assert abs(float(result['rectangles']['2']['probability_decimal'])-shared['favorable_parameter_rectangle_sensitivity'][2]['upper_bound_uniform_correct_probability'])<1e-14
 out=BASE/'outputs'; out.mkdir(exist_ok=True)
 (out/'independent_check.json').write_text(json.dumps(result,indent=2)+'\n')
 print('Independent table transcription: all 16 rows match.')
 print('Principal bound:',result['principal']['probability_decimal'])
 for m,r in result['rectangles'].items(): print('Rectangle',m,r['probability_decimal'],'exactly rejects 99%:',r['reject_99_percent_exact'])
 print('99% reciprocal ratio product:',result['cycle']['99_required_product_exact'])
 print('99% rectangle crossing multiplier:',result['rectangle_multiplier_99_crossing'])
 print('All exact and comparison assertions passed.')
if __name__=='__main__': main()

"""Conditional mixed-substrate competition analysis, not a validated elongation model.

Data: Patel et al. RNA (2021), Table 1, DOI 10.1261/rna.078686.121.
Assume mutually exclusive, dilute NTP capture: hazard E[template,ntp]*[NTP].
The table concerns single additions to preactivated tRNA, not repeated generic
elongation. Published rounded efficiencies, including a censored lower bound,
are preserved. No claim that concentrations can rescue processivity is made.
Python standard library only; run from the job root.
"""
import csv
import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NTS = 'ACGU'
COMP = dict(zip('ACGU', 'UGCA'))
ROWS = list(csv.DictReader((ROOT/'data/thg1_table1_2021.csv').open()))


def efficiencies():
    return {t: {r['ntp']: float(r['eff_uM_inv_min_inv'])
                for r in ROWS if r['template'] == t} for t in NTS}


def probabilities(e, x):
    return {t: e[t][COMP[t]]*x[COMP[t]] /
            sum(e[t][n]*x[n] for n in NTS) for t in NTS}


def perron(b):
    n=len(b)
    # I+B has the same Perron vector and removes period-2 cycling.
    x = [1/n]*n
    for it in range(20000):
        y = [x[i]+sum(b[i][j]*x[j] for j in range(n)) for i in range(n)]
        scale = sum(y)
        y = [v/scale for v in y]
        x = y
        rr = [sum(b[i][j]*x[j] for j in range(n))/x[i] for i in range(n)]
        if max(rr)-min(rr) < 1e-11:
            break
    ratios = [sum(b[i][j]*x[j] for j in range(n))/x[i] for i in range(n)]
    # Collatz-Wielandt lower and upper bounds certify convergence.
    assert max(ratios)-min(ratios) < 1e-9
    return sum(ratios)/n,x,ratios,it+1


def minimax(e):
    # Rows are indexed by the CORRECT incoming nucleotide, not the template.
    b = [[0.0 if i == j else e[COMP[i]][j]/e[COMP[i]][i]
          for j in NTS] for i in NTS]
    rho,x,ratios,iterations=perron(b)
    pools = dict(zip(NTS,x))
    probs = probabilities(e,pools)
    assert max(probs.values())-min(probs.values()) < 1e-9
    return {'rho':rho,'rho_bounds':[min(ratios),max(ratios)],
            'optimal_worst_correct_probability':1/(1+rho),
            'pool_fractions':pools,'correct_probabilities':probs,'iterations':iterations}


def lognormal_around_mean_sd(mean, sd, rng):
    sig2 = math.log1p((sd/mean)**2)
    return rng.lognormvariate(math.log(mean)-sig2/2,math.sqrt(sig2))


def main():
    e = efficiencies()
    equal = probabilities(e,{n:1 for n in NTS})
    fits=[]
    for val in [10,100,1e4,1e8]:
        ee = {t:dict(row) for t,row in e.items()}
        ee['C']['G']=val
        fits.append({'assumed_C_G_efficiency':val,**minimax(ee)})
    # Stronger analytic bound uses just the A<->C off-diagonal 2-cycle,
    # in incoming-nucleotide indexing (templates U and G).
    cycle = math.sqrt((e['U']['C']/e['U']['A'])*(e['G']['A']/e['G']['C']))
    hard_upper = 1/(1+cycle)
    # A,C,U principal submatrix removes the censored template-C row entirely.
    # A nonnegative matrix has spectral radius >= that of a principal submatrix.
    subset='ACU'
    b3=[[0.0 if i==j else e[COMP[i]][j]/e[COMP[i]][i]
         for j in subset] for i in subset]
    rho3,_,ratios3,_=perron(b3)
    # One pool cannot simultaneously give 99% correct at these two templates.
    odds99=0.99/0.01
    requirements99={'A_over_C_min_for_template_U':odds99*e['U']['C']/e['U']['A'],
                    'C_over_A_min_for_template_G':odds99*e['G']['A']/e['G']['C']}
    # Independent lognormal parameter perturbation: sensitivity, not a CI.
    # Covariances and raw replicate data are unavailable; censored row fixed.
    rng=random.Random(20261126)
    samples=[]
    for _ in range(2000):
        ee={t:{} for t in NTS}
        for r in ROWS:
            if r['censoring']!='none':v=10.0
            else:
                k=lognormal_around_mean_sd(float(r['kpol_per_min']),float(r['kpol_sd']),rng)
                kd=lognormal_around_mean_sd(float(r['KD_uM']),float(r['KD_sd']),rng)
                v=k/kd
            ee[r['template']][r['ntp']]=v
        samples.append(minimax(ee)['optimal_worst_correct_probability'])
    samples.sort()
    out={'source_doi':'10.1261/rna.078686.121','equal_pool_correct_probabilities':equal,
         'equal_pool_template_C_probability_is_lower_bound':True,
         'censored_row_sensitivity':fits,'two_cycle_rho_lower_bound':cycle,
         'two_cycle_probability_upper_bound':hard_upper,
         'censoring_independent_three_ntp_upper_bound':{
             'incoming_ntps':subset,'matrix':b3,'rho':rho3,
             'rho_bounds':[min(ratios3),max(ratios3)],
             'optimal_worst_correct_probability_upper_bound':1/(1+rho3)},
         'requirements_for_99_percent_pairwise_correctness':requirements99,
         'product_of_incompatible_ratio_requirements':math.prod(requirements99.values()),
         'parameter_perturbation_not_confidence_interval':{
             'n':len(samples),'seed':20261126,'p025':samples[49],
             'median':samples[999],'p975':samples[1949]},
         'warning':'Conditional capture model on tRNA single-addition data; no processivity or biological validation.'}
    (ROOT/'outputs/thg1_competition.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()

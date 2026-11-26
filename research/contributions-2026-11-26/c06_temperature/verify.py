#!/usr/bin/env python3
"""C06 verifier: re-extract pinned primary inputs, check all ledgers and invariants."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
import csv,hashlib,json,math,tempfile
import audit

def read(p):
    with p.open() as f:return list(csv.DictReader(f))

with tempfile.TemporaryDirectory(prefix='c06_verify_') as d:
    out=Path(d);summary=audit.run(out)
    checked=[]
    for p in sorted(out.iterdir()):
        committed=ROOT/'results'/p.name
        assert p.read_bytes()==committed.read_bytes(), 'Re-extracted output differs: '+p.name
        checked.append(p.name)
    points=read(out/'primary_temperature_points.csv')
    pairs=read(out/'measured_temperature_pairs.csv')
    # Independent arithmetic: compare log increments and telescoping ratios.
    for species in sorted({r['species'] for r in points}):
        pp={float(r['temp_C']):float(r['kcat_s_minus1']) for r in points if r['species']==species}
        gg={(float(r['low_C']),float(r['high_C'])):r for r in pairs if r['species']==species}
        assert set(pp)=={15,25,35} and len(gg)==3
        assert math.isclose(float(gg[15,25]['ratio'])*float(gg[25,35]['ratio']),float(gg[15,35]['ratio']),rel_tol=1e-14)
        for (a,b),r in gg.items():
            logq=(math.log(pp[b])-math.log(pp[a]))*10/(b-a)
            assert math.isclose(math.log(float(r['effective_Q10'])),logq,abs_tol=1e-14)
    # Fixed primary-table sentinels guard accidental column/temperature shifts.
    expected={('Avena sativa',15):1.3,('Avena sativa',35):3.9,('Zea mays',15):2.4,('Zea mays',25):4.1,('Zea mays',35):5.7,('Cucurbita maxima',15):.5}
    index={(r['species'],float(r['temp_C'])):float(r['kcat_s_minus1']) for r in points}
    for k,v in expected.items():assert index[k]==v,(k,index[k],v)
    candidates=read(out/'workbook_candidate_pairs.csv')
    eligible=[r for r in candidates if r['same_origin']=='True' and r['low_vC'] and r['high_vC']]
    assert len(eligible)==1 and eligible[0]['species']=='Fragilariopsis cylindrus'
    assert 'in_vivo_estimate' in eligible[0]['exclusion']
    assert summary['overlap_enzymes_with_verified_series']==0
    # Analytic rank proof: [1,-.5] rows share null [.5,1]; [0,1] breaks it.
    assert .5-.5*1==0 and 0*.5+1*1!=0
    assert summary['strict_augmented_rank']==1 and summary['conditional_universal_Q10_rank']==2
    contribution_path=ROOT/'contribution.json'
    if contribution_path.exists():
        contribution=json.loads(contribution_path.read_text())
        assert contribution['id']=='C06'
        assert contribution['status']=='agent-produced; awaiting coordinator verification'
        for claim in contribution['numerical_claims']:
            value=summary
            for key in claim['summary_path'].split('.'):
                value=value[key]
            assert value==claim['value'],claim
    record={'contribution_id':'C06','status':'PASS','all_input_sha256_checked':summary['input_files_checked'],'byte_identical_reextracted_outputs':checked,'primary_table_sentinels':len(expected),'species_series_arithmetic_verified':len(index)//3,'no_network':True,'no_LLM_calls':True,'human_expert_review':False}
    (ROOT/'results/verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))

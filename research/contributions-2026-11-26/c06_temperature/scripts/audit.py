#!/usr/bin/env python3
"""C06: offline extraction from immutable primary/local inputs; no network or LLM."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'vendor'))
import argparse, csv, hashlib, itertools, json, math, platform, re
from html.parser import HTMLParser
import numpy as np
import pandas as pd
from pypdf import PdfReader

class Tables(HTMLParser):
    def __init__(self):
        super().__init__(); self.tables=[]; self.table=None; self.row=None; self.cell=None
    def handle_starttag(self, tag, attrs):
        if tag=='table': self.table=[]
        elif self.table is not None and tag=='tr': self.row=[]
        elif self.row is not None and tag in ('td','th'): self.cell=[]
    def handle_data(self, data):
        if self.cell is not None: self.cell.append(data)
    def handle_endtag(self, tag):
        if tag in ('td','th') and self.cell is not None:
            self.row.append(' '.join(''.join(self.cell).split())); self.cell=None
        elif tag=='tr' and self.row is not None:
            self.table.append(self.row); self.row=None
        elif tag=='table' and self.table is not None:
            self.tables.append(self.table); self.table=None

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def records(path):
    with path.open() as f: return list(csv.DictReader(f))
def number(x):
    try: return float(x) if math.isfinite(float(x)) else None
    except (ValueError,TypeError): return None
def write(out, name, rows, columns=None):
    columns = columns or list(rows[0])
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
def norm(x):
    if pd.isna(x): return ''
    s=str(x).strip()
    return str(float(s)) if number(s) is not None else s

def run(out):
    out.mkdir(parents=True,exist_ok=True); inp=ROOT/'inputs'
    manifest=json.loads((inp/'manifest.json').read_text())
    for item in manifest:
        assert sha(ROOT/item['path'])==item['sha256'], 'Input changed: '+item['path']
    access=[]
    for line in (inp/'retrievals.jsonl').read_text().splitlines():
        r=json.loads(line);p=inp/r['file']
        kind='transport_error'
        if 'error' not in r:
            data=p.read_bytes()
            assert sha(p)==r['sha256']
            kind='verified_PDF_payload' if data.startswith(b'%PDF') else 'verified_DOCX_payload' if data.startswith(b'PK') else 'XML_error_payload' if b'errorBean' in data else 'HTML_or_other_nondata_payload'
        access.append(dict(file=r['file'],url=r['url'],retrieved_utc=r.get('retrieved_utc',''),http_status=r.get('http_status',''),content_type=r.get('content_type',''),outcome=kind,error=r.get('error',''),sha256=r.get('sha256','')))
    write(out,'source_access_ledger.csv',access)
    # Reconcile every supplied row/cell to the actual original, including exclusion flags.
    w=records(inp/'workbook_rows_unfiltered.csv')
    x=pd.read_excel(inp/'rubisco2019_original.xlsx',sheet_name='kinetic_data')
    assert len(w)==len(x)
    for i,r in enumerate(w):
        for c in x.columns:
            assert norm(r[c])==norm(x.iloc[i][c]), (i,c,r[c],x.iloc[i][c])
    rowledger=[]
    for i,r in enumerate(w):
        rr=dict(r); rr['excel_row']=i+2
        rr['temperature_status']='missing' if number(r['temp_C']) is None else 'reported'
        rr['turnover_status']='missing_or_nonpositive' if not number(r['vC']) or number(r['vC'])<=0 else 'positive_compiled'
        rr['primary_numeric_status']='compilation_only_not_independently_verified_here'
        rowledger.append(rr)
    write(out,'workbook_row_ledger.csv',rowledger)
    # Coarse species grouping deliberately broadens the candidate search. Never use it
    # as proof of exact identity. Every cross-temperature pair gets an exclusion reason.
    cand=[]
    for a,b in itertools.combinations(rowledger,2):
        if a['species']!=b['species']:continue
        ta,tb=number(a['temp_C']),number(b['temp_C'])
        if ta is None or tb is None or ta==tb:continue
        if ta>tb:a,b=b,a;ta,tb=tb,ta
        reasons=[]
        if a['short_ref']!=b['short_ref']:reasons.append('different_assay_origin')
        if number(a['vC']) is None or number(b['vC']) is None:reasons.append('missing_turnover')
        if any('in-vivo' in r['note'].lower() for r in (a,b)):reasons.append('in_vivo_estimate_not_in_vitro_turnover')
        for c in ['isoform','variant','heterologous_expression','pH']:
            if norm(a[c])!=norm(b[c]):reasons.append(c+'_mismatch')
            elif not a[c]:reasons.append(c+'_unknown')
        if a['mutant']!=b['mutant']:reasons.append('mutant_flag_mismatch')
        cand.append(dict(low_row=a['excel_row'],high_row=b['excel_row'],species=a['species'],low_source=a['short_ref'],high_source=b['short_ref'],low_C=ta,high_C=tb,low_vC=a['vC'],high_vC=b['vC'],same_origin=a['short_ref']==b['short_ref'],exclusion=';'.join(reasons),eligible=False,ratio='',effective_Q10=''))
    assert all(c['exclusion'] for c in cand)
    write(out,'workbook_candidate_pairs.csv',cand)
    inventory=[]
    for ref in sorted({r['short_ref'] for r in w}):
        group=[r for r in w if r['short_ref']==ref]
        inventory.append(dict(source=ref,rows=len(group),temperatures=';'.join(map(str,sorted({number(r['temp_C']) for r in group if number(r['temp_C']) is not None}))),missing_temperature=sum(number(r['temp_C']) is None for r in group),positive_turnover=sum(number(r['vC']) is not None and number(r['vC'])>0 for r in group),primary_flags=';'.join(sorted({r['primary'] for r in group})),public_expansion='extracted_original_TableI_and_TableS1' if ref=='Hermida-Carrera 2016' else 'supplement_retrieval_failed_no_numeric_series' if ref=='Orr 2016' else 'not_expanded_beyond_supplied_workbook'))
    write(out,'source_inventory.csv',inventory)
    # Primary 2016 crop Table I and source-author supplemental Table S1.
    p=Tables();p.feed((inp/'hermida2016.html').read_text())
    table=p.tables[0]
    assert table[0][0]=='Species' and 'kcat' in ''.join(table[0])
    source25={r[0]:r for r in table if len(r)==8 and '±' in r[3] and not ('average' in r[0])}
    pdf=PdfReader(inp/'hermida_supp_ljmu_correct.pdf')
    text='\n'.join(pdf.pages[i].extract_text(extraction_mode='layout') for i in [0,1])
    (out/'hermida_supp_extracted.txt').write_text(text)
    supp={}
    for line in text.splitlines():
        nums=re.findall(r'(\d+(?:\.\d+)?)\s*±\s*(\d+(?:\.\d+)?)',line)
        if len(nums)==12 and 'average' not in line:
            species=line[:re.search(r'\d',line).start()].strip()
            supp[species]=nums
    assert len(supp)==20 and len(source25)==20 and set(supp)==set(source25)
    cultivar={'Avena sativa':'Forridena','Beta vulgaris':'Detroit','Brassica oleracea':'var italica Calabres','Capsicum annuum':'Picante','Coffea arabica':'Catuaí Vermelho IAC 44','Cucurbita maxima':'Totanera','Glycine max':'Ransom','Hordeum vulgare':'ssp vulgare Morex','Ipomoea batatas':'Rosa de Málaga','Lactuca sativa':'Cogollo de Tudela','Manihot esculenta':'not specified by source','Medicago sativa':'Aragón','Oryza sativa':'Bomba','Phaseolus vulgaris':'Contender','Saccharum × officinarum':'officinarum × spontaneum hybrid; cultivar unstated','Solanum lycopersicum':'Roma VF','Solanum tuberosum':'Erlanger','Spinacia oleracea':'Butterfly','Triticum aestivum':'Cajeme','Zea mays':'Carella'}
    points=[]
    for species in sorted(supp):
        v25=re.fullmatch(r'(\d+(?:\.\d+)?)\s*±\s*(\d+(?:\.\d+)?)',source25[species][3]);assert v25
        for temp,(value,se),loc in [(15,supp[species][4],'Table S1, kcatc 15°C'),(25,v25.groups(),'Table I, kcatc'),(35,supp[species][5],'Table S1, kcatc 35°C')]:
            points.append(dict(study='Hermida-Carrera 2016',doi='10.1104/pp.16.01846',species=species,cultivar=cultivar[species],identity_scope='same source cultivar/native leaf holoenzyme population; molecular isoform mixture unresolved',temp_C=temp,kcat_s_minus1=float(value),reported_SE=float(se),error_kind='SE; rounded zero retained',n='3–6 turnover biological replicates in methods; table general caption 3–9; exact n/covariance unavailable',pH=8,assay='14CO2 acid-stable fixation; Michaelis–Menten Vmax',normalization='14C-CABP binding active sites',gas='N2 for Kc and turnover; 21% O2 also used for Kcair',origin='present-source measurements',source_location=loc))
    write(out,'primary_temperature_points.csv',points)
    pairs=[]
    for species in sorted(supp):
        g=[p for p in points if p['species']==species]
        for a,b in itertools.combinations(g,2):
            ratio=b['kcat_s_minus1']/a['kcat_s_minus1'];q=ratio**(10/(b['temp_C']-a['temp_C']))
            pairs.append(dict(study=a['study'],species=species,cultivar=a['cultivar'],low_C=a['temp_C'],high_C=b['temp_C'],low_kcat=a['kcat_s_minus1'],high_kcat=b['kcat_s_minus1'],ratio=ratio,effective_Q10=q,exact_25_to_30=False,within_source_eligible=True,calibration_bridge_eligible=False,calibration_exclusion='no calibration enzyme identity match; interval differs; cross-assay transport unsupported',independence='three correlated contrasts per one cultivar series; not independent replicates',uncertainty='not propagated: covariance and exact n unavailable'))
    write(out,'measured_temperature_pairs.csv',pairs)
    library=records(inp/'library2025_original.csv')
    overlaps=[r for r in library if number(r['from_lab'])==1 and number(r['from_literature'])==1]
    assert {r['identifier'] for r in overlaps}=={r['identifier'] for r in records(inp/'calibration_pairs.csv')}
    assert all(r['prok_euk']=='prok' for r in overlaps)
    write(out,'original_calibration_overlaps.csv',overlaps)
    coverage=[];matrix=[]
    for r in overlaps:
        ts=[float(t) for t in r['Tdeg_lit'].split('/')];assert len(set(ts))==1
        matrix.append([1,-(30-ts[0])/10])
        nominal=[z for z in cand if z['species'] in r['species'].split('/') and z['low_vC'] and z['high_vC']]
        species_matches=set(supp).intersection(r['species'].split('/'))
        assert not species_matches, 'New species overlap requires exact construct audit'
        coverage.append(dict(accession=r['identifier'],author_species=r['species'],form=r['form'],literature_temperatures=r['Tdeg_lit'],coupled_C=30,verified_same_enzyme_series=0,plant_series_identity_matches=len(species_matches),nominal_compiled_cross_temperature_candidates=len(nominal),candidate_rows=';'.join(f"{z['low_row']}:{z['high_row']}" for z in nominal),candidate_disposition='all differ in source or endpoint/identity support; see workbook pair ledger',assay_mismatch='2025 NADH coupled assay vs primary 14CO2/CABP; no same-condition comparison',gas_mismatch='2025 4% CO2/0.4% O2 vs N2 fitted Vmax in plant series',pH_match='not established for exact 2025 constructs',normalization_match='CABP concept shared; quantification implementations differ',temperature_support='plant measurements 15/25/35; no measured 30°C point',supported_bridge=False))
    write(out,'calibration_overlap_coverage.csv',coverage)
    m=np.array(matrix);base_rank=int(np.linalg.matrix_rank(m))
    # Actual plant contrast adds information only if its beta is stipulated to equal
    # the microbial beta. With distinct beta_plant the microbial null vector survives.
    augmented=np.vstack([m]+[[0,(p['high_C']-p['low_C'])/10] for p in pairs])
    separate=np.zeros((len(matrix)+len(pairs),3));separate[:len(matrix),:2]=m
    separate[len(matrix):,2]=[(p['high_C']-p['low_C'])/10 for p in pairs]
    null=np.array([.5,1.,0.]);assert np.allclose(separate@null,0)
    rankrows=[dict(scenario='original_overlaps_prior_result',parameters=2,rank=base_rank,exact_calibration_bridges=0,interpretation='known baseline; not C06 novelty'),dict(scenario='append_verified_exact_calibration_bridges',parameters=2,rank=base_rank,exact_calibration_bridges=0,interpretation='no eligible new row'),dict(scenario='force_universal_plant_and_microbe_Q10',parameters=2,rank=int(np.linalg.matrix_rank(augmented)),exact_calibration_bridges=0,interpretation='conditional algebraic rank gain requires unverified cross-enzyme transfer'),dict(scenario='allow_separate_plant_Q10',parameters=3,rank=int(np.linalg.matrix_rank(separate)),exact_calibration_bridges=0,interpretation='microbial alpha and beta remain confounded; explicit null vector [0.5,1,0]')]
    write(out,'rank_tests.csv',rankrows)
    summary=dict(contribution_id='C06',workbook_rows=len(w),workbook_sources=len(inventory),workbook_missing_temperature=sum(number(r['temp_C']) is None for r in w),workbook_non25_reported_temperature=sum(number(r['temp_C']) is not None and number(r['temp_C'])!=25 for r in w),workbook_cross_temperature_candidates=len(cand),workbook_same_origin_candidates=sum(c['same_origin'] for c in cand),workbook_positive_turnover_same_origin_candidates=sum(c['same_origin'] and bool(c['low_vC']) and bool(c['high_vC']) for c in cand),workbook_eligible_in_vitro_pairs=0,primary_series=20,primary_points=len(points),primary_correlated_pairs=len(pairs),primary_source_studies=1,exact_25_to_30_pairs=0,overlap_enzymes=len(overlaps),overlap_enzymes_with_verified_series=0,strict_augmented_rank=base_rank,conditional_universal_Q10_rank=int(np.linalg.matrix_rank(augmented)),effective_Q10_by_interval={f'{low}-{high}':dict(n_series=20,minimum=min(p['effective_Q10'] for p in pairs if p['low_C']==low and p['high_C']==high),median=float(np.median([p['effective_Q10'] for p in pairs if p['low_C']==low and p['high_C']==high])),maximum=max(p['effective_Q10'] for p in pairs if p['low_C']==low and p['high_C']==high)) for low,high in [(15,25),(25,35),(15,35)]},input_files_checked=len(manifest))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'environment.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,pypdf='7.1.1',protocol_sha256=sha(ROOT/'protocol.md')),indent=2)+'\n')
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'results');args=ap.parse_args()
    print(json.dumps(run(args.out),indent=2))

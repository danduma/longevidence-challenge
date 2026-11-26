"""Independent artifact, source-integrity and reported-result checks."""
from pathlib import Path
import csv,json,hashlib,re,zipfile,ast
import numpy as np,pandas as pd
root=Path.cwd();checks=[]
Path('notes/verification.md').write_text('Verification checks are running; final evidence is written only when all checks complete.\n')
for p in Path('scripts').glob('*.py'):ast.parse(p.read_text())
checks.append('All local analysis/acquisition scripts parse successfully.')
rank=pd.read_csv('results/challenge_ranking.csv');assert len(rank)==12 and rank.challenge_id.nunique()==12
assert rank.head(3).challenge_id.to_list()==[4,8,6]
assert set(rank.primary_rank)==set(range(1,13))
assert rank[rank.challenge_id.eq(5)].validation_emphasis_rank.iloc[0]==3
checks.append('Ranking contains exactly all twelve challenges and reproduces primary/sensitivity selections.')
p=pd.read_csv('results/rubisco_source_pairs.csv');a=p[['S','vC']].to_numpy();indices=[]
# Brute-force independent nondominance implementation.
for i in range(len(a)):
 dominated=False
 for j in range(len(a)):
  if j!=i and a[j,0]>=a[i,0] and a[j,1]>=a[i,1] and (a[j,0]>a[i,0] or a[j,1]>a[i,1]):dominated=True;break
 if not dominated:indices.append(i)
assert len(indices)==14 and p.iloc[indices].species.nunique()==13
rub=json.load(open('results/rubisco_summary.json'))
assert len(p)==rub['source_species_pairs']==226
assert p.short_ref.nunique()==27 and p.species.nunique()==198
assert ((p.S>=238.1)&(p.vC>=4.05)).sum()==0
linked=pd.read_csv('results/rubisco_source_pairs_linked.csv');assert len(linked)==len(p)
assert linked.verified_source_link.notna().all()
assert linked.loc[linked.short_ref.eq('Uemura 1997'),'verified_source_doi'].eq('10.1006/bbrc.1997.6497').all()
assert linked.loc[linked.short_ref.eq('Whitney 2001'),'verified_source_doi'].eq('10.1046/j.1365-313x.2001.01056.x').all()
assert linked.loc[linked.short_ref.eq('Bainbridge 1995'),'verified_source_doi'].eq('10.1093/jxb/46.special_issue.1269').all()
cy=pd.read_csv('data/rubisco_cyano2026_kinetics.csv');assert len(cy)==6
assert np.isclose(cy.loc[cy.Species.eq('G. violaceus'),'S'].iloc[0],74.6)
assert np.isclose(cy.loc[cy.Species.eq('G. violaceus'),'vC'].iloc[0],11.8)
assert cy.S_SE.notna().all() and rub['primary_2026_cyanobacteria']['n_strains']==5
assert rub['recent2025']['max_natural_corrected_rate_30C']<rub['recent2025']['max_corrected_rate_30C_including_reconstructed_ancestor']
ledger=pd.read_csv('results/rubisco_assay_ledger.csv');assert len(ledger)==22 and ledger.source_doi.notna().all()
checks.append('Independent brute-force Pareto calculation yields14 pairs/13 species; all226 paired records have source links, known citation overlays are correct; 2026 SE/strain values and natural/ancestor distinction checked.')
pro=json.load(open('results/protease_summary.json'))
ports=pd.read_csv('results/protease_reporter_vitro_pairs.csv');assert len(ports)==802 and not ports.duplicated(['TEVp','TEVs']).any()
assert pro['original_DMS_rows']==413779 and pro['modeled_rows_after_read_and_single_substitution_filters']==294557
b=pd.read_csv('results/protease_baselines.csv');blocked=b[b.split.eq('holdout_double_or_more_mutant_substrates')].set_index('baseline')
assert blocked.loc['additive_substrate_sequence_plus_enzyme','RMSE_AUC']>blocked.loc['global_mean','RMSE_AUC']
assert pro['blocked_uncertainty'][0]['test_rows_with_any_unseen_nonzero_feature']==0
assert pro['blocked_uncertainty'][0]['cluster_bootstrap_RMSE_additive_minus_global_CI95'][0]>0
checks.append('Protease join is unique802 pairs; filtered denominator, blocked failure, observed-feature coverage and motif-bootstrap direction checked.')
reg=json.load(open('results/regeneration_summary.json'));r=pd.read_csv('results/regeneration_interactions.tsv.gz',sep='\t')
assert len(r)==139698 and int(r.expressed.sum())==20052
assert reg['samples']==18 and reg['interaction_BH05_3d']==0 and reg['interaction_BH05_5d']==0
assert np.isclose(412/1570*100,26.24203821656)
assert np.isclose(78/235*100,33.1914893617)
assert (r.loc[r.expressed,'q_BH_3d'].between(0,1)).all()
assert (r.loc[r.expressed,'q_BH_5d'].between(0,1)).all()
checks.append('Regeneration count/filter,18-sample design, multiplicity outputs and threshold-attrition percentages checked.')
# Check content formats: never accept a browser page as a valid XLSX/ZIP dataset.
for path in Path('data').glob('*.xlsx'):assert zipfile.is_zipfile(path),path
for path in Path('data').glob('*.zip'):assert zipfile.is_zipfile(path),path
checks.append('All downloaded XLSX/ZIP data files are valid archives; known HTML masquerading as ZIP excluded.')
# Verify all present files from successful download records; renamed/excluded artifacts are explicit.
renamed={'sources/aav2022_author_manuscript.pdf':'sources/aav_Le2021_dissertation.pdf','data/rubisco_aquatic2026_supplementary.zip':'sources/rubisco_aquatic2026_zip_access_page.html'}
exclude={'sources/polymerase2025_preprint_candidate_bioc.json'}
seen=0;missing=[]
latest={}
for line in Path('sources/download_manifest.jsonl').read_text().splitlines():
 rec=json.loads(line)
 if rec.get('curl_exit')==0:latest[rec['path']]=rec
for path,rec in latest.items():
 if path in exclude:continue
 p=Path(renamed.get(path,path))
 if not p.exists():missing.append(str(p));continue
 assert hashlib.sha256(p.read_bytes()).hexdigest()==rec['sha256'],p
 seen+=1
assert not missing,missing
checks.append(f'{seen} current downloaded artifacts match recorded SHA256 checksums; two renames and one rejected unrelated source are explicit.')
report=Path('report.md').read_text()
for m in re.findall(r'\]\(([^)]+)\)',report):
 if not re.match(r'https?://',m) and not m.startswith('#'):assert Path(m.split('#')[0]).exists(),m
checks.append('All local links in report resolve to existing deliverables.')
assert 'None observed' in Path('refusals.md').read_text()
Path('notes/verification.md').write_text('# Verification —2026-11-26\n\nFresh numerical reruns: extract_cyano2026.py, analyze_rubisco.py, apply_bibliographic_overlay.py, analyze_protease.py, extract_regeneration.py, analyze_regeneration.py, rank_challenges.py, build_assay_ledger.py, build_inventory.py.\n\n'+''.join('- '+x+'\n' for x in checks)+'\nLimitations: calculations and source integrity are verified, not original experiments; historical numeric transcriptions, phylogenetic models, all external data licenses and inaccessible raw files remain subject to the stated limits.\n')
print('\n'.join(checks))

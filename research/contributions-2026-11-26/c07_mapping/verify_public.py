"""Check public hash-pinned repair and derived exposure arithmetic, without source papers."""
from pathlib import Path
import json,csv,hashlib
R=Path(__file__).resolve().parent
for line in (R/'manifest.sha256').read_text().splitlines():
 digest,name=line.split(maxsplit=1);assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
old=json.loads((R/'evidence/uniprot/Q07N61.json').read_text());new=json.loads((R/'evidence/uniprot/Q6N0W9.json').read_text());pdb=json.loads((R/'evidence/pdb_4LF1.json').read_text())
assert 'BisA53' in old['organism']['scientificName'] and 'CGA009' in new['organism']['scientificName']
assert pdb['entity_src_gen'][0]['gene_src_strain']=='CGA009'
assert pdb['entity_poly']['pdbx_seq_one_letter_code_can']=='MGSSHHHHHHSSGLVPRGSH'+new['sequence']['value']
assert len(new['sequence']['value'])==461 and old['sequence']['value']!=new['sequence']['value']
fasta=''.join((R/'results/replacement_native_core_v1.fasta').read_text().splitlines()[1:]);assert fasta==new['sequence']['value']
ops=json.loads((R/'results/mapping_patch_v1.json').read_text())['operations']
assert sum(x['action']=='replace_native_core' for x in ops)==1 and sum(x['action']=='quarantine' for x in ops)==9
hist=(R/'evidence/unisave_A0A510B6V9_v17.txt').read_text();assert 'OS   Avena nuda' in hist and 'NCBI_TaxID=4497' in hist
with (R/'results/exposure_records.csv').open() as f:records=list(csv.DictReader(f))
with (R/'results/exposure_impact.csv').open() as f:impact=list(csv.DictReader(f))
for row in impact:
 data=[x for x in records if all(x[k]==row[k] for k in ['target','arm','convention'])]
 assert len(data)==int(row['prediction_records']) and sum(x['exposed']=='True' for x in data)==int(row['exposed_records'])
print('C07 public PASS: exact accession/deposition native-core repair, patch counts and derived exposure arithmetic. Primary-assay transcription and full input audit require separate source acquisition; see README.')

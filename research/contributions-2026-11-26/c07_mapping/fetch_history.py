#!/usr/bin/env python3
from fetch import get,ROOT
import json,concurrent.futures

def one(p):
 a=p.stem;log=[]
 x=get(('history/'+a+'_versions.json','https://rest.uniprot.org/unisave/'+a+'?format=json'));log.append(x)
 try:
  d=json.loads((ROOT/x['file']).read_text());v=d['results'][0]['entryVersion']
  log.append(get(('history/'+a+'_v'+str(v)+'.txt','https://rest.uniprot.org/unisave/'+a+'?format=txt&versions='+str(v))))
 except Exception as e:log.append({'accession':a,'parse_error':str(e)})
 return log
if __name__=='__main__':
 ps=[p for p in (ROOT/'evidence/uniprot').glob('*.json') if json.loads(p.read_text()).get('entryType')=='Inactive']
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:r=[x for y in pool.map(one,ps) for x in y]
 (ROOT/'evidence/history_retrieval.json').write_text(json.dumps(r,indent=2));print('history downloads',len(r),'success',sum(x.get('http_status')==200 for x in r))

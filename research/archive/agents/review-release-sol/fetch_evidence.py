#!/usr/bin/env python3
"""Read-only public evidence retrieval; failures are separate from model refusals."""
from pathlib import Path
import urllib.request,json,hashlib
O=Path('web_evidence');O.mkdir(exist_ok=True)
URLS={
'paper':'https://www.mdpi.com/3042-6723/1/2/7',
'notes':'https://www.mdpi.com/3042-6723/1/2/7/notes',
'issues_api':'https://api.github.com/repos/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/issues?state=all&per_page=100',
'commit':'https://api.github.com/repos/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/commits/ade87f378aa25c24a24c19915722aa04be92936e',
'pin_tree':'https://api.github.com/repos/jamesyoung93/Rubisco_Protein_Language_Model_Active_Learning/git/trees/ade87f378aa25c24a24c19915722aa04be92936e?recursive=1',
'zenodo':'https://zenodo.org/api/records/19546367',
'challenge':'https://millenniumproblems.bio/'}
def main():
 log=[]
 for name,url in URLS.items():
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'IndependentResearchAudit/1.0'})
   with urllib.request.urlopen(req,timeout=20) as r:b=r.read();status=r.status
   (O/(name+('.json' if name not in ['paper','notes','challenge'] else '.html'))).write_bytes(b)
   log.append(dict(name=name,url=url,status=status,bytes=len(b)))
  except Exception as e:log.append(dict(name=name,url=url,error=str(e)))
 (O/'retrieval_log.json').write_text(json.dumps(log,indent=2));print(log)
if __name__=='__main__':main()

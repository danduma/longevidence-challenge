import urllib.request,pathlib,json
P=pathlib.Path('sources'); log=[]
for suffix in ['4_ESM.xlsx','5_ESM.zip','1_ESM.pdf']:
 name='41467_2025_60622_MOESM'+suffix
 url='https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-60622-7/MediaObjects/'+name
 try:
  r=urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=25);h=dict(r.headers); print(name,h.get('Content-Length'),flush=True)
  log.append(dict(url=url,headers=h))
  if int(h.get('Content-Length','0'))>100_000_000:continue
  (P/name).write_bytes(urllib.request.urlopen(url,timeout=60).read())
 except Exception as e:
  log.append(dict(url=url,error=str(e)));print(e,flush=True)
(P/'source_data_access.json').write_text(json.dumps(log,indent=2))

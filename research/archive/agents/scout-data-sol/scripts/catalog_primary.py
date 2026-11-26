import sys,urllib.parse,json,concurrent.futures
from pathlib import Path
sys.path.insert(0,'scripts');from acquire import fetch
queries={
'origins2026':'DOI:10.1126/science.adt2760',
'origins2025':'DOI:10.1038/s41557-025-01830-y',
'ribozyme_design2025':'DOI:10.1038/s41467-025-63151-5',
'cryo2023':'DOI:10.1038/s41467-023-38824-8',
'reverse2026':'DOI:10.1038/s41587-026-03061-z',
'reverse_steps2026':'TITLE:"From peptides to DNA: All required steps can be catalyzed"',
'quad2021':'DOI:10.1038/s41467-021-25948-y',
'quad2025':'DOI:10.1093/nar/gkae1292',
'quad_mito2024':'DOI:10.1021/acssynbio.4c00630',
'limb2023':'DOI:10.1038/s41467-023-41944-w',
'limb2025':'DOI:10.1038/s41467-025-59497-5',
'aav2022':'DOI:10.1021/acssynbio.2c00451',
'aav2025':'DOI:10.1016/j.ymben.2025.03.011',
'protease_cys2026':'DOI:10.1038/s41586-026-11143-w',
'protease_metallo2026':'DOI:10.1038/s41586-026-11167-2',
'cleavenet2026':'DOI:10.1038/s41467-025-67226-1',
'binder2026':'DOI:10.1038/s41589-026-02199-w',
'delivery2026':'DOI:10.1038/s41467-026-70054-6',
'capa2024':'DOI:10.1021/acschembio.4c00211',
'peptide_amp2026':'DOI:10.1002/anie.202518911',
'polymerase2026':'DOI:10.1261/rna.080815.125',
'nitrogenase2025':'EXT_ID:PMC12425478',
'nfembed2025':'EXT_ID:PMC12417089',
'rubisco2019':'DOI:10.1021/acs.biochem.9b00237',
'rubisco2021':'DOI:10.1093/molbev/msab079',
'rubisco_aquatic2026':'DOI:10.1093/plphys/kiag507',
'rubisco_cyano2026':'DOI:10.1093/plphys/kiag259',
'rubisco_engineering2026':'DOI:10.1038/s41467-026-73783-w',
}
def run(item):
 name,q=item
 p=fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/search?query='+urllib.parse.quote(q)+'&format=json&resultType=core','sources/catalog_'+name+'.json')
 if not p:return
 j=json.loads(p.read_text()); r=j.get('resultList',{}).get('result',[])
 if r:
  x=r[0]; print(name,x.get('title'),x.get('doi'),x.get('firstPublicationDate'),x.get('pmcid'),flush=True)
  if x.get('pmcid') and x.get('isOpenAccess')=='Y':fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/'+x['pmcid']+'/fullTextXML','sources/'+name+'.xml')
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:list(e.map(run,queries.items()))

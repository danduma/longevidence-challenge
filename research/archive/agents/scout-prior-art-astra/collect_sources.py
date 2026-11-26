"""Resolve the article DOIs against Crossref, recording titles and publication dates."""
from concurrent.futures import ThreadPoolExecutor
from html import unescape
import json
from pathlib import Path
import re
from urllib.request import urlopen

DOIS = """
10.1021/acs.biochem.9b00237 10.1093/molbev/msab079 10.1093/jxb/erab356
10.1038/s44318-025-00419-y 10.1073/pnas.2321050121 10.1038/s41586-024-08455-0
10.1073/pnas.2501433122 10.1016/j.bbrc.2025.151940 10.1093/oxfordjournals.pcp.a028949
10.1006/bbrc.1997.6497 10.1093/jxb/erm283 10.1104/pp.16.01846 10.1104/pp.16.00750
10.1038/s41467-025-63483-2 10.1038/s41467-023-38824-8 10.1073/pnas.2516848123 10.1038/1731136a0
10.7554/eLife.76941 10.1038/s41467-021-25948-y 10.1093/nar/gkz745 10.1021/acssynbio.4c00630
10.1126/science.adt2760 10.1038/s41467-026-69531-9 10.1038/s41587-026-03061-z 10.1126/science.aed1656
10.1126/sciadv.abj2164 10.1038/s41467-026-72066-8 10.1021/acssynbio.2c00451 10.1038/s41586-026-11167-2
10.1038/ncomms15090 10.1038/s41586-025-09721-5 10.1038/382525a0 10.1021/ja026024i
10.1038/s41467-026-73783-w 10.64898/2026.10.07.756990 10.1080/15476286.2017.1372076 10.1126/sciadv.1501397 10.1038/srep20086 10.1038/s41467-025-67423-y
""".split()


def lookup(doi):
    url = 'https://api.crossref.org/works/' + doi
    for attempt in range(3):
        try:
            with urlopen(url, timeout=12) as response:
                record = json.load(response)['message']
            break
        except OSError:
            if attempt == 2:
                raise
    assert record['DOI'].lower() == doi.lower()
    date = record.get('published-online', record['published'])['date-parts'][0]
    title = unescape(re.sub('<[^>]+>', '', record['title'][0]))
    return {'doi': doi, 'title': title, 'date': '-'.join(f'{n:02d}' for n in date),
            'first_author': record.get('author',[{}])[0].get('family')}


if __name__ == '__main__':
    output = Path('bibliography.json')
    cached = {r['doi'].lower(): r for r in json.loads(output.read_text())} if output.exists() else {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for r in pool.map(lookup, [doi for doi in DOIS if doi.lower() not in cached]):
            cached[r['doi'].lower()] = r
    records = [cached[doi.lower()] for doi in DOIS]
    Path('bibliography.json').write_text(json.dumps(records, indent=2, ensure_ascii=False)+'\n')

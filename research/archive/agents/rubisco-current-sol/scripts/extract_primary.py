from pathlib import Path
import sys,json,hashlib
from html.parser import HTMLParser
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'vendor'))
from pypdf import PdfReader
class Text(HTMLParser):
 def __init__(self):super().__init__();self.parts=[]
 def handle_starttag(self,t,a):
  if t in ['p','h1','h2','h3','h4','tr','li']:self.parts.append('\n')
 def handle_endtag(self,t):
  if t in ['p','h1','h2','h3','h4','tr','li']:self.parts.append('\n')
 def handle_data(self,s):self.parts.append(s)
for f in (ROOT/'raw').glob('*.html'):
 if f.stat().st_size<10000:continue
 p=Text();p.feed(f.read_text(errors='replace'));(ROOT/'notes'/f'{f.stem}_readable.txt').write_text(''.join(p.parts))
for name in ['whitney2001','cyano2023','kubien2008']:
 path=ROOT/'raw'/f'{name}.pdf'
 if not path.exists() or not path.read_bytes().startswith(b'%PDF'):continue
 r=PdfReader(path)
 for i,p in enumerate(r.pages):
  (ROOT/'notes'/f'{name}_page{i+1:02}.txt').write_text(p.extract_text(extraction_mode='layout'))

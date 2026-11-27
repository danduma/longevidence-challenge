"""Read worksheet names and cells using only the Python standard library."""
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'


def worksheets(path):
    with ZipFile(path) as archive:
        def xml(name):
            return ET.fromstring(archive.read(name))

        shared = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            shared = [''.join(t.text or '' for t in si.findall('.//s:t', NS))
                      for si in xml('xl/sharedStrings.xml')]
        targets = {link.get('Id'): link.get('Target')
                   for link in xml('xl/_rels/workbook.xml.rels')}
        for sheet in xml('xl/workbook.xml').findall('.//s:sheet', NS):
            target = targets[sheet.attrib[REL]]
            member = target.lstrip('/') if target.startswith('/') else 'xl/' + target
            cells = {}
            for cell in xml(member).findall('.//s:c', NS):
                value = cell.findtext('s:v', namespaces=NS)
                if cell.get('t') == 's' and value is not None:
                    value = shared[int(value)]
                if cell.get('t') == 'inlineStr':
                    value = ''.join(t.text or '' for t in cell.findall('.//s:t', NS))
                if value is not None:
                    cells[cell.attrib['r']] = value
            yield sheet.attrib['name'], cells


if __name__ == '__main__':
    for path in sorted(Path('sources').glob('*.xlsx')):
        print(path)
        for name, cells in worksheets(path):
            print(name, len(cells), list(cells.items())[:12])

"""Read actual 2025 source workbook without external dependencies.

Source: Gangwar et al., DOI 10.1038/s41467-025-63483-2, Fig. 7.
Repeated probes are not treated as independent biological/experimental n.
The adiabatic two-region calculation is a feasibility screen, not validation.
"""
import csv
import json
import statistics
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def sheet_cells(path, sheet_index):
    with zipfile.ZipFile(path) as z:
        ss=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',NS)]
        r=ET.fromstring(z.read(f'xl/worksheets/sheet{sheet_index}.xml'))
        out={}
        for c in r.findall('.//m:c',NS):
            v=c.find('m:v',NS)
            if v is None:continue
            s=v.text
            if c.attrib.get('t')=='s':s=ss[int(s)]
            out[c.attrib['r']]=s
        return out


def main():
    cells=sheet_cells(ROOT/'data/cryo_liter_2025_source.xlsx',4)
    rows=[];summary={}
    cols={'Left':['M','N','O'],'Center':['P','Q','R'],'Right':['S','T','U']}
    for vol,row,fe,dr in [('1L',6,10.7,13),('2L',7,4.6,14)]:
        reps=[]
        for rep in range(3):
            rates=[]
            for probe,cc in cols.items():
                rate=float(cells[f'{cc[rep]}{row}']);rates.append(rate)
                rows.append({'volume':vol,'Fe_mg_per_mL':fe,'replicate':rep+1,'probe':probe,'warming_C_per_min':rate})
            reps.append(statistics.mean(rates))
        dt=[float(cells[f'{col}{dr}']) for col in ['M','N','O','P','Q','R']]
        rates=[r['warming_C_per_min'] for r in rows if r['volume']==vol]
        summary[vol]={'replicate_means_C_per_min':reps,'mean_C_per_min':statistics.mean(reps),
                      'sd_between_3_replicate_means':statistics.stdev(reps),
                      'minimum_probe_rate_C_per_min':min(rates),'maximum_probe_rate_C_per_min':max(rates),
                      'rate_per_Fe_concentration':statistics.mean(reps)/fe,
                      'delta_T_values_C':dt,'delta_T_mean_C':statistics.mean(dt),'delta_T_max_C':max(dt),
                      'min_rate_divided_by_nominal_M22_CWR_0_4':min(rates)/0.4}
    with (ROOT/'data/cryo_fig7_rates.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    # Cold region must reach -40 C from -150 C; hot region must stay below 0 C.
    # If r_hot/r_cold > 150/110, changing a common scalar RF duty cannot help
    # when transport is negligible and local source ratios are constant.
    ratio_limit=150/110
    examples=[{'rate_ratio_hot_over_cold':r,'hot_temperature_when_cold_reaches_minus40_C':-150+r*110,
               'meets_selected_hot_ceiling':r<=ratio_limit} for r in [1,1.2,1.36,1.5,2]]
    out={'source_doi':'10.1038/s41467-025-63483-2','summary':summary,
         'ratio_2L_to_1L_rate_per_Fe':summary['2L']['rate_per_Fe_concentration']/summary['1L']['rate_per_Fe_concentration'],
         'adiabatic_two_region_screen':{'T_start_C':-150,'cold_exit_C':-40,'hot_ceiling_C':0,
            'maximum_hot_cold_rate_ratio':ratio_limit,'examples':examples,
            'assumptions':'selected endpoints; negligible conduction; no latent heat; constant local source ratios and heat capacity; common duty schedule'},
         'warning':'Bulk M22 solution, not organs. Nominal CWR is not a voxel-specific measured threshold. No viability inference.'}
    (ROOT/'outputs/cryo_source_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()

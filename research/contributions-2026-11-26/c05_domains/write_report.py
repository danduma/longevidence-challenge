from pathlib import Path
import csv,json
O=Path('results');rows=list(csv.DictReader((O/'domains.csv').open()));v=json.loads((O/'verification.json').read_text());s=json.loads((O/'summary.json').read_text())
def get(case,domain):return next(r for r in rows if r['case']==case and r['domain']==domain)
text='''# C05 — tissue versus bath constraint attribution

C05 found an informative null for domain-dependent feasibility: **requiring the surrounding bath as well as myocardium changed no thermal decision or paired schedule ranking in the predefined comparison.** Attribution still matters: myocardium set every whole-specimen cold endpoint; surrounding CPA set every hot endpoint and most restricted-window rate minima. Excluding the bath changes reported margins and locations, but does not rescue a failing myocardial screen here.

Status: agent-produced; awaiting coordinator verification. Executed 2026-11-26; no human expert review. Computational thermal results only.

The new contribution is the systematic compartment ledger, full-cell failure/crossing accounting, and the tested domain-comparison null for frozen schedules. Reproducing the previously accepted ramp/flat result is a prerequisite, not claimed as new. Equal-iron spatial testing, staged heating, and temperature-window concepts are prior art. The rejected measured-heart triplicate inference is absent.

The frozen [protocol](protocol.md), recorded before new scoring, defines the exact question, comparisons and falsification. “Original” means the released **margin ramp**, not the earlier non-margin optimum. Flat retains that ramp's surrounding iron field and uses its time-average commanded duty, 0.9568032351. Each delivers 66.97622646 duty-seconds over 70 s; absorbed energy is not forced equal. Myocardium is uniformly nominal 1.47 or high 1.83 mg Fe/mL; chamber CPA is 10. Both low-temperature SAR readings are tested. No schedule or loading was tuned to these results.

The conservative finite-volume midpoint RK2 solver integrates wall, surrounding CPA, myocardium, and chamber CPA in every run. Primary screening requires local instantaneous rate ≥50 °C/min throughout 0–70 s and every terminal temperature between −35 and 0 °C. Secondary diagnostics restrict the rate screen to local −100 to −40 °C or −123 to −35 °C, require complete crossings of the respective interval, and keep the same endpoint constraints. Stage states, midpoint states and the final RHS contribute to instantaneous minima. The complete trajectories, not mean temperature, determine crossing times by linear interpolation.

C05's primary comparison is label 2 (myocardium) versus labels 1+2+3 (whole specimen). The retained diagnostic domains are wall 0, bath 1, chamber 3, and heart contents 2+3. Upstream `specimen_mask` explicitly excludes wall, while `active_mask` includes it. This distinction comes from the [pinned Vitrify geometry](https://github.com/Marc-Dvci/Vitrify/blob/1d498ee5c583e0112a762e60892d576ec2e27f57/src/vitrify/geometry.py), and the thermal implementation uses the specimen mask for reported minima. The zero-iron wall starts at zero warming rate and would fail an all-time ≥50 rule in every case if included. It remains in heat transport and energy accounting; excluding it from specimen screening is an explicit source definition, not removal of an adverse simulation domain.

[Joshi & Rabin 2023](https://doi.org/10.1371/journal.pone.0290063) models CPA-filled chambers and a heart immersed in surrounding CPA, with different iron concentrations by domain. It assigns CPA properties to tissue because tissue-specific properties were unavailable, and discusses differing tissue/pure-CPA critical rates. Therefore applying the same numerical threshold to both domains is a declared benchmark assumption, not proof of their biological equivalence. The present geometry is Vitrify's synthetic matched-dimension/volume phantom, not the article's MRI anatomy. [Sources, permissions and downloads](sources.md).

Discovery used all eight combinations at 0.75 mm, dt≤0.01 s. The predefined 0.5-mm subset used both schedules at high/clamped and nominal/extrapolated settings, at dt≤0.01 and ≤0.005 s: 16 deterministic integrations total, not 16 independent experiments. Every integration has six domain summaries. Spatial cells are quadrature volumes, not biological replicates. There are four discovery schedule pairs and four additional pairs from spatial/time checks, not eight independent observations.

Discovery results below give myocardial / whole-specimen minima in °C/min. Cold is myocardial (also the specimen minimum); hot is whole-specimen surrounding CPA. All-time and window decisions are identical under both domain definitions. The −123 to −35 °C decisions also coincide.

| Loading / SAR | Schedule | All-time, muscle / whole | −100 to −40, muscle / whole | Cold / hot °C | All-time / window |
|---|---|---:|---:|---:|---|
'''
for r in rows:
 if r['mm']!='0.75' or r['domain']!='myocardium':continue
 w=get(r['case'],'whole_specimen')
 text+=f"| {r['mean']} / {r['mode']} | {r['design']} | {float(r['all_rate']):.3f} / {float(w['all_rate']):.3f} | {float(r['window_rate']):.3f} / {float(w['window_rate']):.3f} | {float(w['cold']):.3f} / {float(w['hot']):.3f} | {'pass' if r['all_pass']=='True' else 'fail'} / {'pass' if r['window_pass']=='True' else 'fail'} |\n"
text+='''
No decision differs in any of the three screens over all 16 integrations (also true for adding only surrounding CPA to heart contents 2+3; see `results/bath_addition.csv`); no paired rate-ranking sign differs over all eight numerical pairs. Across discovery, whole-specimen all-time minima lie in myocardium in five runs and surrounding CPA in three; all restricted-window minima lie in surrounding CPA. Every cold endpoint lies in myocardium and every hot endpoint lies in surrounding CPA. Chamber CPA sets none of these whole-specimen extrema, although its own endpoints/transits are retained. These are descriptions of this finite design set, not universal rules.

At nominal/extrapolated/flat, discovery's myocardial all-time minimum is 54.180882 versus whole-specimen 52.767221 °C/min, the largest domain difference (1.413660). At 0.5 mm, dt=.005 s, those minima are 53.211295 and 52.681423, respectively. Rankings still coincide: flat wins the all-time comparison at nominal/extrapolated; original wins both restricted-temperature comparisons. High/extrapolated retains the prior exception in which original also wins all-time.

The refined high/clamped case shows the domains separately (0.5 mm, dt=.005 s). Rates and local transit ranges are °C/min; endpoints are °C. Transit rates are interval width divided by each cell's crossing duration, not the instantaneous minimum.

| Schedule / domain | All-time min | −100 to −40 min | −123 to −35 min | Terminal cold / hot | −100 to −40 transit min–max |
|---|---:|---:|---:|---:|---:|
'''
for d in ['original','flat']:
 for domain in ['myocardium','surrounding_CPA','chamber_CPA','whole_specimen']:
  r=get(f'0.5_{d}_1.83_clamped_0.005',domain)
  text+=f"| {d} / {domain} | {float(r['all_rate']):.4f} | {float(r['window_rate']):.4f} | {float(r['glass_rate']):.4f} | {float(r['cold']):.4f} / {float(r['hot']):.4f} | {float(r['window_transit_rate_min']):.4f}–{float(r['window_transit_rate_max']):.4f} |\n"
text+='''
Locations are not unique when symmetry or a uniform start ties cells. The detailed table records representative voxel indices, physical voxel-center coordinates, extremum time/temperature, tie counts and rate-tie bounding indices. Coordinates use the phantom array origin, in mm; they are not anatomical coordinates. For refined high/clamped/original:

- The 44.435757 °C/min all-time minimum occurs at −150 °C, t=0, throughout all 3,700 myocardial cells. One representative is (4.75, 8.25, 5.75) mm; it is not a uniquely identified cold spot.
- The −100 to −40 minimum is 75.764081 °C/min in surrounding CPA at (11.75, 8.75, 1.25) mm, t=61.3525 s, T=−40.002311 °C. The broader-window minimum is 71.479443 there at t=65.43 s, T=−35.000914 °C.
- The terminal minimum is −32.996510 °C in myocardium at (11.75, 8.75, 2.75) mm; the maximum is −12.845650 °C in surrounding CPA at (3.25, 8.75, 6.25) mm, both at 70 s.

Exact bottleneck labels are not uniformly mesh-stable. In nominal/extrapolated/flat at 0.5 mm, dt=.005 s, the −100 to −40 minimum is myocardial (71.613167), just below surrounding CPA (71.614939). At 0.75 mm it is in surrounding CPA. Thus C05 supports the observed decision concordance but does not assert that the bath always sets the restricted-rate minimum. The refined high/clamped flat bath–muscle gap is also tiny (0.000430 °C/min); fine distinctions in tied or nearly tied locations should not be interpreted as anatomical discoveries.

Failure volumes and crossing accounting are essential to the null. At 0.75 mm myocardium is 0.461109375 mL and specimen 2.19965625 mL. Both nominal/clamped schedules have all myocardial cells below 50 °C/min at some point, so their all-time union failure volume is exactly 0.461109375 mL. Only myocardial cells fail the cold endpoint: 0.013921875 mL (33 cells) for original and 0.010968750 mL (26 cells) for flat. These same cells have not reached −35 °C, and are explicitly incomplete in the −123 to −35 interval. They are not omitted to produce a favorable transit statistic. All have reached −40 °C, so their narrow-window endpoint failure remains separate from rate clearance.

High/clamped/original also fails the all-time rate over all myocardium (0.461109375 mL at discovery; 0.462500 mL at 0.5 mm), but clears the endpoint and both restricted-rate screens. High/clamped/flat has zero failed specimen volume for all three screens. All other discovery combinations have zero failed specimen volume. Surrounding and chamber CPA have zero all-time/window/endpoint failure volumes in this set. No specimen cell exceeds 0 °C at the terminal time. Wall failures are separately tabulated, not counted as specimen failures.

Every nonwall cell in every run completes the −100 to −40 crossing. The only incomplete −123 to −35 specimen crossings are the nominal/clamped discovery cells just listed. Wall crossing accounting is also retained. Complete-only transit ranges are labeled alongside missing counts and volumes in the table; they must not be read as full-domain traversal when an exit is missing.

Verification was executed with the shared scientific Python using pinned local inputs, no model calls or network dependency:

```sh
./verify.sh
# Portable alternative:
PYTHON=/absolute/path/to/python ./verify.sh
```

The verifier compiles the local C++17 solver, checks protocol/input hashes, runs analytic checks, reproduces the central baseline before new scoring, executes all model comparisons, and regenerates domain/ranking ledgers. The central high/clamped pair matches the pinned independent RK2 rates/endpoints to maximum absolute error 3.65×10⁻¹³ (tolerance 0.02). This is a code baseline reproduction, not independent experimental validation.

The largest relative enthalpy/source-energy imbalance across full-model runs is 1.230×10⁻⁷, including wall energy. The scalar nonlinear-heating quadrature differs by 0.0000635 s at dt=.001 s; an independent −123 °C crossing calculation differs by 2.49×10⁻¹⁰ s. The unit check validates q=1000 C_Fe SAR u and the conversion to °C/min. Zero-source uniform storage is preserved and symmetric heated cells remain equal. All observed rates are nonnegative; crossing order and complete/missing counts reconcile.

Halving the full-model time step at 0.5 mm changes any reported instantaneous minimum by at most 0.002872 °C/min and either endpoint by at most 0.00001481 °C. None of the predefined subset's domain decisions or rankings changes with mesh or time step. Endpoint temperatures do change materially with voxelization (e.g. high/clamped/original cold −34.752 at 0.75 mm versus −32.997 at 0.5 mm); C05 does not claim continuum convergence. Nominal/clamped was discovery-only, as predefined, so its cold-failure volumes have no new 0.5-mm confirmation here.

The full evidence is in [domains.csv](results/domains.csv), [paired comparisons](results/paired.csv), [domain decisions](results/domain_comparison.csv), [domain rankings](results/domain_rankings.csv), [constraint attribution](results/attribution.csv), [time refinement](results/time_refinement.csv), [baseline checks](results/baseline_check.csv), and [verification](results/verification.json). Per-cell CSVs and exact solver input snapshots remain local under `results/`; these permit independent re-aggregation of any mask. Six overlapping domains do not constitute independent replicates. Versions and source hashes accompany the tables in [provenance.json](provenance.json) and `manifest.sha256`.

C05 does not establish that it is appropriate to ignore surrounding CPA biologically. The bath still conducts heat into tissue and can matter mechanically; neither mechanism is tested by changing an output mask. Source-defined whole-specimen screening remains the primary benchmark, with myocardium-only as an explicitly isolated thermal diagnostic. The uniform loading, phantom geometry, SAR extrapolation/clamp, CPA property substitution, adiabatic exterior and numerical criteria are model assumptions. Co-registered tissue loading/CPA/thermometry and outcome measurements are unavailable. No ice fraction, stress, fracture, toxicity, survival, organ recovery or probability of failure follows. No human expert has reviewed this new work. [Source-access failures and development exceptions](refusals.md) are recorded separately from scientific negative results.
'''
Path('report.md').write_text(text)
claims=[dict(name='domain_decision_changes',value=0,unit='of 16 deterministic integrations per screen',evidence='results/domain_comparison.csv'),dict(name='domain_ranking_changes',value=0,unit='of 8 numerical schedule pairs per rate screen',evidence='results/domain_rankings.csv'),dict(name='max_myocardial_minus_specimen_all_time_minimum',value=s['max_myocardium_minus_whole_all_rate'],unit='C/min',evidence='results/summary.json'),dict(name='max_energy_relative_error',value=v['max_energy_rel_error'],unit='relative',evidence='results/verification.json'),dict(name='max_time_refinement_rate_change',value=v['max_time_rate_delta'],unit='C/min',evidence='results/time_refinement.csv'),dict(name='max_time_refinement_endpoint_change',value=v['max_time_endpoint_delta'],unit='C',evidence='results/time_refinement.csv')]
c=dict(id='C05',title='Tissue versus bath constraint attribution in the pinned Vitrify phantom',question='Does requiring surrounding CPA as well as myocardium change frozen schedule rankings or thermal feasibility, and which domains set the constraints?',method='Predefined 2 schedules x 2 loadings x 2 SAR readings at .75 mm; high/clamped and nominal/extrapolated pairs at .5 mm with .01/.005 s exact-horizon conservative RK2; six domain masks with full-cell minima, endpoint, crossing and failure-volume diagnostics.',result='Informative null: no domain-dependent decision or ranking changes. Myocardium sets all cold endpoints; surrounding CPA sets all hot endpoints and most restricted-window minima. Exact rate-bottleneck attribution can change with voxelization.',verification_command='./verify.sh',numerical_claims=claims,primary_sources=[dict(title='Vitrify pinned release',url='https://github.com/Marc-Dvci/Vitrify/tree/1d498ee5c583e0112a762e60892d576ec2e27f57',version='1d498ee5c583e0112a762e60892d576ec2e27f57',retrieved='2026-11-26'),dict(title='Joshi and Rabin 2023',url='https://doi.org/10.1371/journal.pone.0290063',version='published 2023-08-16',retrieved='2026-11-26')],limitations=['Synthetic geometry, uniform myocardial loading, substituted CPA properties, unresolved SAR reading and adiabatic boundary.','Deterministic parameter scenarios and overlapping domain masks are not independent replicates.','Near-tied spatial attribution and endpoint values depend on voxelization; only predefined subset refined.','No measured loading/CPA/temperature maps, ice, mechanics or biological validation; isolated tissue screen is not biological safety.','No human expert review; no universal novelty claim.'],status='agent-produced; awaiting coordinator verification')
Path('contribution.json').write_text(json.dumps(c,indent=2)+'\n')

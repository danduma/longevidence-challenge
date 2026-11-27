# Evidence needed to reopen questions

These are reviewed data-gap records, not contributions toward the target. Suggested producers are capability/interest matches, not partners; no outreach or laboratory work has been performed.

## rubisco-assay-identifiability

The Rubisco calibration needs same-preparation assay pairs at 25 and 30 °C

**Later progress:** the pRR1 clone-lineage follow-up supplied supporting source-identity evidence after this gap was reviewed. It does not establish the identity of the historical assay stock or supply crossed assay-temperature measurements. The original gap and its review are preserved below; see the [Rubisco audit](audits/rubisco/report.md) for the qualified lineage result.

### Potential producers

- **Milo Lab, Weizmann Institute of Science** — Published the screen and Rubrum control normalization. Fit: Likely best placed to recover original construct and assay identifiers and reproduce the coupled arm. Confidence: High for published screen capability. Fit is inferred. [Supporting source](https://pmc.ncbi.nlm.nih.gov/articles/PMC12663927/).
- **Oliver Mueller-Cajar and Tsai collaborators** — Published both microbial turnover formats and a radiometric active-site recalibration. Fit: Likely able to link the corrected radiometric method to a matched comparison. Confidence: High for demonstrated methods. Fit is inferred. [Supporting source](https://link.springer.com/article/10.1038/s44318-025-00419-y).

**Controls**

Match dissolved gases, saturation, activation, pH at assay temperature and Mg conditions. Record raw rates and a shared independent active-site estimate, retain the unnormalized Rr control and test for readout linearity or coupling limitations.

**Existing evidence**

Eleven reported literature overlaps at 25 °C versus the 30 °C screen. The geometric comparison is 1.436 before and 2.129 after fixed-Q10 adjustment. Data 2 provides 180 matching rows without a crossed assay-temperature field. McNevin provides a within-format 1.66-fold temperature contrast for a historical recombinant Rubrum preparation.

**Experimental unit**

Independent preparation of a verified construct, split across all conditions. Wells and titration points are technical observations within a preparation.

**Inconclusive**

Unresolved sequence identity, saturation or active-site discrepancies, or paired intervals spanning the prespecified material differences.

**Minimum design**

Measure direct and coupled turnover at 25 and 30 °C on split aliquots of independently prepared material. Include Rr and a sequence-verified nonreference anchor such as PCC 6301. Four cells identify the interaction. Three genuinely matched cells identify additive contrasts under a no-interaction assumption.

**Missing data**

Per-preparation raw carboxylation readouts and active-site estimates for both formats at both temperatures, with verified construct and condition correspondence.

**Precision basis**

Use pilot paired log-contrast variation across independent preparations and a prespecified smallest rate difference that would affect benchmark interpretation.

**Question**

What are the 25→30 °C changes within each format and the direct-to-coupled ratios at each temperature for the same microbial enzyme preparations?

**Refutes**

A material ratio-of-ratios departure from one or material construct heterogeneity would refute a shared temperature-independent correction over the sampled conditions. A precise within-format temperature ratio incompatible with sqrt(2.2) would refute that temperature adjustment locally.

**Reopen trigger**

Public per-preparation crossed turnover data or a verified exact-construct fixed-format temperature bridge with documented gas and normalization settings.

**Supports**

Tightly estimated compatible direct/coupled ratios at both temperatures, stable across tested preparations and constructs, would support a common format correction over that sampled range. Compare each within-format temperature ratio with sqrt(2.2).

**Why insufficient**

The published diagonal also varies gases, preparation and active-site measurement. The McNevin series constrains temperature within its own radiometric preparation, and the Davidi/Tsai format comparison changes temperature and protocol simultaneously.

**Work possible now**

Resolve the pRR1 lineage from public clone and sequence records and inspect the unretrieved numerical Synechococcus and Thermosynechococcus records if openly available.

**Review limitation:** The roughly 1.66-fold McNevin contrast comes from digitized points in one historical radiometric preparation. It lacks preparation-level uncertainty and a verified construct and protocol match to the 2025 reference. Davidi/Tsai still couples a format change with a temperature change. Even a confirmed pRR1 sequence match would strengthen only the within-format comparison. The [published pRR1 clone lineage](https://onlinelibrary.wiley.com/doi/pdfdirect/10.1016/0014-5793%2890%2980879-N) gives the offered follow-up a concrete starting point.

## cryo-measurement-bridge

Public warming data constrain a warm EMG308 reference but leave the matched-organ model inputs unresolved

### Potential producers

- **ATP-Bio / University of Minnesota organ nanowarming team** — Organ imaging, thermometry and solution calorimetry across kidney and liter-scale studies Fit: A specimen-keyed release could connect measurements within the team’s organ-banking work Confidence: High inferred capability fit. Existence of a matched package is unverified. [Supporting source](https://atp-bio.org/press-release-landmark-study-demonstrates-first-successful-long-term-cryopreservation-and-nanowarming-of-human-scale-kidneys/).
- **Until** — Publicly describes organ perfusion, nanoparticle delivery and rewarming development Fit: A paired model-validation release aligns with the stated technology direction Confidence: Moderate inferred fit. Specific assay and paired dataset capabilities are unverified. [Supporting source](https://www.untillabs.com/science-technology).

**Controls**

Calibration phantoms, unspiked Fe backgrounds, CPA-only RF and thermal controls, independent heat-transfer calibration, and matched particle preparation and medium for SAR and the warmed organ.

**Existing evidence**

Relative SWIFT with a cross-study dry-weight anchor, bulk MPI, adjacent-fluid and sparse core thermometry, a warm EMG308 vial reference and separate phase-specific solution calibrations.

**Experimental unit**

Whole organ per warming run, grouped by donor where paired. Voxels, probes and time samples are repeated measurements.

**Inconclusive**

Map, probe, SAR or boundary uncertainty spans the predicted differences, or the released records cannot be joined at specimen level.

**Minimum design**

Start with existing records if available. Couple a calibrated prewarming image to core, peripheral and bag thermometry in common coordinates with run metadata. Independently calibrate SAR and boundary terms. Hold out complete organs and compare against uniform loading at the same total Fe and a mass-matched two-region reconstruction when region volumes are known.

**Missing data**

A specimen-keyed absolute prewarming Fe map in mg/mL of bulk voxel, registered segmentation and probe locations, loaded regional dry mass per volume with Fe background, matched particle/CPA/field SAR across the observed temperature range, and measured bag loading, geometry, initial state and boundary history.

**Precision basis**

Propagate Fe assay, map, registration, probe, SAR and boundary uncertainty. Use the separation between candidate predictions for prespecified window crossings and gradients to set required precision.

**Question**

Can an absolute prewarming iron map and an independently calibrated heating law predict spatial thermometry on held-out intact organs?

**Refutes**

With adequate precision the map-based prediction misses measured gradients or crossings or fails to improve on the mass-matched baseline under the same fixed SAR and boundary terms.

**Reopen trigger**

Release of the 2021 segmented volumes and model assignments would reopen mass closure. A specimen-keyed package joining absolute loading, independent matched SAR, registered thermometry and boundary records would reopen the organ benchmark.

**Supports**

On held-out organs the registered map plus independently fixed SAR and boundary terms predicts prespecified spatial gradients and temperature crossings beyond propagated uncertainty and outperforms a mass-matched uniform baseline.

**Why insufficient**

Bulk Fe does not resolve tissue compartments. MRI and CT signals need particle-specific absolute calibration. Temperature depends on concentration, SAR, heat capacity and heat exchange. The phase-specific sIONP assay cannot calibrate the EMG308 warming traces. The matched EMG308 reference is a warm vial measurement, and the linked models are for cooling.

**Work possible now**

The matched warm EMG308 anchor supports the fixed-input bag diagnostic executed here. The released S22A equal-concentration traces can test whether temperature-window net heating transfers across bag volumes. That comparison remains conditional on heat capacity and cannot isolate particle SAR or validate the organ model.

**Review limitation:** The residual comes from an audit-selected −20 to 0 °C window in two example bag runs. The six probes are repeated measurements. Approximate heat capacity, vial-to-bag SAR transfer, direct RF heating and boundary exchange can all contribute. The authors' Fig. S20 estimates near-zero eddy heating of comparable magnitude under surrogate assumptions, so the residual cannot identify a heating term. The gate is limited to inspected public records. The porcine supplement remained unverified, and the [porcine preprint](https://www.biorxiv.org/content/10.64898/2026.10.07.756990v1.full) says raw data are available on request. This agent audit establishes neither human review nor biological validation.


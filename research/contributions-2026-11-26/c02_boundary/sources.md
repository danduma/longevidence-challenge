# C02 · Sources and reproduction

Sources inspected 2026-11-26. Code and inputs in pinned/ are byte copies from the local public challenge checkout at 01b6be05300e5e7878e51c77ab2f4ccb7658259d. Locate research/cryopreservation/independent.cpp, reproduce.py, vendor/ and inputs/ in that revision and LICENSE at its root. Compare downloaded bytes to manifest.sha256. The project MIT license covers the independent solver, including its local C02 instrumentation. Geometry and the six margin-design inputs derive from Vitrify at 1d498ee5c583e0112a762e60892d576ec2e27f57 under Apache-2.0. Both notices are included. The input JSON also contains earlier outputs; C02 reads only design_parameters. [Challenge snapshot](https://github.com/danduma/longevidence-challenge/tree/01b6be05300e5e7878e51c77ab2f4ccb7658259d/research/cryopreservation) · [Original geometry and design](https://github.com/Marc-Dvci/Vitrify/tree/1d498ee5c583e0112a762e60892d576ec2e27f57)

The original article XML was inspected from the shared snapshot dated 2026-11-26. Public XML can be downloaded via the links below. Article files are unnecessary to execute and are omitted. SHA-256 values identify the inspected copies.

- Gao et al., *Vitrification and Rewarming of Magnetic Nanoparticle-Loaded Rat Hearts*, published online 2021. Methods give the 22 °C and 15 W/(m² K) simulation boundary and separately describe experimental thermometry. XML SHA-256 `01161db73cb046c18fd5e39744aebe29a28f93e200c6b3e6969f44eab0826304`. [Primary paper](https://doi.org/10.1002/admt.202100873) · [Europe PMC XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9164386/fullTextXML)
- Joshi & Rabin (2023), *Thermomechanical stress analyses of nanowarming-assisted recovery from cryopreservation by vitrification in human heart and rat heart models*. Table 2 supplies constitutive fits and Table 3 gives modeled boundaries. XML SHA-256 `446d7960c218d3d16c0c2670cf30fb3d4cdbdf350f1366942f23380a50465472`. [Primary paper](https://doi.org/10.1371/journal.pone.0290063) · [Europe PMC XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10431620/fullTextXML)

Consulted the corrected millennium/FINDINGS.md, cryo-robust-astra/report.md and review-cryo-sol/report.md in the existing research/millennium tree. The latter's independent baseline table supplies C02's six numerical acceptance targets. Earlier results at h=0, 5 and 10 are prior work.

To reproduce, copy the C02 deliverables listed in manifest.sha256 and the manifest itself into a writable directory, preserving relative paths. Use Python 3.14.10 with NumPy 2.5.5, a C++17 compiler and shasum. The original execution used arm64 and Apple clang 21.0.0. Run `python3 -B verify.py` there. The verifier recompiles and computes both meshes in a temporary child directory, checks the frozen baselines first, then compares both complete numerical tables. Four local solver processes run concurrently by default; C02_WORKERS sets the count. All required numerical inputs are bundled. No model or network service is used by the verifier. For the additional boundary-unit fixture, run the compiler command in run.py (or run run.py) in this directory, then run `python3 -B unit_check.py`. For a faster stored-table and checksum audit run `python3 -B verify.py --tables-only`.

Local prior-report SHA-256 digests, relative to research/millennium, are recorded below.

```text
faa72dcb79787fd0f3533d247928b0b496325d314eabf47ccfe927b6d6b461af  FINDINGS.md
a0707e051a61d19b09aec90ba8bbef34b1c78bc5bc19c4937b80d98df5f33198  agents/cryo-robust-astra/report.md
536d3f28bdf3b834e9272115ca8f7c1c68479a5496cb7a2218bf0e9f6060703a  agents/review-cryo-sol/report.md
```

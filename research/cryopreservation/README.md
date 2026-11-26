# Independent thermal benchmark

Requires Python 3.11+, a C++17 compiler, NumPy and SciPy. No model account, SSH key, Docker, JAX or private data is needed.

```sh
python3 -m venv .venv
.venv/bin/pip install -r research/cryopreservation/requirements.txt
.venv/bin/python research/cryopreservation/reproduce.py
```

Four full-field cases independently integrate the nonlinear thermal equations using exact-horizon RK2. At 0.50 mm and uniform 1.83 mg Fe/mL, clamped SAR, the constant-duty all-time minimum is about 50.319 °C/min versus the ramp's 44.436; the ramp's local −100 to −40 °C minimum is higher. The second pair at nominal loading/extrapolated SAR also reverses the rate-screen ranking. Output is written under `results/`.

The six released design coefficients and analytic phantom are pinned to Vitrify commit `1d498ee5c583e0112a762e60892d576ec2e27f57`. Geometry is reproduced under Apache-2.0 (vendor/VITRIFY-LICENSE); equations follow Joshi & Rabin, DOI 10.1371/journal.pone.0290063. The separate solver and integration are this project's work. The review history and original extended results are in `../archive/agents/review-cryo-sol/` and `reference/`.

This is a conditional model result. The high-loading pass has a small margin, assumes an adiabatic exterior, changes the original fixed loading, and provides no ice/stress/toxicity/viability validation. The measured-triplicate empirical inference was rejected and is not part of this benchmark. General spatial loading, weakest-region metrics and staged heating are prior art. See the corrected findings and original primary references.

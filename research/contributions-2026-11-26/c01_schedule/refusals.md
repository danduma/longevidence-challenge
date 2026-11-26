# C01 — refusals, access failures and computational corrections

No model refusals occurred. No other models, subagents or people were contacted, and no external posting, commits, deployment or biological work was performed.

Source access failures (2026-11-26): Europe PMC model XML timed out; the already pinned, CC BY original model XML from the designated prior worker was copied locally as an explicit fallback. Europe PMC PMC8823202 fullTextXML returned HTTP 500; the PMC browser page returned a reCAPTCHA. These are source-service failures, not model refusals. The retrieved Gao heart paper independently states VS55 Tg = -123 C, and supports the temperature-window definitions. `inputs/retrievals.json` records individual outcomes and timestamps. Pinned upstream README and heart XML retrieval succeeded. Live search returned Vitrify; no exhaustive new priority search is claimed.

Numerical check development: the first scalar crossing-check fixture used Fe=3 mg/mL and did not reach -40 C within 70 s, so an assertion requiring that crossing failed. Inspection showed the solver correctly marked the crossing incomplete. The validation fixture was changed to Fe=6 mg/mL to test completed crossings against independent quadrature. This changes neither the family nor any scored scenario. The completed checks and actual errors are in `results/numerical_checks.json`.

No measured-heart triplicate failure inference is used. Unavailable anatomical geometry, local delivery maps, measured SAR histories, ice kinetics, mechanics and recovery outcomes remain unsupported; they were not substituted with invented observations.

After all 66 trajectories completed, the initial verifier's onset-identity assertion exposed a rounded heat capacity in the Python analytic metadata (1,442,178 instead of the exact 1,442,177.75 J/m³/K). The C++ solver already used the exact correlations. `run.py` and the derived analytic metadata were corrected to the exact formula; thermal scores were unchanged. The corrected verifier passes. `results/run.log` retains the initial assertion failure; `results/final_checks.log` records subsequent verification. This is a check/metadata correction, not a source failure, model refusal or discarded outcome.

Execution-scope correction: a temporary solver-construction helper was initially created under `/tmp`, contrary to the requested job-only write scope. It was moved into `inputs/construct_solver.py`, made portable, and no task artifact remains at that temporary path. Shared inputs were never changed. This helper is not called during numerical reproduction.

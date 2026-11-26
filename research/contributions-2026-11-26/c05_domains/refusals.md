# C05 access, refusals and execution exceptions

No model refusals; no other model calls, subagents, people contacted, external posts, commits, deployments, or laboratory/animal/human work.

Source-access failures: browser tool could not fetch Europe PMC XML (internal error) or the immutable GitHub tree (restricted URL). Direct Python Europe PMC XML retrieval timed out. These are source/tool-access failures, not model refusals. Live raw GitHub geometry/workflow downloads succeeded. The PLOS publisher page and direct publisher XML download succeeded; original XML is pinned locally. See `sources/retrieval.json`.

Development checks found and corrected two erroneous verifier assumptions: the isolated heated cell's all-time minimum need not equal its initial rate (the unit check now uses the actual initial RHS), and not every cell necessarily completes the -123 to -35 C interval (the verifier now checks explicit complete/missing accounting). A strict scalar quadrature check at dt=.005 s missed its 0.0002-s tolerance by 0.000121 s; the scalar check was refined to .001 s. These do not change the preregistered full-model grids, schedules, criteria or time steps. Final fresh execution and numerical tolerances are recorded in `results/verification.json`. No unexplained execution failure is used as a scientific result.

Unavailable evidence: matched spatial myocardial loading/CPA maps, direct phantom thermometry, ice kinetics, mechanical solutions for both schedules and biological outcomes are absent. They limit interpretation, not the ability to perform this thermal attribution comparison. No human expert review of C05 has occurred.

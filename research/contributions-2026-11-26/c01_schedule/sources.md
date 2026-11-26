# C01 sources and reproduction

Geometry and design: Marc Davici, Vitrify, commit1d498ee5c583e0112a762e60892d576ec2e27f57, Apache-2.0; geometry code and notice are in vendor/. Original released parameters are in inputs/grid_convergence_margin.json. https://github.com/Marc-Dvci/Vitrify/tree/1d498ee5c583e0112a762e60892d576ec2e27f57

Thermal properties: Joshi and Rabin (2023), DOI10.1371/journal.pone.0290063, CC BY4.0, pinned article XML inputs/model.xml. Gao et al., DOI10.1002/admt.202100873, CC BY4.0, pinned XML inputs/heart.xml. These supply nominal assumptions, not empirical validation of this schedule family.

The earlier independent RK2 code and baseline are Longevidence MIT software; its license is retained. Input and code byte hashes, original source access outcomes, versions and URLs are in provenance.json and inputs/retrievals.json. original-worker-provenance.json preserves hashes before portable interpreter-path edits. sources_download.py performs optional source acquisition, which is not needed for offline reproduction and may retrieve new bytes at mutable endpoints.

Install requirements.txt, then run `PYTHON=python3 sh reproduce.sh` with a C++17 compiler. All66 cases, scalar checks and reports are regenerated from bundled inputs without model/network access. The source archive contains the earlier baseline as a control, not as this contribution's novelty. No ice, mechanics or recovery is established.

# Vitrify

**End-to-end inverse design of organ nanowarming across JAX and FEniCSx.**

Tesseract Hackathon 2026 · Track 02: Multi-physics & coupled systems

[![CI](https://github.com/Marc-Dvci/Vitrify/actions/workflows/ci.yml/badge.svg)](https://github.com/Marc-Dvci/Vitrify/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Demo film](https://img.shields.io/badge/demo-5%20min%20film-red)](https://youtu.be/Nqu64OvAIhg)

📺 **[Five-minute demo film](https://youtu.be/Nqu64OvAIhg)** — what the problem
is, why the gradient has to cross a container boundary, and what came out.

![Baseline, designed and difference stress fields](docs/figures/hero.png)

Vitrify designs physically realizable nanoparticle loading shells and RF duty
schedules for rewarming a vitrified organ. Its objective is mechanical — avoid
fracture — while its hard constraints are thermal: warm faster than the
critical warming rate everywhere, and finish inside the safe temperature
window. Neither solver can see both.

The result is one differentiable function spanning two incompatible scientific
stacks:

| Tesseract | Solver | Differentiation |
|---|---|---|
| `vitrify-thermal` | JAX finite volume, nonlinear transient conduction | exact forward/reverse AD through a checkpointed time loop |
| `vitrify-mechanics` | FEniCSx 0.11 + PETSc, 3-D Maxwell thermomechanics | exact analytic **tangent-linear model and adjoint**, sharing the forward factorisation |

Tesseract itself is forked, to let a tangent carry a batch of directions:
[`batched-jvp`](https://github.com/Marc-Dvci/tesseract-core/tree/batched-jvp),
which turns a Jacobian from twelve container calls into two.

![Architecture](docs/figures/architecture.png)

## Headline result

Against the physiologically achievable Case II protocol from Joshi & Rabin
(2023), a multi-start SQP driven entirely by end-to-end gradients produces:

| Metric | Case II | Designed |
|---|---:|---:|
| Peak tension | 1.234 MPa | **1.176 MPa** (−4.7%) |
| Peak compression | 1.384 MPa | **0.995 MPa** (−28.1%) |
| Slowest warming rate | 50.82 °C/min | **50.01 °C/min** |
| Final temperature range | −40.25 to −15.70 °C | **−35.00 to −7.47 °C** |
| Feasible? | **no** — cold spot below the melting margin | **yes** |

Nothing is warm-started by hand. Every run begins at a uniform shell at a
published loading and a constant duty cycle; the same gradients first *restore
feasibility* — the published protocol is infeasible on this phantom — and then
reduce stress. Seven of nine starts produced a feasible design and four land
within 0.7% of each other.

![Optimisation trace](docs/figures/optimization_trace.png)

The recommendation is one sentence: **load the surrounding shell close to
saturation and throttle the early RF power to about 82%.** Every converged
start agreed on that early duty cycle to within a few percent. Slowing the
first seconds of rewarming flattens the gradient across the
myocardium–chamber interface while the glass is still stiff enough to hold
stress, and the shell loading buys back the warming rate the throttle costs.

![Designed protocol](docs/figures/designed_protocol.png)

![Designed shell, rotating](docs/figures/design_rotation.gif)

![Stress through the dangerous window](docs/figures/stress_evolution.gif)

## Why Tesseract is load-bearing

The thermal solver is a JAX program with a checkpointed scan over thousands of
explicit steps. The mechanics solver is a FEniCSx/PETSc program that assembles
and factorises a changing equilibrium system at every saved thermal state.
Their dependency stacks do not share a practical runtime — the mechanics image
inherits PETSc, MPI and dolfinx from `dolfinx/dolfinx:stable` and has no JAX at
all — and their differentiation strategies are different.

Without Tesseract the options are to rewrite a validated mechanics code in JAX,
to abandon end-to-end gradients, or to hand-roll RPC and derivative plumbing.
Here both solvers stay native to their own tools while the workflow exposes a
single JVP, and the optimiser linearises the objective *and* all three physical
constraints through it.

## The mechanics side differentiates itself

The mechanics endpoint is not a finite difference. Every coefficient of the
Maxwell update, the assembled operator, the load vector and the smoothed
principal-stress objective carries an analytic derivative, so the tangent
problem

```text
K(T) u̇ = ḟ(T, Ṫ) − K̇(T, Ṫ) u
```

shares its operator with the forward problem: once PETSc has factorised `K(T)`,
the tangent is one more back-substitution rather than a second equilibrium
trajectory.

That also fixes the validation. A finite-difference endpoint checked against
finite differences agrees to ten digits and proves nothing, because both sides
are nearly the same expression. With an analytic tangent, finite differences
become an independent instrument again — and they converge to it at **second
order over three decades** (fitted order 2.00) before round-off takes over.

![Gradient validation](docs/figures/gradient_validation.png)

Across the whole two-container composition the best agreement is 1.0×10⁻⁸.
Finite differences are the unreliable side there: the published material
correlations are piecewise, so most step sizes push some cells across a fitted
seam and measure a secant across a kink.

## Validation before optimisation

Both solvers were gated against published results before the optimiser existed,
and both gates ship as `experiments/validate_forward_models.py`.

![Forward validation](docs/figures/forward_validation.png)

Thermal, at 0.5 mm: domain-average rewarming rates for the three unambiguous
published cases land within 7.5%, with the ordering reproduced. Mechanics, with
no fitted coefficient: within the pre-declared factor of two on absolute stress,
with the decisive signature correct — slow myocardium in tension
(+0.79 MPa mean maximum principal stress) against fast chambers in compression
(−0.42 MPa), and peak tension rising when loading moves from the unattainable
uniform case to the physiologically achievable one.

## What the result does and does not survive

Two checks that a coupled optimum on a coarse grid deserves.

**The optimiser cannot game the horizon.** The mechanics Tesseract sees only
the first 12 s of rewarming, which creates a way to lower the objective without
lowering stress: push the peak past the window. Re-evaluated with the window
more than doubled, both designs' peak value and peak timing are unchanged to
machine precision.

**The improvement is resolution-dependent in magnitude, not in sign.** The
design runs at 1.4 mm, where the optimum sits exactly on the 50 °C/min ice
constraint. The slowest cell in the domain is the least converged statistic a
coarse grid produces, so both designs were re-evaluated under refinement:

![Grid convergence](docs/figures/grid_convergence.png)

The reduction is positive at every level but ranges from 4.7% on the design
grid to about 2.6% at the two finest levels — so **2.6% is the number to
believe**, and the 4.7% headline is a design-grid figure. At 1.1 mm the
designed protocol also dips 0.21 °C/min below the critical rate: a design that
sits exactly on a constraint inherits that constraint's discretisation error.

The remedy is to design against a raised floor rather than the bare physical
one, and it works. Re-optimising with `--rate-margin-c-per-min 1.5`:

| spacing | design at the bare floor | design with a 1.5 °C/min margin |
|---:|---:|---:|
| 1.40 mm | 50.01 °C/min, feasible | 51.50 °C/min, feasible |
| 1.10 mm | 49.79 °C/min, **infeasible** | 51.27 °C/min, feasible |
| 0.90 mm | 50.79 °C/min, feasible | 52.32 °C/min, feasible |
| 0.75 mm | 50.79 °C/min, feasible | 52.31 °C/min, feasible |

The margined design is feasible at every resolution tested, and costs about
0.3 percentage points of stress reduction (4.37% rather than 4.71% on the design
grid) to buy it. That is the trade a practitioner should be making, and the
gradients make it explicit rather than accidental.

## A reproducibility finding in the source data

The paper states that nanoparticle SAR is held constant below −80 °C. Its own
properties then imply a 287.5 °C/min adiabatic ceiling, yet its tables report up
to 410 °C/min. Extending the measured fit instead gives 1013 W/gFe at −150 °C, a
421 °C/min ceiling, and reproduces the published table.

Both readings are implemented, and each is optimised separately. The outcome is
sharper than a robustness check:

![SAR robustness](docs/figures/sar_robustness.png)

Under the text-faithful clamp there is **no feasible protocol anywhere in the
design space** — not even maximum loading at full RF power reaches the critical
warming rate, missing it by 8.4 °C/min. One reading gives a designable problem;
the other says the experiment cannot work at all. That is independent evidence
about which reading the published simulations used.

## Cost, and differentiating in the other direction

The mechanics Tesseract exposes both a tangent-linear model and an **adjoint**.
The adjoint pulls the stress objective back through the entire equilibrium
trajectory in one reverse sweep, returning its sensitivity to every entry of the
temperature history at once; the thermal adjoint then carries that field to the
design variables in one more sweep. A whole Jacobian is five container
derivative calls — one mechanics reverse sweep and four thermal ones — against
twelve for the forward-mode equivalent.

Both modes are kept, because each is the other's only independent check. They
agree to **1.8×10⁻¹⁴** row by row, and the optimiser returns a design identical
to thirteen digits whichever it uses.

There is a third way to build the same Jacobian, which the released Tesseract
runtime cannot express at all. It needed a patch to Tesseract itself, described
in the next section.

| | forward | **batched** | adjoint |
|---|---:|---:|---:|
| container derivative calls per Jacobian | 12 | **2** | 5 |
| seconds per Jacobian (median of 7) | 5.84 | **1.36** | 2.33 |
| worst relative disagreement with forward | — | 8.5×10⁻¹⁵ | 1.8×10⁻¹⁴ |

The ratio matters more than the numbers. Forward mode costs one pass per design
variable, reverse mode one per Jacobian row, so the gap widens with every
variable added; at the per-voxel shell this method extends to, forward mode
simply stops being an option. The full run is 374 forward evaluations at 1.4 mm
(802 specimen cells) on one workstation.

## The Tesseract fork: batched JVP tangents

A JVP call carries exactly one tangent direction, because the runtime validates
every tangent against the exact shape of the input it belongs to. Getting `k`
directional derivatives is therefore `k` round trips — and each one repeats
everything the mechanics solver does *before* it looks at the direction:
assembling and factorising the equilibrium operator at every time step. That
work is the whole cost, and it is direction-independent.

[`Marc-Dvci/tesseract-core@batched-jvp`](https://github.com/Marc-Dvci/tesseract-core/tree/batched-jvp)
([diff](https://github.com/pasteurlabs/tesseract-core/compare/main...Marc-Dvci:tesseract-core:batched-jvp))
lets a tangent carry one extra leading axis of size `k`, with the results
carrying the same axis. Batching is detected against the *concrete* shape of the
input, so it can never be confused with an unbatched call, and the unbatched
path is unchanged. Mixing batched and unbatched tangents in one call is rejected
rather than broadcast. Tesseracts opt in — one that ignores the extra axis
returns a wrongly shaped result and the runtime rejects it — and Tesseracts
built on the JAX recipe get it for free, because `jax_jvp` maps the batch axis
with `vmap`. The thermal Tesseract here gained batched tangents with no change
at all.

This project installs the fork rather than the release, because `tesseract build`
copies the runtime from the installed package into both images.

Batched forward mode is **4.3× faster** than column-by-column forward mode, and
1.7× faster than the adjoint at six design variables, for the same Jacobian to
the round-off floor. It changes no result in this repository: the design, the
figures and every quoted number are unchanged, and `experiments/optimize_protocol`
reproduces the committed optimum bit for bit. What it changes is which regime is
affordable. Details, including what was deliberately left out, are in
[docs/fork.md](docs/fork.md).

![Every evaluation in the feasibility plane](docs/figures/feasibility_map.png)

## Reproduce

Requirements: Python 3.12, Docker, and `uv`. The install pulls the
[forked Tesseract runtime](https://github.com/Marc-Dvci/tesseract-core/tree/batched-jvp)
from GitHub, which is what puts batched tangents inside both images.

```bash
uv venv --python 3.12 .venv
uv pip install -e ".[dev,figures]"

tesseract build tesseracts/thermal
tesseract build tesseracts/mechanics

pytest -q                                        # kernels, correlations, results

python -m experiments.validate_forward_models    # both gates, vs published data
python -m experiments.validate_composed_gradient # tangent-linear vs finite differences
python -m experiments.optimize_protocol --sar-mode extrapolated \
    --output results/optimization_extrapolated.json
python -m experiments.optimize_protocol --sar-mode clamped \
    --output results/optimization_clamped.json
python -m experiments.evaluate_robustness
python -m experiments.check_objective_horizon
python -m experiments.grid_convergence
python -m experiments.compare_gradient_modes      # forward vs batched vs adjoint
python -m experiments.render_fields
python -m experiments.make_figures
```

Each experiment starts and tears down both services itself. Figures rebuild from
the committed JSON and `results/render_fields.npz` without Docker.

## Repository map

```text
src/vitrify/             cited properties, geometry, thermal kernel, Maxwell material
tesseracts/thermal/      JAX Tesseract API and image configuration
tesseracts/mechanics/    FEniCSx/PETSc Tesseract API, forward and tangent-linear solver
coupler/                 cross-container composition and the design parameterisation
experiments/             validation, optimisation, convergence, rendering, figures
tests/                   unit, correlation, derivative and committed-result tests
docs/                    validation record, technical writeup, the fork, figures
results/                 machine-readable outputs behind every number quoted here
```

## Scope and sources

The geometry is a volume- and dimension-matched analytic rat-heart phantom, not
a patient scan. All comparisons are simulation-to-simulation; no wet-lab claim
is made. The 3.2 MPa tensile strength is treated as a hard ceiling, which is
optimistic rather than conservative, since higher strain rates fail at lower
stress. Loading design is restricted to the surrounding solution because free
nanoparticle grading inside myocardium is not physiologically realizable. The
thermal-to-mechanical coupling is one-way; the genuine nonlinearity is inside
the thermal solve, where the nanoparticle source depends on the field it
creates.

Primary source: Joshi P, Rabin Y (2023), *Thermomechanical stress analyses of
nanowarming-assisted recovery from cryopreservation by vitrification in human
heart and rat heart models*, PLoS ONE 18(8): e0290063.

See [the validation record](docs/validation.md), [the technical
writeup](docs/writeup.md) and [the Tesseract fork](docs/fork.md).

Apache-2.0. See [LICENSE](LICENSE).

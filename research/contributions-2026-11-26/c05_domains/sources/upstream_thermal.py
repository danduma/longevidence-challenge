# Copyright 2026 Marc Davici
# SPDX-License-Identifier: Apache-2.0
"""Transient conduction with nanowarming, on a structured voxel grid.

Solves, per Joshi & Rabin (2023) Eqs. (1)-(3),

    C(T) dT/dt = div( k(T) grad T ) + q,        q = SAR(T) * C_Fe

by explicit finite volume with harmonic-mean face conductivities, subject to a
combined convective/radiative condition at the container's outer surface,

    -k dT/dn = U ( T_surface - T_environment(t) ).

Everything is differentiable with respect to the nanoparticle loading field and
the protocol parameters, which is the entire point: those are the design
variables.

The source term depends on temperature, so the heat source and the field it
produces are mutually determined. That is what stops per-particle superposition
from being valid and forces the nonlinear solve.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import partial

import jax
import jax.numpy as jnp
import numpy as np

from vitrify import geometry as geo
from vitrify.properties import CPA, POLYETHYLENE, nanowarming_source

__all__ = ["ThermalModel", "ProtocolBC", "SimulationResult"]

_AXES = (0, 1, 2)


@dataclass(frozen=True)
class ProtocolBC:
    """Boundary condition on the container's outer surface during rewarming.

    `u_w_per_m2k` of 0.0 is the adiabatic case used for Cases I and II of the
    published comparison. `t_env_c` may be a constant or a callable of time in
    seconds, allowing a ramped chamber temperature.
    """

    u_w_per_m2k: float = 0.0
    t_env_c: float = -150.0
    t_env_ramp_c_per_s: float = 0.0

    def environment_temperature(self, t_s):
        return self.t_env_c + self.t_env_ramp_c_per_s * t_s

    @property
    def is_adiabatic(self) -> bool:
        return self.u_w_per_m2k == 0.0


@dataclass
class SimulationResult:
    """Outputs of a rewarming simulation.

    `temperature_history` is (n_samples, nx, ny, nz) in degrees Celsius,
    sampled at `sample_times_s`. Reduction statistics are computed on the full
    time grid, not the samples, so they do not depend on sampling density.

    Two versions of the ice constraint are returned. The exact minimum is what
    feasibility is judged on. The softmin is a differentiable surrogate whose
    derivative averages over every cell and instant within a few degrees per
    minute of the worst one, instead of tracking whichever single cell happens
    to be the argmin. A hard min has a perfectly valid derivative almost
    everywhere, but it jumps when the argmin moves, and that discontinuity is
    what stalls a line search near an active constraint.
    """

    temperature_history: jnp.ndarray
    sample_times_s: jnp.ndarray
    min_warming_rate_c_per_min: jnp.ndarray
    smooth_min_warming_rate_c_per_min: jnp.ndarray
    final_min_temperature_c: jnp.ndarray
    final_max_temperature_c: jnp.ndarray
    mean_warming_rate_c_per_min: jnp.ndarray


class ThermalModel:
    """Precomputed grid operators for one phantom.

    Static structure (masks, neighbour connectivity, face areas) is built once
    in NumPy; only the fields that carry gradients stay in JAX.
    """

    def __init__(
        self, phantom: geo.Phantom, cpa: CPA, *, sar_low_temperature: str = "clamped"
    ):
        self.phantom = phantom
        self.cpa = cpa
        # Which reading of the unmeasured sub -80 C SAR behaviour to use; see
        # `properties.sar_emg308`. Absolute rates depend on it by up to 47%.
        self.sar_low_temperature = sar_low_temperature
        self.dx = float(phantom.spacing_m)

        labels = phantom.labels
        self.shape = labels.shape

        active = phantom.active_mask
        self.active = jnp.asarray(active, dtype=jnp.float64)

        # Cells that carry nanoparticles: CPA-filled regions, never the wall.
        self.is_cpa = jnp.asarray(phantom.specimen_mask, dtype=jnp.float64)
        self.is_wall = jnp.asarray(labels == geo.CONTAINER_WALL, dtype=jnp.float64)
        self.specimen_mask_np = phantom.specimen_mask

        # Per-direction neighbour masks. A face is internal when both cells are
        # active; it is an exterior face when the cell is active and its
        # neighbour is not.
        internal = []
        exterior_face_count = np.zeros(self.shape, dtype=np.float64)
        for axis in _AXES:
            for shift in (+1, -1):
                nbr_active = _shift(active, axis, shift, fill=False)
                internal.append(jnp.asarray(active & nbr_active, dtype=jnp.float64))
                exterior_face_count += (active & ~nbr_active).astype(np.float64)
        self.internal_face = internal
        self.exterior_face_count = jnp.asarray(exterior_face_count)

        # Subdomain masks for assigning the loading field.
        self.mask_chamber = jnp.asarray(labels == geo.CHAMBER, dtype=jnp.float64)
        self.mask_myocardium = jnp.asarray(labels == geo.MYOCARDIUM, dtype=jnp.float64)
        self.mask_surrounding = jnp.asarray(
            labels == geo.SURROUNDING_CPA, dtype=jnp.float64
        )

    # ------------------------------------------------------------------
    # Field helpers
    # ------------------------------------------------------------------

    def uniform_loading(
        self, chamber: float, myocardium: float, surrounding: float
    ) -> jnp.ndarray:
        """Build a piecewise-constant loading field, in mg Fe per mL."""
        return (
            chamber * self.mask_chamber
            + myocardium * self.mask_myocardium
            + surrounding * self.mask_surrounding
        )

    def _volumetric_heat_capacity(self, t_c: jnp.ndarray) -> jnp.ndarray:
        cpa = self.cpa.volumetric_heat_capacity(t_c)
        wall = POLYETHYLENE.volumetric_heat_capacity(t_c)
        return jnp.where(self.is_wall > 0.0, wall, cpa)

    def _conductivity(self, t_c: jnp.ndarray) -> jnp.ndarray:
        cpa = self.cpa.thermal_conductivity(t_c)
        wall = POLYETHYLENE.thermal_conductivity(t_c)
        return jnp.where(self.is_wall > 0.0, wall, cpa)

    # ------------------------------------------------------------------
    # Time stepping
    # ------------------------------------------------------------------

    def stable_dt(self, safety: float = 0.4) -> float:
        """Explicit stability limit, dt < rho c dx^2 / (6 k), with margin.

        Swept over the whole protocol temperature range and over both materials,
        because the binding combination (lowest heat capacity against highest
        conductivity) does not occur at the same temperature for each.
        """
        sweep = jnp.linspace(-180.0, 20.0, 400)
        rho_c = jnp.minimum(
            jnp.min(self.cpa.volumetric_heat_capacity(sweep)),
            jnp.min(POLYETHYLENE.volumetric_heat_capacity(sweep)),
        )
        k_max = jnp.maximum(
            jnp.max(self.cpa.thermal_conductivity(sweep)),
            jnp.max(POLYETHYLENE.thermal_conductivity(sweep)),
        )
        return float(safety * rho_c * self.dx**2 / (6.0 * k_max))

    def _rhs(
        self,
        t_field: jnp.ndarray,
        loading: jnp.ndarray,
        bc: ProtocolBC,
        t_s,
        duty_cycle,
    ):
        """dT/dt for every active cell."""
        k = self._conductivity(t_field)
        rho_c = self._volumetric_heat_capacity(t_field)

        flux_sum = jnp.zeros_like(t_field)
        idx = 0
        for axis in _AXES:
            for shift in (+1, -1):
                t_nbr = _jshift(t_field, axis, shift)
                k_nbr = _jshift(k, axis, shift)
                # Harmonic mean is the correct face conductivity for a jump in
                # material properties; the arithmetic mean leaks heat across the
                # CPA/polyethylene interface.
                k_face = 2.0 * k * k_nbr / (k + k_nbr + 1e-30)
                face_flux = k_face * (t_nbr - t_field) / self.dx
                flux_sum += self.internal_face[idx] * face_flux
                idx += 1

        if not bc.is_adiabatic:
            t_env = bc.environment_temperature(t_s)
            flux_sum += self.exterior_face_count * bc.u_w_per_m2k * (t_env - t_field)

        q = (
            nanowarming_source(
                t_field, loading, low_temperature=self.sar_low_temperature
            )
            * self.is_cpa
            * duty_cycle
        )
        return self.active * (flux_sum / self.dx + q) / rho_c

    @partial(jax.jit, static_argnums=(0, 4, 5, 6, 7, 8, 9, 10))
    def _run(
        self,
        loading: jnp.ndarray,
        t0: jnp.ndarray,
        power_schedule: jnp.ndarray,
        dt: float,
        n_steps: int,
        n_chunks: int,
        chunk_size: int,
        bc: ProtocolBC,
        n_specimen: int,
        rate_smoothing: float,
    ):
        """Two-level scan: the inner one is rematerialised, the outer one stores.

        Storing every step of a 3D field would cost hundreds of megabytes and
        most of it would be discarded at output time. The inner scan is wrapped
        in `jax.checkpoint` so reverse mode recomputes a chunk rather than
        keeping it, and only chunk-final states are stacked.
        """
        scale = self.specimen_mask_np.size / n_specimen

        # Softmin is accumulated in streaming log-sum-exp form: carry the
        # running minimum and the exponential sum rescaled to it, so nothing
        # overflows across millions of (cell, instant) pairs.
        def step(carry, i):
            t_field, rate_min, soft_shift, soft_sum, rate_sum = carry
            valid = i < n_steps
            schedule_index = jnp.minimum(i, n_steps - 1)
            dtdt = self._rhs(
                t_field,
                loading,
                bc,
                schedule_index * dt,
                power_schedule[schedule_index],
            )
            t_next = jnp.where(valid, t_field + dt * dtdt, t_field)

            # Worst instantaneous warming rate anywhere in the specimen: this is
            # the ice constraint. It is a min over both space and time, so it
            # cannot be recovered from sampled output and has to be carried.
            rate = jnp.where(self.specimen_mask_np, dtdt, jnp.inf)
            step_min = jnp.min(rate)
            rate_min = jnp.where(valid, jnp.minimum(rate_min, step_min), rate_min)

            shift = jnp.minimum(soft_shift, step_min)
            step_sum = jnp.sum(
                jnp.where(
                    self.specimen_mask_np,
                    jnp.exp(-(dtdt - shift) / rate_smoothing),
                    0.0,
                )
            )
            soft_sum = jnp.where(
                valid,
                soft_sum * jnp.exp(-(soft_shift - shift) / rate_smoothing)
                + step_sum,
                soft_sum,
            )
            soft_shift = jnp.where(valid, shift, soft_shift)

            rate_sum = rate_sum + jnp.where(
                valid,
                jnp.mean(jnp.where(self.specimen_mask_np, dtdt, 0.0)),
                0.0,
            )
            return (t_next, rate_min, soft_shift, soft_sum, rate_sum), None

        @jax.checkpoint
        def chunk(carry, c):
            steps = c * chunk_size + jnp.arange(chunk_size)
            carry, _ = jax.lax.scan(step, carry, steps)
            return carry, carry[0]

        (t_final, rate_min, soft_shift, soft_sum, rate_sum), history = jax.lax.scan(
            chunk,
            (t0, jnp.inf, jnp.inf, 0.0, 0.0),
            jnp.arange(n_chunks),
        )
        soft_min = soft_shift - rate_smoothing * jnp.log(soft_sum)
        return t_final, rate_min, soft_min, rate_sum / n_steps * scale, history

    def simulate_rewarming(
        self,
        loading: jnp.ndarray,
        *,
        t_start_c: float = -150.0,
        duration_s: float = 40.0,
        bc: ProtocolBC | None = None,
        dt: float | None = None,
        n_samples: int = 40,
        duty_cycle: jnp.ndarray | None = None,
        rate_smoothing_c_per_min: float = 1.0,
    ) -> SimulationResult:
        """Rewarm from uniform cryogenic storage.

        Cases I-IV of the published comparison share an identical cooling and
        storage history, so the rewarming phase can be simulated on its own from
        a uniform initial temperature - which is what makes them comparable.
        """
        bc = bc or ProtocolBC()
        dt = dt or self.stable_dt()
        n_steps = int(np.ceil(duration_s / dt))
        chunk_size = max(1, int(np.ceil(n_steps / n_samples)))
        n_chunks = int(np.ceil(n_steps / chunk_size))

        if duty_cycle is None:
            duty_cycle = jnp.ones((1,), dtype=jnp.float64)
        else:
            duty_cycle = jnp.asarray(duty_cycle, dtype=jnp.float64).reshape(-1)
            if duty_cycle.size == 0:
                raise ValueError("duty_cycle must contain at least one control point")
        control_x = jnp.linspace(0.0, 1.0, duty_cycle.size)
        step_x = jnp.linspace(0.0, 1.0, n_steps)
        power_schedule = jnp.interp(step_x, control_x, duty_cycle)

        t0 = jnp.full(self.shape, t_start_c, dtype=jnp.float64) * self.active
        n_specimen = int(np.count_nonzero(self.specimen_mask_np))
        t_final, rate_min, soft_min, rate_mean, history = self._run(
            loading,
            t0,
            power_schedule,
            dt,
            n_steps,
            n_chunks,
            chunk_size,
            bc,
            n_specimen,
            float(rate_smoothing_c_per_min) / 60.0,
        )

        spec = jnp.asarray(self.specimen_mask_np)
        times = jnp.minimum(
            (jnp.arange(n_chunks) + 1) * chunk_size,
            n_steps,
        ) * dt
        return SimulationResult(
            temperature_history=history,
            sample_times_s=times,
            min_warming_rate_c_per_min=rate_min * 60.0,
            smooth_min_warming_rate_c_per_min=soft_min * 60.0,
            mean_warming_rate_c_per_min=rate_mean * 60.0,
            final_min_temperature_c=jnp.min(jnp.where(spec, t_final, jnp.inf)),
            final_max_temperature_c=jnp.max(jnp.where(spec, t_final, -jnp.inf)),
        )


# ----------------------------------------------------------------------
# Shift helpers
# ----------------------------------------------------------------------


def _shift(arr: np.ndarray, axis: int, shift: int, fill=False) -> np.ndarray:
    out = np.full_like(arr, fill)
    src = [slice(None)] * arr.ndim
    dst = [slice(None)] * arr.ndim
    if shift > 0:
        src[axis] = slice(shift, None)
        dst[axis] = slice(None, -shift)
    else:
        src[axis] = slice(None, shift)
        dst[axis] = slice(-shift, None)
    out[tuple(dst)] = arr[tuple(src)]
    return out


def _jshift(arr: jnp.ndarray, axis: int, shift: int) -> jnp.ndarray:
    """Shift with edge padding; exterior faces are masked out anyway."""
    return jnp.roll(arr, -shift, axis=axis)

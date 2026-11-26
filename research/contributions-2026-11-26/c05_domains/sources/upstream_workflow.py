# Copyright 2026 Marc Davici
# SPDX-License-Identifier: Apache-2.0
"""Compose the JAX thermal and FEniCSx mechanics Tesseracts.

The two containers are joined into one differentiable function of the design
vector. Forward mode is used across the boundary because the design space is
six-dimensional while the field crossing the boundary has thousands of
components: a directional derivative costs one thermal tangent and one
mechanics tangent regardless of how many temperature voxels there are.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import numpy.typing as npt
from tesseract_core import Tesseract

from vitrify.geometry import CHAMBER, MYOCARDIUM, SURROUNDING_CPA, Phantom

SarMode = Literal["clamped", "extrapolated"]

#: Physiological loading that cannot be designed: nanoparticles stay in the
#: vasculature, so chamber and myocardium concentrations are fixed.
CHAMBER_LOADING_MG_PER_ML = 10.0
MYOCARDIUM_LOADING_MG_PER_ML = 1.47


@dataclass(frozen=True)
class WorkflowResult:
    stress_objective_mpa: float
    peak_tension_mpa: float
    peak_compression_mpa: float
    peak_time_s: float
    min_warming_rate_c_per_min: float
    smooth_min_warming_rate_c_per_min: float
    mean_warming_rate_c_per_min: float
    final_min_temperature_c: float
    final_max_temperature_c: float
    mean_surrounding_loading_mg_per_ml: float
    loading_field_mg_per_ml: np.ndarray
    duty_cycle: np.ndarray
    peak_max_principal_pa: np.ndarray
    stress_trace_pa: np.ndarray
    stress_history_pa: np.ndarray = field(default_factory=lambda: np.zeros(0))
    temperature_history_c: np.ndarray = field(default_factory=lambda: np.zeros(0))
    sample_times_s: np.ndarray = field(default_factory=lambda: np.zeros(0))

    def constraint_margins(self) -> dict[str, float]:
        """Signed feasibility margins; every one must be non-negative."""
        return {
            "ice_rate": self.min_warming_rate_c_per_min - MIN_WARMING_RATE_C_PER_MIN,
            "cold_spot": self.final_min_temperature_c - FINAL_MIN_TEMPERATURE_C,
            "hot_spot": FINAL_MAX_TEMPERATURE_C - self.final_max_temperature_c,
        }

    @property
    def is_feasible(self) -> bool:
        return all(value >= 0.0 for value in self.constraint_margins().values())


@dataclass(frozen=True)
class WorkflowTangent:
    stress_objective_mpa: float
    min_warming_rate_c_per_min: float
    smooth_min_warming_rate_c_per_min: float
    final_min_temperature_c: float
    final_max_temperature_c: float
    mean_surrounding_loading_mg_per_ml: float


#: Critical warming rate for VS55 (Joshi & Rabin 2023) and the published
#: final-temperature window, 3 C clear of the melting point.
MIN_WARMING_RATE_C_PER_MIN = 50.0
FINAL_MIN_TEMPERATURE_C = -35.0
FINAL_MAX_TEMPERATURE_C = 0.0


class DesignParameterization:
    """Three radial loading bands and three RF duty-cycle control points."""

    n_loading_parameters = 3
    n_duty_parameters = 3
    size = n_loading_parameters + n_duty_parameters
    lower_bounds = np.array([0.0, 0.0, 0.0, 0.2, 0.2, 0.2])
    upper_bounds = np.array([10.0, 10.0, 10.0, 1.0, 1.0, 1.0])

    def __init__(self, phantom: Phantom):
        self.phantom = phantom
        shape = np.asarray(phantom.shape)
        grid = np.meshgrid(
            *(np.arange(n, dtype=np.float64) + 0.5 for n in shape),
            indexing="ij",
        )
        centre = shape / 2.0
        # The heart's measured bounding dimensions define an ellipsoidal radial
        # coordinate. The three overlapping bands form a smooth printable shell.
        semi_axes = np.asarray(phantom.spec.heart_dims_mm) / (
            2.0 * phantom.spacing_m * 1.0e3
        )
        radius = np.sqrt(
            sum(
                ((grid[axis] - centre[axis]) / semi_axes[axis]) ** 2
                for axis in range(3)
            )
        )
        centres = np.array([1.05, 1.45, 2.0])
        widths = np.array([0.32, 0.42, 0.65])
        raw = np.stack(
            [
                np.exp(-0.5 * ((radius - centre_r) / width) ** 2)
                for centre_r, width in zip(centres, widths, strict=True)
            ]
        )
        raw_sum = raw.sum(axis=0, keepdims=True)
        self.loading_basis = raw / np.maximum(raw_sum, 1.0e-15)
        self.surrounding = phantom.labels == SURROUNDING_CPA
        self.loading_basis *= self.surrounding[None]
        self.band_mass_weights = self.loading_basis[:, self.surrounding].mean(axis=1)

    def fields(self, parameters: npt.ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        values = np.asarray(parameters, dtype=np.float64)
        if values.shape != (self.size,):
            raise ValueError(f"expected {self.size} parameters, got {values.shape}")
        concentrations = values[: self.n_loading_parameters]
        duty_cycle = values[self.n_loading_parameters :]
        loading = np.tensordot(concentrations, self.loading_basis, axes=1)
        loading = np.where(
            self.phantom.labels == CHAMBER, CHAMBER_LOADING_MG_PER_ML, loading
        )
        loading = np.where(
            self.phantom.labels == MYOCARDIUM, MYOCARDIUM_LOADING_MG_PER_ML, loading
        )
        return loading, duty_cycle

    def field_tangents(self, direction: npt.ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        tangent = np.asarray(direction, dtype=np.float64)
        loading = np.tensordot(
            tangent[: self.n_loading_parameters], self.loading_basis, axes=1
        )
        duty = tangent[self.n_loading_parameters :]
        return loading, duty

    def mean_surrounding_loading(self, parameters: npt.ArrayLike) -> float:
        values = np.asarray(parameters, dtype=np.float64)
        return float(self.band_mass_weights @ values[: self.n_loading_parameters])

    def mean_surrounding_loading_jvp(self, direction: npt.ArrayLike) -> float:
        tangent = np.asarray(direction, dtype=np.float64)
        return float(self.band_mass_weights @ tangent[: self.n_loading_parameters])

    def pull_back(
        self, loading_cotangent: np.ndarray, duty_cotangent: np.ndarray
    ) -> np.ndarray:
        """Map a gradient on the fields back onto the six design parameters.

        The chamber and myocardium concentrations are constants rather than
        design variables, and the radial basis is already zero outside the
        surrounding solution, so the contraction below silently ignores them.
        """
        loading = np.tensordot(
            self.loading_basis,
            np.asarray(loading_cotangent, dtype=np.float64),
            axes=([1, 2, 3], [0, 1, 2]),
        )
        return np.concatenate([loading, np.asarray(duty_cotangent, dtype=np.float64)])

    def uniform(self, surrounding_mg_per_ml: float, duty: float = 1.0) -> np.ndarray:
        """The published piecewise-constant protocols, in design coordinates."""
        return np.array([surrounding_mg_per_ml] * 3 + [duty] * 3)


class ComposedWorkflow:
    """A single differentiable function spanning two served Tesseracts."""

    #: Jacobian rows, in the order the optimizer indexes them: the objective,
    #: then the three constraint quantities it linearises.
    tangent_fields = (
        "stress_objective_mpa",
        "smooth_min_warming_rate_c_per_min",
        "final_min_temperature_c",
        "final_max_temperature_c",
    )
    #: Thermal outputs whose gradients come from a single reverse sweep each.
    _scalar_rows = (
        "smooth_min_warming_rate_c_per_min",
        "final_min_temperature_c",
        "final_max_temperature_c",
    )

    def __init__(
        self,
        phantom: Phantom,
        *,
        duration_s: float = 70.0,
        n_samples: int = 40,
        mechanics_window_s: float = 12.0,
        rate_smoothing_c_per_min: float = 1.0,
        stream_logs: bool = False,
    ):
        self.phantom = phantom
        self.grid_spacing_mm = float(phantom.spacing_m * 1.0e3)
        self.parameterization = DesignParameterization(phantom)
        self.duration_s = duration_s
        self.n_samples = n_samples
        self.mechanics_window_s = mechanics_window_s
        self.rate_smoothing_c_per_min = rate_smoothing_c_per_min
        self.thermal = Tesseract.from_image(
            "vitrify-thermal:latest", stream_logs=stream_logs, timeout=1800.0
        )
        self.mechanics = Tesseract.from_image(
            "vitrify-mechanics:latest", stream_logs=stream_logs, timeout=1800.0
        )
        self.n_forward_calls = 0
        self.n_tangent_calls = 0
        self.n_adjoint_calls = 0
        #: Derivative endpoints actually invoked on a container, which is the
        #: quantity that costs money: a batched call is one of these, not k.
        self.n_container_derivative_calls = 0
        self._last_parameters: np.ndarray | None = None
        self._last_sar_mode: str | None = None
        self._thermal_inputs: dict | None = None
        self._thermal_outputs: dict | None = None
        self._mechanics_inputs: dict | None = None

    def __enter__(self):
        self.thermal.serve()
        self.mechanics.serve()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.teardown()

    def teardown(self):
        self.mechanics.teardown()
        self.thermal.teardown()

    def _thermal_input(
        self, loading: np.ndarray, duty_cycle: np.ndarray, sar_mode: SarMode
    ) -> dict:
        return {
            "loading_mg_per_ml": loading,
            "duty_cycle": duty_cycle,
            "sar_mode": sar_mode,
            "duration_s": self.duration_s,
            "n_samples": self.n_samples,
            "u_w_per_m2k": 0.0,
            "environment_temperature_c": -150.0,
            "grid_spacing_mm": self.grid_spacing_mm,
            "rate_smoothing_c_per_min": self.rate_smoothing_c_per_min,
        }

    def forward(
        self,
        parameters: npt.ArrayLike,
        *,
        sar_mode: SarMode = "extrapolated",
        keep_fields: bool = True,
    ) -> WorkflowResult:
        values = np.asarray(parameters, dtype=np.float64)
        loading, duty = self.parameterization.fields(values)
        thermal_inputs = self._thermal_input(loading, duty, sar_mode)
        thermal_outputs = self.thermal.apply(thermal_inputs)
        times = np.asarray(thermal_outputs["sample_times_s"])
        mechanics_indices = np.flatnonzero(times <= self.mechanics_window_s)
        if len(mechanics_indices) < 2:
            raise RuntimeError("mechanics window contains fewer than two samples")
        mechanics_inputs = {
            "temperature_history_c": np.asarray(
                thermal_outputs["temperature_history_c"]
            )[mechanics_indices],
            "sample_times_s": times[mechanics_indices],
            "grid_spacing_mm": self.grid_spacing_mm,
        }
        mechanics_outputs = self.mechanics.apply(mechanics_inputs)
        self.n_forward_calls += 1

        self._last_parameters = values.copy()
        self._last_sar_mode = sar_mode
        self._thermal_inputs = thermal_inputs
        self._thermal_outputs = thermal_outputs
        self._mechanics_inputs = mechanics_inputs

        empty = np.zeros(0)
        return WorkflowResult(
            stress_objective_mpa=float(mechanics_outputs["stress_objective_pa"]) / 1e6,
            peak_tension_mpa=float(mechanics_outputs["peak_tension_pa"]) / 1e6,
            peak_compression_mpa=float(mechanics_outputs["peak_compression_pa"]) / 1e6,
            peak_time_s=float(mechanics_outputs["peak_time_s"]),
            min_warming_rate_c_per_min=float(
                thermal_outputs["min_warming_rate_c_per_min"]
            ),
            smooth_min_warming_rate_c_per_min=float(
                thermal_outputs["smooth_min_warming_rate_c_per_min"]
            ),
            mean_warming_rate_c_per_min=float(
                thermal_outputs["mean_warming_rate_c_per_min"]
            ),
            final_min_temperature_c=float(thermal_outputs["final_min_temperature_c"]),
            final_max_temperature_c=float(thermal_outputs["final_max_temperature_c"]),
            mean_surrounding_loading_mg_per_ml=(
                self.parameterization.mean_surrounding_loading(values)
            ),
            loading_field_mg_per_ml=loading,
            duty_cycle=duty,
            peak_max_principal_pa=np.asarray(
                mechanics_outputs["peak_max_principal_pa"]
            ),
            stress_trace_pa=np.asarray(mechanics_outputs["stress_trace_pa"]),
            stress_history_pa=(
                np.asarray(mechanics_outputs["stress_history_pa"])
                if keep_fields
                else empty
            ),
            temperature_history_c=(
                np.asarray(thermal_outputs["temperature_history_c"])
                if keep_fields
                else empty
            ),
            sample_times_s=times if keep_fields else empty,
        )

    def jvp(
        self,
        parameters: npt.ArrayLike,
        direction: npt.ArrayLike,
        *,
        sar_mode: SarMode = "extrapolated",
    ) -> WorkflowTangent:
        values = np.asarray(parameters, dtype=np.float64)
        if (
            self._last_parameters is None
            or self._last_sar_mode != sar_mode
            or not np.array_equal(values, self._last_parameters)
        ):
            self.forward(values, sar_mode=sar_mode)
        loading_tangent, duty_tangent = self.parameterization.field_tangents(direction)
        thermal_tangent = self.thermal.jacobian_vector_product(
            self._thermal_inputs,
            ["loading_mg_per_ml", "duty_cycle"],
            [
                "temperature_history_c",
                "min_warming_rate_c_per_min",
                "smooth_min_warming_rate_c_per_min",
                "final_min_temperature_c",
                "final_max_temperature_c",
            ],
            {"loading_mg_per_ml": loading_tangent, "duty_cycle": duty_tangent},
        )
        mechanics_count = len(self._mechanics_inputs["sample_times_s"])
        mechanics_tangent = self.mechanics.jacobian_vector_product(
            self._mechanics_inputs,
            ["temperature_history_c"],
            ["stress_objective_pa"],
            {
                "temperature_history_c": np.asarray(
                    thermal_tangent["temperature_history_c"]
                )[:mechanics_count]
            },
        )
        self.n_tangent_calls += 1
        self.n_container_derivative_calls += 2
        return WorkflowTangent(
            stress_objective_mpa=float(mechanics_tangent["stress_objective_pa"]) / 1e6,
            min_warming_rate_c_per_min=float(
                thermal_tangent["min_warming_rate_c_per_min"]
            ),
            smooth_min_warming_rate_c_per_min=float(
                thermal_tangent["smooth_min_warming_rate_c_per_min"]
            ),
            final_min_temperature_c=float(thermal_tangent["final_min_temperature_c"]),
            final_max_temperature_c=float(thermal_tangent["final_max_temperature_c"]),
            mean_surrounding_loading_mg_per_ml=(
                self.parameterization.mean_surrounding_loading_jvp(direction)
            ),
        )

    def batched_jvp(
        self,
        parameters: npt.ArrayLike,
        directions: npt.ArrayLike,
        *,
        sar_mode: SarMode = "extrapolated",
    ) -> list[WorkflowTangent]:
        """All ``k`` directional derivatives, in one call per container.

        The unbatched :meth:`jvp` sends one direction per call, so a Jacobian
        costs ``k`` round trips through each container and, in the mechanics,
        ``k`` complete equilibrium trajectories. With the ``batched-jvp`` fork
        of tesseract-core the tangent may carry a leading axis, so this is two
        calls whatever ``k`` is: the thermal solver vmaps the batch, and the
        mechanics solver assembles and factorises once per step and then
        back-substitutes ``k`` times.
        """
        values = np.asarray(parameters, dtype=np.float64)
        stack = np.atleast_2d(np.asarray(directions, dtype=np.float64))
        if (
            self._last_parameters is None
            or self._last_sar_mode != sar_mode
            or not np.array_equal(values, self._last_parameters)
        ):
            self.forward(values, sar_mode=sar_mode)

        field_tangents = [
            self.parameterization.field_tangents(direction) for direction in stack
        ]
        loading_tangents = np.stack([loading for loading, _ in field_tangents])
        duty_tangents = np.stack([duty for _, duty in field_tangents])

        thermal_tangent = self.thermal.jacobian_vector_product(
            self._thermal_inputs,
            ["loading_mg_per_ml", "duty_cycle"],
            [
                "temperature_history_c",
                "min_warming_rate_c_per_min",
                "smooth_min_warming_rate_c_per_min",
                "final_min_temperature_c",
                "final_max_temperature_c",
            ],
            {"loading_mg_per_ml": loading_tangents, "duty_cycle": duty_tangents},
        )
        mechanics_count = len(self._mechanics_inputs["sample_times_s"])
        mechanics_tangent = self.mechanics.jacobian_vector_product(
            self._mechanics_inputs,
            ["temperature_history_c"],
            ["stress_objective_pa"],
            {
                "temperature_history_c": np.asarray(
                    thermal_tangent["temperature_history_c"]
                )[:, :mechanics_count]
            },
        )
        self.n_tangent_calls += 1
        self.n_container_derivative_calls += 2

        stress = np.asarray(mechanics_tangent["stress_objective_pa"])
        return [
            WorkflowTangent(
                stress_objective_mpa=float(stress[index]) / 1e6,
                min_warming_rate_c_per_min=float(
                    np.asarray(thermal_tangent["min_warming_rate_c_per_min"])[index]
                ),
                smooth_min_warming_rate_c_per_min=float(
                    np.asarray(
                        thermal_tangent["smooth_min_warming_rate_c_per_min"]
                    )[index]
                ),
                final_min_temperature_c=float(
                    np.asarray(thermal_tangent["final_min_temperature_c"])[index]
                ),
                final_max_temperature_c=float(
                    np.asarray(thermal_tangent["final_max_temperature_c"])[index]
                ),
                mean_surrounding_loading_mg_per_ml=(
                    self.parameterization.mean_surrounding_loading_jvp(stack[index])
                ),
            )
            for index in range(len(stack))
        ]

    def _ensure(self, values: np.ndarray, sar_mode: SarMode) -> None:
        if (
            self._last_parameters is None
            or self._last_sar_mode != sar_mode
            or not np.array_equal(values, self._last_parameters)
        ):
            self.forward(values, sar_mode=sar_mode, keep_fields=False)

    def _thermal_pull_back(self, outputs: dict) -> np.ndarray:
        return self.parameterization.pull_back(
            np.asarray(outputs["loading_mg_per_ml"]),
            np.asarray(outputs["duty_cycle"]),
        )

    def jacobian(
        self,
        parameters: npt.ArrayLike,
        *,
        sar_mode: SarMode = "extrapolated",
        mode: Literal["adjoint", "forward", "batched"] = "adjoint",
    ) -> np.ndarray:
        """Full ``(len(tangent_fields), 6)`` Jacobian of the composed model.

        Adjoint mode sweeps backwards once per row: the mechanics adjoint turns
        the stress objective into a sensitivity field over the whole
        temperature history, the thermal adjoint carries that field back to the
        design variables in one more pass, and the three constraint rows are
        thermal-only reverse sweeps. Four passes, whatever the design dimension.

        Forward mode is kept because it is the independent check on all of
        that, and because for a design space smaller than the number of rows it
        would be the cheaper choice.

        Batched mode is forward mode with the six directions sent together, so
        the mechanics pays for one set of factorisations instead of six. It
        needs the ``batched-jvp`` fork of tesseract-core in both images.
        """
        values = np.asarray(parameters, dtype=np.float64)
        if mode in ("forward", "batched"):
            identity = np.eye(self.parameterization.size)
            if mode == "batched":
                tangents = self.batched_jvp(values, identity, sar_mode=sar_mode)
            else:
                tangents = [
                    self.jvp(values, direction, sar_mode=sar_mode)
                    for direction in identity
                ]
            columns = [
                [getattr(tangent, name) for name in self.tangent_fields]
                for tangent in tangents
            ]
            return np.asarray(columns).T

        self._ensure(values, sar_mode)
        n_samples = len(np.asarray(self._thermal_outputs["sample_times_s"]))
        mechanics_count = len(self._mechanics_inputs["sample_times_s"])

        # The stress objective, back through both containers.
        mechanics_cotangent = self.mechanics.vector_jacobian_product(
            self._mechanics_inputs,
            ["temperature_history_c"],
            ["stress_objective_pa"],
            {"stress_objective_pa": 1.0},
        )
        field = np.asarray(mechanics_cotangent["temperature_history_c"])
        # The mechanics only sees the early window; the rest of the history
        # carries no stress sensitivity.
        padded = np.zeros((n_samples, *field.shape[1:]))
        padded[:mechanics_count] = field
        objective = self.thermal.vector_jacobian_product(
            self._thermal_inputs,
            ["loading_mg_per_ml", "duty_cycle"],
            ["temperature_history_c"],
            {"temperature_history_c": padded},
        )
        rows = [self._thermal_pull_back(objective) / 1e6]

        # The constraints never leave the thermal container.
        for name in self._scalar_rows:
            outputs = self.thermal.vector_jacobian_product(
                self._thermal_inputs,
                ["loading_mg_per_ml", "duty_cycle"],
                [name],
                {name: 1.0},
            )
            rows.append(self._thermal_pull_back(outputs))
        self.n_adjoint_calls += 1
        # One mechanics reverse sweep, one thermal sweep for the objective, and
        # one thermal sweep per constraint row.
        self.n_container_derivative_calls += 2 + len(self._scalar_rows)
        return np.asarray(rows)

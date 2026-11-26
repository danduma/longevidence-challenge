# Copyright 2026 Marc Davici
# SPDX-License-Identifier: Apache-2.0
"""Voxelised heart-in-container phantoms.

The published models are built on an MRI reconstruction of a human heart from
the Visible Heart Laboratory, which is not redistributed with the paper. This
module builds an analytic phantom instead: a prolate ellipsoidal myocardial
shell containing two chamber cavities, immersed in CPA inside a pillow-shaped
cryobag, with outer dimensions *and* subdomain volumes matched to the published
values (Joshi & Rabin 2023, Table 1).

Both are matched because either alone is easy to hit and misleading on its own:
a box with the right bounding dimensions overstates the cryobag volume by 50%,
and a shape with the right volume can have quite the wrong surface area, which
is what sets the convective boundary condition.

The substitution is deliberate and is why this project validates against the
*ordering and scale* of the published stress results rather than claiming to
reproduce absolute values voxel for voxel. Point-for-point agreement from a
different geometry would be a red flag, not a success.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "Subdomain",
    "Phantom",
    "PhantomSpec",
    "build_phantom",
    "rat_heart_cryobag",
    "RAT_HEART_CRYOBAG",
]

EXTERIOR = -1
CONTAINER_WALL = 0
SURROUNDING_CPA = 1
MYOCARDIUM = 2
CHAMBER = 3


class Subdomain:
    """Named label values, for readability at call sites."""

    EXTERIOR = EXTERIOR
    CONTAINER_WALL = CONTAINER_WALL
    SURROUNDING_CPA = SURROUNDING_CPA
    MYOCARDIUM = MYOCARDIUM
    CHAMBER = CHAMBER


@dataclass(frozen=True)
class PhantomSpec:
    """Target dimensions, in millimetres and millilitres."""

    name: str
    container_dims_mm: tuple[float, float, float]
    container_volume_ml: float
    heart_dims_mm: tuple[float, float, float]
    heart_volume_ml: float
    wall_thickness_mm: float = 1.0
    # Fraction of total heart volume that is chamber cavity rather than muscle.
    # Micro-CT of the rat heart showed 100% CPA loading in the ventricles; the
    # split itself is a modelling choice, reported by `Phantom.report()`.
    chamber_fraction: float = 0.35


# Joshi & Rabin 2023, Table 1.
RAT_HEART_CRYOBAG = PhantomSpec(
    name="rat heart in cryobag",
    container_dims_mm=(23.4, 17.5, 12.5),
    container_volume_ml=2.2,
    heart_dims_mm=(17.5, 13.0, 11.0),
    heart_volume_ml=0.69,
)


@dataclass
class Phantom:
    """A labelled voxel grid plus its metrics."""

    labels: np.ndarray  # int8, shape (nx, ny, nz)
    spacing_m: float
    spec: PhantomSpec
    pillow_exponent: float

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.labels.shape

    @property
    def cell_volume_m3(self) -> float:
        return self.spacing_m**3

    def volume_ml(self, *labels: int) -> float:
        n = int(np.count_nonzero(self.mask(*labels)))
        return n * (self.spacing_m * 1e3) ** 3 / 1e3  # mm^3 -> mL

    def mask(self, *labels: int) -> np.ndarray:
        m = np.zeros(self.labels.shape, dtype=bool)
        for lab in labels:
            m |= self.labels == lab
        return m

    @property
    def specimen_mask(self) -> np.ndarray:
        """Everything that must be vitrified: CPA and tissue, not the wall."""
        return self.mask(SURROUNDING_CPA, MYOCARDIUM, CHAMBER)

    @property
    def active_mask(self) -> np.ndarray:
        """Cells the heat solver integrates: specimen plus container wall."""
        return self.mask(CONTAINER_WALL, SURROUNDING_CPA, MYOCARDIUM, CHAMBER)

    def report(self) -> dict[str, float]:
        return {
            "cells_total": float(self.labels.size),
            "cells_active": float(np.count_nonzero(self.active_mask)),
            "spacing_mm": self.spacing_m * 1e3,
            "pillow_exponent": self.pillow_exponent,
            "container_contents_ml": self.volume_ml(
                SURROUNDING_CPA, MYOCARDIUM, CHAMBER
            ),
            "heart_ml": self.volume_ml(MYOCARDIUM, CHAMBER),
            "myocardium_ml": self.volume_ml(MYOCARDIUM),
            "chamber_ml": self.volume_ml(CHAMBER),
            "wall_ml": self.volume_ml(CONTAINER_WALL),
        }


def _superellipsoid(
    grid: tuple[np.ndarray, np.ndarray, np.ndarray],
    centre_mm: tuple[float, float, float],
    semi_axes_mm: tuple[float, float, float],
    exponent: float,
) -> np.ndarray:
    """|dx/a|^p + |dy/b|^p + |dz/c|^p <= 1.

    p = 2 is an ellipsoid; p -> infinity approaches the bounding box. A cryobag
    lies between the two, which is exactly the free parameter we need.
    """
    x, y, z = grid
    cx, cy, cz = centre_mm
    ax, ay, az = semi_axes_mm
    return (
        np.abs((x - cx) / ax) ** exponent
        + np.abs((y - cy) / ay) ** exponent
        + np.abs((z - cz) / az) ** exponent
    ) <= 1.0


def _solve_exponent(
    grid, centre, semi_axes, cell_mm3: float, target_ml: float
) -> float:
    """Bisect the pillow exponent so the voxelised interior hits `target_ml`.

    Solved on the voxelisation rather than analytically, so the reported volume
    is the volume the solver actually integrates over.
    """
    target_mm3 = target_ml * 1e3

    def volume(p: float) -> float:
        count = np.count_nonzero(_superellipsoid(grid, centre, semi_axes, p))
        return float(count) * cell_mm3

    lo, hi = 1.2, 12.0
    if volume(hi) < target_mm3:
        return hi  # box-shaped and still too small; nothing more to give
    if volume(lo) > target_mm3:
        return lo
    for _ in range(48):
        mid = 0.5 * (lo + hi)
        if volume(mid) < target_mm3:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def build_phantom(spec: PhantomSpec, spacing_mm: float = 0.4) -> Phantom:
    """Voxelise `spec` at the given spacing.

    Two shape parameters are solved rather than assumed: the cryobag's pillow
    exponent (to match container volume at fixed bounding dimensions) and a
    scalar shrink on the myocardial ellipsoid (to match heart volume, since a
    heart does not fill its bounding ellipsoid).
    """
    lx, ly, lz = spec.container_dims_mm
    nx, ny, nz = (int(round(v / spacing_mm)) for v in (lx, ly, lz))

    xs = (np.arange(nx) + 0.5) * spacing_mm
    ys = (np.arange(ny) + 0.5) * spacing_mm
    zs = (np.arange(nz) + 0.5) * spacing_mm
    grid = np.meshgrid(xs, ys, zs, indexing="ij")
    cell_mm3 = spacing_mm**3
    centre = (lx / 2.0, ly / 2.0, lz / 2.0)

    w = spec.wall_thickness_mm
    inner_semi = (lx / 2.0 - w, ly / 2.0 - w, lz / 2.0 - w)
    outer_semi = (lx / 2.0, ly / 2.0, lz / 2.0)

    p = _solve_exponent(grid, centre, inner_semi, cell_mm3, spec.container_volume_ml)

    interior = _superellipsoid(grid, centre, inner_semi, p)
    outer = _superellipsoid(grid, centre, outer_semi, p)

    labels = np.full((nx, ny, nz), EXTERIOR, dtype=np.int8)
    labels[outer] = CONTAINER_WALL
    labels[interior] = SURROUNDING_CPA

    # Myocardium: ellipsoid shrunk to match the published heart volume.
    hx, hy, hz = spec.heart_dims_mm
    bounding_semi = (hx / 2.0, hy / 2.0, hz / 2.0)
    target_mm3 = spec.heart_volume_ml * 1e3
    bounding_mm3 = (
        4.0 / 3.0 * np.pi * bounding_semi[0] * bounding_semi[1] * bounding_semi[2]
    )
    scale = float(np.cbrt(target_mm3 / bounding_mm3))
    heart_semi = tuple(s * scale for s in bounding_semi)

    heart = _superellipsoid(grid, centre, heart_semi, 2.0) & interior
    labels[heart] = MYOCARDIUM

    # Two chamber cavities offset along the long axis.
    chamber_scale = float(np.cbrt(spec.chamber_fraction / 2.0))
    cham_semi = tuple(s * chamber_scale for s in heart_semi)
    offset = heart_semi[0] * 0.40
    for sign in (-1.0, +1.0):
        c = (centre[0] + sign * offset, centre[1], centre[2])
        cav = _superellipsoid(grid, c, cham_semi, 2.0) & heart
        labels[cav] = CHAMBER

    return Phantom(
        labels=labels, spacing_m=spacing_mm * 1e-3, spec=spec, pillow_exponent=p
    )


def rat_heart_cryobag(spacing_mm: float = 0.4) -> Phantom:
    return build_phantom(RAT_HEART_CRYOBAG, spacing_mm=spacing_mm)

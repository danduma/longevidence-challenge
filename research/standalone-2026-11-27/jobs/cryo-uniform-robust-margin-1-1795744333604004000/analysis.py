"""Independent NumPy/SciPy forward audit of Marc Davici's Apache-2.0 Vitrify.

Geometry and released comparator inputs are from the attributed source snapshot.
Thermal equations and weak form are reimplemented from the pinned source.
"""
import json
import math
from itertools import product
from pathlib import Path
import sys

import numpy as np
from scipy.special import logsumexp
from scipy.linalg import cholesky_banded, cho_solve_banded

ROOT = Path(__file__).resolve().parent
UPSTREAM = ROOT / "sources/upstream"
sys.path.insert(0, str(UPSTREAM / "src"))
from vitrify.geometry import rat_heart_cryobag


def read(name):
    return json.loads((UPSTREAM / "results" / name).read_text())


def material(t, wall):
    cp = np.select([t < -120, t <= -119],
                   [2402.3 + 7.6*t, 66747 + 545.24*t],
                   2678.8 + 11.73*t + 0.039*t*t)
    k = np.select([t <= -116, t < -96, t <= -77, t < -34],
                  [0.440 + 0.000816*t, 0.234 - 0.001*t,
                   0.374 + 0.000459*t, 0.374 + 0.0005*t],
                  0.374 + 0.000459*t)
    wall_cp = 1031 + 3.57*t - 0.014*t*t - 6.171e-5*t**3
    wall_rho = 931 - 0.49*t - 0.00091*t*t
    return np.where(wall, wall_cp*wall_rho, cp*(1078-0.43*t)), np.where(wall, 0.39+0.00052*t, k)


def thermal(phantom, loading, duty):
    labels = phantom.labels
    active = labels >= 0
    wall = labels[active] == 0
    heart = labels[active] >= 2
    specimen = ~wall
    index = np.full(labels.shape, -1)
    index[active] = np.arange(active.sum())
    face_i, face_j = [], []
    for axis in range(3):
        left = np.take(index, range(labels.shape[axis]-1), axis)
        right = np.take(index, range(1, labels.shape[axis]), axis)
        pair = (left >= 0) & (right >= 0)
        face_i.extend(left[pair])
        face_j.extend(right[pair])
    i, j = np.array(face_i), np.array(face_j)
    sweep = np.linspace(-180, 20, 400)
    cap, k = material(np.tile(sweep, 2), np.repeat([False, True], 400))
    dx = phantom.spacing_m
    dt = 0.4*cap.min()*dx**2/(6*k.max())
    steps70 = math.ceil(70/dt)
    schedule = np.interp(np.linspace(0, 1, steps70), [0, .5, 1], duty)
    temperatures = [np.full(active.sum(), -150.)]
    rates = []
    window_rates = []
    glass_rates = []
    endpoints = {}
    for step in range(math.ceil(100/dt)):
        t = temperatures[-1]
        cap, k = material(t, wall)
        flux = (2*k[i]*k[j]/(k[i]+k[j]))*(t[j]-t[i])/dx**2
        q = 1000*loading[active]*(691-4.6*(np.minimum(t, 0)+80))
        dtdt = (q*schedule[min(step, steps70-1)] +
                np.bincount(i, flux, minlength=len(t)) -
                np.bincount(j, flux, minlength=len(t))) / cap
        rates.append(dtdt[specimen].min()*60)
        above_glass = heart & (t >= -123)
        glass_rates.append(dtdt[above_glass].min()*60 if above_glass.any() else np.inf)
        temperature_window = above_glass & (t <= -38)
        window_rates.append(dtdt[temperature_window].min()*60 if temperature_window.any() else np.inf)
        next_t = t + dt*dtdt
        temperatures.append(next_t)
        for name, mask in (("heart", heart), ("specimen", specimen)):
            if name not in endpoints and next_t[mask].min() >= -35.:
                # Linear interpolation of each cell, then take the last crossing.
                cold = mask & (t < -35.)
                fraction = np.max((-35-t[cold])/(next_t[cold]-t[cold]))
                event_t = t + fraction*(next_t-t)
                endpoints[name] = ((step+fraction)*dt, event_t)
        if len(endpoints) == 2 and step+1 >= steps70:
            break
    assert len(endpoints) == 2, "Endpoint missed within thermal integration bound"
    return dt, steps70, np.array(temperatures), np.array(rates), np.array(window_rates), np.array(glass_rates), endpoints


class Mechanics:
    """Q1 hex equilibrium, cell-centre stress and DG0 Maxwell history."""

    def __init__(self, phantom):
        self.specimen = phantom.specimen_mask.ravel()
        self.active = phantom.active_mask.ravel()
        shape = np.array(phantom.shape)
        cells = np.indices(tuple(shape)).reshape(3, -1).T
        corners = np.array(list(product((0, 1), repeat=3)))
        node_shape = tuple(shape + 1)
        nodes = np.ravel_multi_index((cells[:, None]+corners).transpose(2, 0, 1), node_shape)
        self.dofs = (nodes[:, :, None]*3 + np.arange(3)).reshape(-1, 24)
        self.ndof = 3*math.prod(node_shape)
        anchor = np.round(shape/2).astype(int)
        pins = [(anchor, c) for c in range(3)]
        pins += [(anchor+[2, 0, 0], c) for c in (1, 2)]
        pins += [(anchor+[0, 2, 0], 2)]
        self.free = np.ones(self.ndof, dtype=bool)
        for point, component in pins:
            self.free[3*np.ravel_multi_index(point, node_shape)+component] = False
        rows = np.repeat(self.dofs, 24, axis=1).ravel()
        cols = np.tile(self.dofs, (1, 24)).ravel()
        self.lower = self.free[rows] & self.free[cols] & (rows >= cols)
        renumber = np.cumsum(self.free)-1
        row = renumber[rows[self.lower]]
        col = renumber[cols[self.lower]]
        self.free_count = int(self.free.sum())
        self.band_rows = int((row-col).max()+1)
        self.band_indices = (row-col)*self.free_count+col

        def strain(point):
            # Affine length scaling cancels if displacement is in cell lengths.
            signs = 2*corners-1
            factors = np.where(corners, point, 1-point)
            grad = np.column_stack([
                signs[:, a]*np.prod(np.delete(factors, a, axis=1), axis=1)
                for a in range(3)
            ])
            gradient = np.einsum("na,bc->banc", grad, np.eye(3)).reshape(3, 3, 24)
            sym = (gradient + gradient.swapaxes(0, 1))/2
            trace = np.trace(sym, axis1=0, axis2=1)
            dev = sym - np.eye(3)[:, :, None]*trace/3
            return dev.reshape(9, 24), trace

        self.dev, self.trace = strain(np.full(3, .5))
        self.kdev = np.zeros((24, 24))
        self.kvol = np.zeros((24, 24))
        for p in product((.5-.5/math.sqrt(3), .5+.5/math.sqrt(3)), repeat=3):
            dev, trace = strain(np.array(p))
            self.kdev += 2*dev.T@dev/8
            self.kvol += np.outer(trace, trace)/8
        self.bulk = np.where(self.specimen, 800e6/1.5, 800e6/1.5*1e-7)
        self.kbulk = self.bulk[:, None, None]*self.kvol

    def solve(self, times, field):
        previous_t = np.full(len(self.specimen), -150.)
        previous_e = np.zeros((len(self.specimen), 9))
        previous_s = np.zeros_like(previous_e)
        previous_time = 0.
        principal = []
        last_shear = None
        for time, values in zip(times, field, strict=True):
            t = previous_t.copy()
            t[self.active] = values
            midpoint = (t+previous_t)/2
            eta = np.select([midpoint >= -100, midpoint >= -140],
                            [1.21e4, 4.2783e-23*np.exp(-.6091*midpoint)], 4.63e14)
            ratio = (time-previous_time)*320e6/eta
            decay = np.exp(-ratio)
            shear = 320e6*(-np.expm1(-ratio))/ratio
            shear[~self.specimen] = 320e6*1e-7
            history = decay[:, None]*previous_s - 2*shear[:, None]*previous_e
            eigenstrain_trace = np.where(self.specimen, 3*1.1e-4*(t+150), 0)
            load = (self.bulk*eigenstrain_trace)[:, None]*self.trace - history@self.dev
            rhs = np.bincount(self.dofs.ravel(), load.ravel(), minlength=self.ndof)
            if last_shear is None or not np.allclose(shear, last_shear, rtol=1e-13, atol=0):
                local = self.kbulk + shear[:, None, None]*self.kdev
                band = np.bincount(self.band_indices, local.ravel()[self.lower],
                                   minlength=self.band_rows*self.free_count)
                factor = cholesky_banded(band.reshape(self.band_rows, self.free_count),
                                         lower=True, check_finite=False)
                last_shear = shear
            displacement = np.zeros(self.ndof)
            displacement[self.free] = cho_solve_banded((factor, True), rhs[self.free], check_finite=False)
            cell_disp = displacement[self.dofs]
            strain = cell_disp@self.dev.T
            previous_s = history + 2*shear[:, None]*strain
            previous_s[~self.specimen] = 0.
            pressure = self.bulk*(cell_disp@self.trace - eigenstrain_trace)
            stress = previous_s.reshape(-1, 3, 3) + pressure[:, None, None]*np.eye(3)
            principal.append(np.linalg.eigvalsh(stress[self.specimen]))
            previous_e = strain
            previous_e[~self.specimen] = 0.
            previous_t, previous_time = t, time
        return np.array(principal)


def summarize(times, principal):
    peak_by_time = principal[:, :, -1].max(axis=1)
    return dict(peak_mpa=float(peak_by_time.max()/1e6),
                peak_time_s=float(times[peak_by_time.argmax()]),
                objective_mpa=float(.02*logsumexp(principal/.02e6)),
                steps=len(times))


def design_fields(phantom, parameters):
    shape = np.asarray(phantom.shape)
    coords = np.indices(phantom.shape) + .5
    reshape = (3, 1, 1, 1)
    axes = np.asarray(phantom.spec.heart_dims_mm)/(2000*phantom.spacing_m)
    radius = np.sqrt(np.sum(((coords-shape.reshape(reshape)/2)/axes.reshape(reshape))**2, axis=0))
    centres, widths = [[1.05, 1.45, 2.0], [.32, .42, .65]]
    bands = np.exp(-.5*((radius-np.array(centres).reshape(reshape)) /
                         np.array(widths).reshape(reshape))**2)
    field = np.tensordot(parameters[:3], bands/bands.sum(axis=0), axes=1)
    return np.select([phantom.labels == 1, phantom.labels == 2, phantom.labels == 3],
                     [field, 1.47, 10.])


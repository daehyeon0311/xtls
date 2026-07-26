"""Broadening and file output for calculated spectra."""

from __future__ import annotations

from pathlib import Path

import numpy as np


GAUSSIAN_HWHM_TO_SIGMA = 1.0 / np.sqrt(2.0 * np.log(2.0))


def gaussian_broaden(y: np.ndarray, energy_step: float, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return y.copy()
    half_width = max(1, int(np.ceil(5 * sigma / abs(energy_step))))
    offsets = np.arange(-half_width, half_width + 1)
    kernel = np.exp(-0.5 * (offsets * energy_step / sigma) ** 2)
    kernel /= np.sum(kernel)
    return np.convolve(y, kernel, mode="same")


def variable_gaussian_broaden(
    y: np.ndarray,
    energy: np.ndarray,
    sigma: float | np.ndarray,
) -> np.ndarray:
    """Gaussian broadening with a constant or energy-dependent sigma.

    This mirrors the X-code `spc.x`/`spclib.x` idea where widths can be given
    as energy-dependent tables. The kernel is normalized at every source point,
    so integrated intensity is stable on a uniform grid.
    """
    energy = np.asarray(energy, dtype=float)
    y = np.asarray(y, dtype=float)
    sigma_array = np.asarray(sigma, dtype=float)
    if energy.shape != y.shape:
        raise ValueError("energy and y must have the same shape")
    if len(energy) < 2:
        return y.copy()
    if sigma_array.ndim == 0:
        if float(sigma_array) <= 0.0:
            return y.copy()
        return gaussian_broaden(y, float(energy[1] - energy[0]), float(sigma_array))
    if sigma_array.shape != energy.shape:
        raise ValueError("array sigma must have the same shape as energy")
    if np.all(sigma_array <= 0.0):
        return y.copy()

    out = np.zeros_like(y, dtype=float)
    step = abs(float(energy[1] - energy[0]))
    for index, width in enumerate(sigma_array):
        if width <= 0.0 or y[index] == 0.0:
            out[index] += y[index]
            continue
        half_width = max(1, int(np.ceil(5.0 * width / step)))
        lo = max(0, index - half_width)
        hi = min(len(y), index + half_width + 1)
        offsets = energy[lo:hi] - energy[index]
        kernel = np.exp(-0.5 * (offsets / width) ** 2)
        kernel /= np.sum(kernel)
        out[lo:hi] += y[index] * kernel
    return out


def interpolate_width(
    energy: np.ndarray,
    base_width: float,
    points: list[tuple[float, float]] | tuple[tuple[float, float], ...] | np.ndarray,
) -> float | np.ndarray:
    """Return constant or linearly interpolated width values on `energy`.

    `points` should contain `(energy, width)` pairs in eV. Outside the table,
    the nearest endpoint width is used, matching the practical behavior of
    XTLS/X-code broadening tables.
    """
    if len(points) == 0:
        return float(base_width)
    table = np.asarray(points, dtype=float)
    if table.ndim != 2 or table.shape[1] != 2:
        raise ValueError("width points must be (energy, width) pairs")
    order = np.argsort(table[:, 0])
    table = table[order]
    return np.interp(np.asarray(energy, dtype=float), table[:, 0], table[:, 1])


def save_xy(path: str | Path, energy: np.ndarray, columns: dict[str, np.ndarray]) -> None:
    names = list(columns)
    data = np.column_stack([energy, *[columns[name] for name in names]])
    header = "energy " + " ".join(names)
    np.savetxt(path, data, header=header, comments="")

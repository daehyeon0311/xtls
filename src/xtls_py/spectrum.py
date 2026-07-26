from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


GAUSSIAN_HWHM_TO_SIGMA = 1.0 / np.sqrt(2.0 * np.log(2.0))


@dataclass(frozen=True)
class SpectrumBlock:
    polarization: str
    initial_index: int
    energy: np.ndarray
    intensity: np.ndarray
    gamma: float


@dataclass(frozen=True)
class XtlsSpectrum:
    mode: str
    polarizations: list[str]
    initial_energies: np.ndarray
    blocks: list[SpectrumBlock]

    def by_polarization(self) -> dict[str, np.ndarray]:
        grouped: dict[str, list[np.ndarray]] = {}
        for block in self.blocks:
            grouped.setdefault(block.polarization, []).append(block.intensity)
        return {key: np.sum(values, axis=0) for key, values in grouped.items()}

    @property
    def energy(self) -> np.ndarray:
        if not self.blocks:
            return np.array([])
        return self.blocks[0].energy

    def isotropic(self) -> np.ndarray:
        grouped = self.by_polarization()
        if all(key in grouped for key in ("x", "y", "z")):
            return (grouped["x"] + grouped["y"] + grouped["z"]) / 3.0
        if grouped:
            return np.sum(list(grouped.values()), axis=0) / len(grouped)
        return np.array([])


def read_xtls_obj(path: str | Path) -> XtlsSpectrum:
    lines = Path(path).read_text(errors="replace").splitlines()
    mode = _read_scalar_after(lines, "MODE", default="").strip().lower()
    polarizations = _read_list_after(lines, "DICRO1") or ["total"]
    initial_energies = np.array(_read_numeric_list_after(lines, "ENERGY"), dtype=float)
    blocks = _read_spcd_blocks(lines, polarizations, max(1, len(initial_energies)))
    return XtlsSpectrum(mode, polarizations, initial_energies, blocks)


def gaussian_broaden(y: np.ndarray, energy_step: float, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return y.copy()
    half_width = max(1, int(np.ceil(5 * sigma / abs(energy_step))))
    offsets = np.arange(-half_width, half_width + 1)
    kernel = np.exp(-0.5 * (offsets * energy_step / sigma) ** 2)
    kernel /= np.sum(kernel)
    return np.convolve(y, kernel, mode="same")


def gaussian_hwhm_broaden(y: np.ndarray, energy_step: float, hwhm: float) -> np.ndarray:
    """Gaussian broadening with XTLS-style HWHM input."""
    return gaussian_broaden(y, energy_step, hwhm * GAUSSIAN_HWHM_TO_SIGMA)


def lorentzian_broaden(y: np.ndarray, energy_step: float, hwhm: float) -> np.ndarray:
    if hwhm <= 0:
        return y.copy()
    half_width = max(1, int(np.ceil(50 * hwhm / abs(energy_step))))
    offsets = np.arange(-half_width, half_width + 1)
    kernel = 1.0 / (1.0 + (offsets * energy_step / hwhm) ** 2)
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


def broaden(
    y: np.ndarray,
    energy: np.ndarray,
    gaussian_sigma: float = 0.0,
    gaussian_hwhm: float = 0.0,
    lorentzian_hwhm: float = 0.0,
) -> np.ndarray:
    if len(energy) < 2:
        return y.copy()
    step = float(energy[1] - energy[0])
    sigma = gaussian_sigma
    if sigma <= 0.0 and gaussian_hwhm > 0.0:
        sigma = gaussian_hwhm * GAUSSIAN_HWHM_TO_SIGMA
    out = gaussian_broaden(y, step, sigma)
    out = lorentzian_broaden(out, step, lorentzian_hwhm)
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


def _read_scalar_after(lines: list[str], marker: str, default: str = "") -> str:
    for index, line in enumerate(lines):
        if line.strip().upper() == marker:
            for value in lines[index + 1 :]:
                stripped = value.strip()
                if stripped:
                    return stripped
    return default


def _read_list_after(lines: list[str], marker: str) -> list[str]:
    for index, line in enumerate(lines):
        if line.strip().upper() != marker:
            continue
        count = int(lines[index + 1].strip())
        labels: list[str] = []
        cursor = index + 2
        while cursor < len(lines) and len(labels) < count:
            stripped = lines[cursor].strip()
            if stripped and stripped not in ("{", "}"):
                labels.append(stripped.lower())
            cursor += 1
        return labels
    return []


def _read_numeric_list_after(lines: list[str], marker: str) -> list[float]:
    for index, line in enumerate(lines):
        if line.strip().upper() != marker:
            continue
        count = int(lines[index + 1].strip())
        values: list[float] = []
        cursor = index + 2
        while cursor < len(lines) and len(values) < count:
            values.extend(_numbers(lines[cursor]))
            cursor += 1
        return values[:count]
    return []


def _read_spcd_blocks(
    lines: list[str],
    polarizations: list[str],
    initial_count: int,
) -> list[SpectrumBlock]:
    blocks: list[SpectrumBlock] = []
    spcd_count = 0
    for index, line in enumerate(lines):
        if line.strip().upper() != "SPCD":
            continue
        header = _numbers(lines[index + 1])
        if len(header) < 4:
            continue
        n_points = int(header[0]) + 1
        energy_min, energy_max, gamma = header[1], header[2], header[3]
        values: list[float] = []
        cursor = index + 2
        while cursor < len(lines) and len(values) < n_points:
            values.extend(_numbers(lines[cursor]))
            cursor += 1
        pol = polarizations[(spcd_count // initial_count) % max(1, len(polarizations))]
        initial_index = spcd_count % initial_count
        blocks.append(
            SpectrumBlock(
                polarization=pol,
                initial_index=initial_index,
                energy=np.linspace(energy_min, energy_max, n_points),
                intensity=np.array(values[:n_points], dtype=float),
                gamma=gamma,
            )
        )
        spcd_count += 1
    return blocks


def _numbers(line: str) -> list[float]:
    out: list[float] = []
    for part in line.replace("D", "E").replace("d", "e").split():
        try:
            out.append(float(part))
        except ValueError:
            pass
    return out

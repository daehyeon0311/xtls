from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class StickSpectrum:
    """Discrete transition spectrum."""

    energy: np.ndarray
    intensity: np.ndarray

    def shifted(self, offset: float) -> "StickSpectrum":
        return StickSpectrum(self.energy + offset, self.intensity.copy())

    def normalized(self, mode: str = "max") -> "StickSpectrum":
        if mode == "none":
            return StickSpectrum(self.energy.copy(), self.intensity.copy())
        if mode == "max":
            scale = float(np.max(self.intensity)) if self.intensity.size else 0.0
        elif mode == "area":
            scale = float(np.sum(self.intensity))
        else:
            raise ValueError("mode must be 'max', 'area', or 'none'")
        if scale == 0.0:
            return StickSpectrum(self.energy.copy(), self.intensity.copy())
        return StickSpectrum(self.energy.copy(), self.intensity / scale)


def exact_transitions(
    h_initial: np.ndarray,
    h_final: np.ndarray,
    transition: np.ndarray,
    n_initial_states: int = 1,
    initial_weights: np.ndarray | None = None,
) -> StickSpectrum:
    """Compute stick transitions between exact eigenstates.

    `transition` must have shape `(dim_final, dim_initial)` and represents the
    operator connecting the initial Hilbert space to the final Hilbert space.
    By default only the initial ground state contributes.
    """
    h_initial = np.asarray(h_initial)
    h_final = np.asarray(h_final)
    transition = np.asarray(transition)
    if h_initial.shape[0] != h_initial.shape[1]:
        raise ValueError("h_initial must be square")
    if h_final.shape[0] != h_final.shape[1]:
        raise ValueError("h_final must be square")
    if transition.shape != (h_final.shape[0], h_initial.shape[0]):
        raise ValueError("transition must have shape (dim_final, dim_initial)")
    if not 1 <= n_initial_states <= h_initial.shape[0]:
        raise ValueError("n_initial_states must fit the initial Hilbert space")

    e_initial, v_initial = np.linalg.eigh(h_initial)
    e_final, v_final = np.linalg.eigh(h_final)

    if initial_weights is None:
        weights = np.zeros(n_initial_states, dtype=float)
        weights[0] = 1.0
    else:
        weights = np.asarray(initial_weights, dtype=float)
        if weights.shape != (n_initial_states,):
            raise ValueError("initial_weights must match n_initial_states")
        total = float(np.sum(weights))
        if total <= 0.0:
            raise ValueError("initial_weights must have positive total weight")
        weights = weights / total

    energies: list[float] = []
    intensities: list[float] = []
    for initial_index in range(n_initial_states):
        initial_vec = v_initial[:, initial_index]
        amplitudes = v_final.conjugate().T @ transition @ initial_vec
        strengths = weights[initial_index] * np.abs(amplitudes) ** 2
        for final_index, strength in enumerate(strengths):
            if strength <= 0.0:
                continue
            energies.append(float(e_final[final_index] - e_initial[initial_index]))
            intensities.append(float(strength))

    order = np.argsort(energies)
    return StickSpectrum(np.array(energies)[order], np.array(intensities)[order])


def sticks_to_curve(
    sticks: StickSpectrum,
    energy_grid: np.ndarray,
    gaussian_sigma: float = 0.0,
    lorentzian_hwhm: float = 0.0,
    normalize: bool = True,
) -> np.ndarray:
    """Broaden a stick spectrum onto a fixed energy grid."""
    energy_grid = np.asarray(energy_grid, dtype=float)
    curve = np.zeros_like(energy_grid)
    for energy, intensity in zip(sticks.energy, sticks.intensity):
        if gaussian_sigma > 0.0:
            curve += intensity * np.exp(-0.5 * ((energy_grid - energy) / gaussian_sigma) ** 2)
        elif lorentzian_hwhm <= 0.0:
            nearest = int(np.argmin(np.abs(energy_grid - energy)))
            curve[nearest] += intensity
        if lorentzian_hwhm > 0.0:
            curve += intensity / (1.0 + ((energy_grid - energy) / lorentzian_hwhm) ** 2)
    max_value = float(np.max(curve)) if curve.size else 0.0
    if normalize and max_value > 0.0:
        curve = curve / max_value
    return curve

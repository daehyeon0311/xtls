from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class LanczosTridiagonal:
    """Tridiagonal Lanczos representation of a Hermitian operator."""

    alpha: np.ndarray
    beta: np.ndarray
    start_norm: float


def lanczos_tridiagonal(
    operator,
    start_vector: np.ndarray,
    max_iter: int = 200,
    tol: float = 1e-12,
    reorthogonalize: bool = False,
) -> LanczosTridiagonal:
    """Build a Lanczos tridiagonal from `start_vector`.

    `operator` may be a dense matrix, sparse matrix, or callable matvec.
    """
    if max_iter <= 0:
        raise ValueError("max_iter must be positive")
    matvec = _as_matvec(operator)
    start = np.asarray(start_vector, dtype=complex)
    start_norm = float(np.linalg.norm(start))
    if start_norm <= tol:
        return LanczosTridiagonal(np.array([]), np.array([]), 0.0)

    q = start / start_norm
    q_prev = np.zeros_like(q)
    beta_prev = 0.0
    basis_vectors: list[np.ndarray] = []
    alpha: list[float] = []
    beta: list[float] = []

    for iteration in range(max_iter):
        z = np.asarray(matvec(q), dtype=complex)
        a = np.vdot(q, z)
        z = z - a * q - beta_prev * q_prev
        if reorthogonalize:
            for old_q in basis_vectors:
                z = z - np.vdot(old_q, z) * old_q
            basis_vectors.append(q.copy())
        b = float(np.linalg.norm(z))
        alpha.append(float(np.real_if_close(a)))
        if b <= tol or iteration == max_iter - 1:
            break
        beta.append(b)
        q_prev = q
        q = z / b
        beta_prev = b

    return LanczosTridiagonal(
        alpha=np.array(alpha, dtype=float),
        beta=np.array(beta, dtype=float),
        start_norm=start_norm,
    )


def continued_fraction_green(
    tridiagonal: LanczosTridiagonal,
    energy: np.ndarray,
    reference_energy: float = 0.0,
    broadening: float | np.ndarray = 0.2,
) -> np.ndarray:
    """Evaluate <b|(E_ref + omega + i eta - H)^-1|b>."""
    broadening_array = np.asarray(broadening, dtype=float)
    if np.any(broadening_array <= 0.0):
        raise ValueError("broadening must be positive")
    energy = np.asarray(energy, dtype=float)
    if broadening_array.ndim > 0 and broadening_array.shape != energy.shape:
        raise ValueError("array broadening must have the same shape as energy")
    if tridiagonal.start_norm == 0.0 or tridiagonal.alpha.size == 0:
        return np.zeros_like(energy, dtype=complex)

    z = reference_energy + energy + 1j * broadening_array
    value = 1.0 / (z - tridiagonal.alpha[-1])
    for index in range(tridiagonal.alpha.size - 2, -1, -1):
        value = 1.0 / (z - tridiagonal.alpha[index] - tridiagonal.beta[index] ** 2 * value)
    return tridiagonal.start_norm**2 * value


def continued_fraction_spectrum(
    operator,
    start_vector: np.ndarray,
    energy: np.ndarray,
    reference_energy: float = 0.0,
    broadening: float | np.ndarray = 0.2,
    max_iter: int = 200,
    tol: float = 1e-12,
    reorthogonalize: bool = False,
    normalize: bool = False,
) -> np.ndarray:
    """Return the spectral function from a Lanczos continued fraction.

    The result is `-Im G(omega) / pi`, where
    `G = <b|(E_ref + omega + i eta - H)^-1|b>`.
    """
    tri = lanczos_tridiagonal(
        operator,
        start_vector,
        max_iter=max_iter,
        tol=tol,
        reorthogonalize=reorthogonalize,
    )
    green = continued_fraction_green(
        tri,
        energy,
        reference_energy=reference_energy,
        broadening=broadening,
    )
    spectrum = -np.imag(green) / np.pi
    if normalize:
        scale = float(np.max(np.abs(spectrum))) if spectrum.size else 0.0
        if scale > 0.0:
            spectrum = spectrum / scale
    return spectrum


def _as_matvec(operator) -> Callable[[np.ndarray], np.ndarray]:
    if callable(operator):
        return operator
    return operator.dot if hasattr(operator, "dot") else lambda vector: operator @ vector

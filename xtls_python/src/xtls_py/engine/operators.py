from __future__ import annotations

import numpy as np

from .basis import FockBasis


def create(state: int, orbital: int) -> tuple[int, int] | None:
    """Apply creation operator c^dagger_orbital to `state`.

    Returns `(sign, new_state)` or `None` when the orbital is already occupied.
    """
    if (state >> orbital) & 1:
        return None
    sign = _fermion_sign(state, orbital)
    return sign, state | (1 << orbital)


def annihilate(state: int, orbital: int) -> tuple[int, int] | None:
    """Apply c_orbital to `state`.

    Returns `(sign, new_state)` or `None` when the orbital is empty.
    """
    if not ((state >> orbital) & 1):
        return None
    sign = _fermion_sign(state, orbital)
    return sign, state & ~(1 << orbital)


def number_matrix(basis: FockBasis, orbital: int) -> np.ndarray:
    """Matrix representation of n_orbital."""
    matrix = np.zeros((len(basis), len(basis)), dtype=float)
    for col, state in enumerate(basis.states):
        matrix[col, col] = (state >> orbital) & 1
    return matrix


def density_density_matrix(basis: FockBasis, interaction: np.ndarray) -> np.ndarray:
    """Build matrix for sum_{i<j} interaction[i,j] n_i n_j.

    This is a first two-body building block. It is not yet the full Slater
    multiplet Coulomb tensor.
    """
    interaction = np.asarray(interaction)
    if interaction.shape != (basis.n_orbitals, basis.n_orbitals):
        raise ValueError("interaction must have shape (n_orbitals, n_orbitals)")
    matrix = np.zeros((len(basis), len(basis)), dtype=float)
    for col, state in enumerate(basis.states):
        energy = 0.0
        for i in range(basis.n_orbitals):
            if not ((state >> i) & 1):
                continue
            for j in range(i + 1, basis.n_orbitals):
                if (state >> j) & 1:
                    energy += float(interaction[i, j])
        matrix[col, col] = energy
    return matrix


def two_body_matrix(basis: FockBasis, tensor: np.ndarray, prefactor: float = 0.5) -> np.ndarray:
    """Build matrix for prefactor * sum_ijkl V[i,j,k,l] c_i^dag c_j^dag c_l c_k.

    The tensor convention follows the common Coulomb form
    <ij|V|kl> c_i^dag c_j^dag c_l c_k. The default prefactor avoids double
    counting when the supplied tensor includes both exchange-related terms.
    """
    tensor = np.asarray(tensor)
    shape = (basis.n_orbitals, basis.n_orbitals, basis.n_orbitals, basis.n_orbitals)
    if tensor.shape != shape:
        raise ValueError("tensor must have shape (n_orbitals, n_orbitals, n_orbitals, n_orbitals)")

    matrix = np.zeros((len(basis), len(basis)), dtype=np.result_type(tensor, complex))
    nonzero = np.argwhere(tensor != 0)
    for col, state in enumerate(basis.states):
        for i, j, k, l in nonzero:
            value = tensor[i, j, k, l]
            first = annihilate(state, k)
            if first is None:
                continue
            sign_k, state_k = first
            second = annihilate(state_k, l)
            if second is None:
                continue
            sign_l, state_kl = second
            third = create(state_kl, j)
            if third is None:
                continue
            sign_j, state_jkl = third
            fourth = create(state_jkl, i)
            if fourth is None:
                continue
            sign_i, new_state = fourth
            row = basis.index.get(new_state)
            if row is not None:
                matrix[row, col] += prefactor * value * sign_k * sign_l * sign_j * sign_i
    return _real_if_close(matrix)


def one_body_matrix(basis: FockBasis, h: np.ndarray) -> np.ndarray:
    """Build many-body matrix for sum_ij h[i,j] c_i^dagger c_j."""
    h = np.asarray(h)
    if h.shape != (basis.n_orbitals, basis.n_orbitals):
        raise ValueError("h must have shape (n_orbitals, n_orbitals)")
    matrix = np.zeros((len(basis), len(basis)), dtype=np.result_type(h, complex))
    for col, state in enumerate(basis.states):
        for j in range(basis.n_orbitals):
            ann = annihilate(state, j)
            if ann is None:
                continue
            sign_ann, intermediate = ann
            for i in range(basis.n_orbitals):
                if h[i, j] == 0:
                    continue
                cre = create(intermediate, i)
                if cre is None:
                    continue
                sign_cre, new_state = cre
                row = basis.index.get(new_state)
                if row is not None:
                    matrix[row, col] += h[i, j] * sign_ann * sign_cre
    return _real_if_close(matrix)


def one_body_between_bases(
    final_basis: FockBasis,
    initial_basis: FockBasis,
    h: np.ndarray,
) -> np.ndarray:
    """Build matrix for sum_ij h[i,j] c_i^dagger c_j between two bases."""
    h = np.asarray(h)
    if final_basis.n_orbitals != initial_basis.n_orbitals:
        raise ValueError("bases must use the same orbital space")
    if final_basis.n_electrons != initial_basis.n_electrons:
        raise ValueError("one-body operators conserve total electron count")
    if h.shape != (initial_basis.n_orbitals, initial_basis.n_orbitals):
        raise ValueError("h must have shape (n_orbitals, n_orbitals)")

    matrix = np.zeros((len(final_basis), len(initial_basis)), dtype=np.result_type(h, complex))
    for col, state in enumerate(initial_basis.states):
        for j in range(initial_basis.n_orbitals):
            ann = annihilate(state, j)
            if ann is None:
                continue
            sign_ann, intermediate = ann
            for i in range(initial_basis.n_orbitals):
                if h[i, j] == 0:
                    continue
                cre = create(intermediate, i)
                if cre is None:
                    continue
                sign_cre, new_state = cre
                row = final_basis.index.get(new_state)
                if row is not None:
                    matrix[row, col] += h[i, j] * sign_ann * sign_cre
    return _real_if_close(matrix)


def _fermion_sign(state: int, orbital: int) -> int:
    lower_mask = (1 << orbital) - 1
    lower_count = (state & lower_mask).bit_count()
    return -1 if lower_count % 2 else 1


def _real_if_close(matrix: np.ndarray) -> np.ndarray:
    if np.allclose(matrix.imag, 0):
        return matrix.real
    return matrix

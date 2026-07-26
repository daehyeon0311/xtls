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


def two_body_matrix(basis: FockBasis, tensor: np.ndarray, prefactor: float = 0.5) -> np.ndarray:
    """Build matrix for prefactor * sum_ijkl V[i,j,k,l] c_i^dag c_j^dag c_l c_k.

    The tensor convention follows the common Coulomb form
    <ij|V|kl> c_i^dag c_j^dag c_l c_k. The default prefactor avoids double
    counting when the supplied tensor includes both exchange-related terms.

    This is the dense counterpart of `two_body_sparse`; the two share one
    kernel, so they agree by construction.
    """
    from .sparse import two_body_sparse

    return _dense(two_body_sparse(basis, tensor, prefactor=prefactor))


def one_body_matrix(basis: FockBasis, h: np.ndarray) -> np.ndarray:
    """Build many-body matrix for sum_ij h[i,j] c_i^dagger c_j."""
    from .sparse import one_body_sparse

    return _dense(one_body_sparse(basis, h))


def one_body_between_bases(
    final_basis: FockBasis,
    initial_basis: FockBasis,
    h: np.ndarray,
) -> np.ndarray:
    """Build matrix for sum_ij h[i,j] c_i^dagger c_j between two bases."""
    from .sparse import one_body_between_bases_sparse

    return _dense(one_body_between_bases_sparse(final_basis, initial_basis, h))


def _dense(matrix) -> np.ndarray:
    return _real_if_close(np.asarray(matrix.todense()))


def _fermion_sign(state: int, orbital: int) -> int:
    lower_mask = (1 << orbital) - 1
    lower_count = (state & lower_mask).bit_count()
    return -1 if lower_count % 2 else 1


def _real_if_close(matrix: np.ndarray) -> np.ndarray:
    if np.allclose(matrix.imag, 0):
        return matrix.real
    return matrix

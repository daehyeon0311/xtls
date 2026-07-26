from __future__ import annotations

import numpy as np

from .basis import FockBasis
from .operators import annihilate, create


def one_body_sparse(basis: FockBasis, h):
    """Sparse matrix for sum_ij h[i,j] c_i^dagger c_j."""
    sparse = _require_scipy_sparse()
    h = np.asarray(h)
    if h.shape != (basis.n_orbitals, basis.n_orbitals):
        raise ValueError("h must have shape (n_orbitals, n_orbitals)")
    rows: list[int] = []
    cols: list[int] = []
    data: list[complex] = []
    nonzero = np.argwhere(h != 0)
    for col, state in enumerate(basis.states):
        for i, j in nonzero:
            ann = annihilate(state, j)
            if ann is None:
                continue
            sign_ann, intermediate = ann
            cre = create(intermediate, i)
            if cre is None:
                continue
            sign_cre, new_state = cre
            row = basis.index.get(new_state)
            if row is not None:
                rows.append(row)
                cols.append(col)
                data.append(h[i, j] * sign_ann * sign_cre)
    return sparse.coo_matrix((data, (rows, cols)), shape=(len(basis), len(basis))).tocsr()


def one_body_between_bases_sparse(final_basis: FockBasis, initial_basis: FockBasis, h):
    """Sparse one-body matrix connecting two fixed sectors."""
    sparse = _require_scipy_sparse()
    h = np.asarray(h)
    if final_basis.n_orbitals != initial_basis.n_orbitals:
        raise ValueError("bases must use the same orbital space")
    if final_basis.n_electrons != initial_basis.n_electrons:
        raise ValueError("one-body operators conserve total electron count")
    if h.shape != (initial_basis.n_orbitals, initial_basis.n_orbitals):
        raise ValueError("h must have shape (n_orbitals, n_orbitals)")
    rows: list[int] = []
    cols: list[int] = []
    data: list[complex] = []
    nonzero = np.argwhere(h != 0)
    for col, state in enumerate(initial_basis.states):
        for i, j in nonzero:
            ann = annihilate(state, j)
            if ann is None:
                continue
            sign_ann, intermediate = ann
            cre = create(intermediate, i)
            if cre is None:
                continue
            sign_cre, new_state = cre
            row = final_basis.index.get(new_state)
            if row is not None:
                rows.append(row)
                cols.append(col)
                data.append(h[i, j] * sign_ann * sign_cre)
    return sparse.coo_matrix((data, (rows, cols)), shape=(len(final_basis), len(initial_basis))).tocsr()


def two_body_sparse(basis: FockBasis, tensor, prefactor: float = 0.5, threshold: float = 0.0):
    """Sparse matrix for prefactor * sum_ijkl V[i,j,k,l] c_i^dag c_j^dag c_l c_k."""
    sparse = _require_scipy_sparse()
    tensor = np.asarray(tensor)
    shape = (basis.n_orbitals, basis.n_orbitals, basis.n_orbitals, basis.n_orbitals)
    if tensor.shape != shape:
        raise ValueError("tensor must have shape (n_orbitals, n_orbitals, n_orbitals, n_orbitals)")
    if threshold > 0.0:
        nonzero = np.argwhere(np.abs(tensor) > threshold)
    else:
        nonzero = np.argwhere(tensor != 0)
    rows: list[int] = []
    cols: list[int] = []
    data: list[complex] = []
    for col, state in enumerate(basis.states):
        for i, j, k, l in nonzero:
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
                rows.append(row)
                cols.append(col)
                data.append(prefactor * tensor[i, j, k, l] * sign_k * sign_l * sign_j * sign_i)
    return sparse.coo_matrix((data, (rows, cols)), shape=(len(basis), len(basis))).tocsr()


def diagonal_sparse(values):
    sparse = _require_scipy_sparse()
    return sparse.diags(np.asarray(values), format="csr")


def lowest_eigenpairs(matrix, k: int = 6):
    """Return the lowest eigenpairs for dense or sparse Hermitian matrices."""
    if k <= 0:
        raise ValueError("k must be positive")
    n = matrix.shape[0]
    if k >= n:
        values, vectors = np.linalg.eigh(matrix.toarray() if hasattr(matrix, "toarray") else matrix)
        return values, vectors
    try:
        from scipy.sparse.linalg import eigsh
    except ImportError:
        values, vectors = np.linalg.eigh(matrix.toarray() if hasattr(matrix, "toarray") else matrix)
        return values[:k], vectors[:, :k]
    values, vectors = eigsh(matrix, k=k, which="SA")
    order = np.argsort(values)
    return values[order], vectors[:, order]


def _require_scipy_sparse():
    try:
        import scipy.sparse as sparse
    except ImportError as exc:
        raise ImportError("scipy is required for sparse matrix support") from exc
    return sparse

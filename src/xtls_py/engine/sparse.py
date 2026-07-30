from __future__ import annotations

import numpy as np

from .basis import FockBasis


def one_body_sparse(basis: FockBasis, h):
    """Sparse matrix for sum_ij h[i,j] c_i^dagger c_j."""
    sparse = _require_scipy_sparse()
    h = np.asarray(h)
    if h.shape != (basis.n_orbitals, basis.n_orbitals):
        raise ValueError("h must have shape (n_orbitals, n_orbitals)")
    lookup = _StateLookup(basis)
    states = lookup.states
    columns = np.arange(len(basis))

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data: list[np.ndarray] = []
    for i, j in np.argwhere(h != 0):
        moved, source, signs = _apply_chain(states, columns, (("a", j), ("c", i)))
        row, keep = lookup.index_of(moved)
        if row.size:
            rows.append(row)
            cols.append(source[keep])
            data.append(h[i, j] * signs[keep])
    return _assemble(sparse, rows, cols, data, (len(basis), len(basis)))


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
    lookup = _StateLookup(final_basis)
    states = np.asarray(initial_basis.states, dtype=np.int64)
    columns = np.arange(len(initial_basis))

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data: list[np.ndarray] = []
    for i, j in np.argwhere(h != 0):
        moved, source, signs = _apply_chain(states, columns, (("a", j), ("c", i)))
        row, keep = lookup.index_of(moved)
        if row.size:
            rows.append(row)
            cols.append(source[keep])
            data.append(h[i, j] * signs[keep])
    return _assemble(sparse, rows, cols, data, (len(final_basis), len(initial_basis)))


def annihilation_between_bases_sparse(final_basis: FockBasis, initial_basis: FockBasis, orbital: int):
    """Sparse matrix for the single annihilation operator `c_orbital`.

    This is the photoemission transition operator: it maps an `N`-electron
    basis onto an `N-1`-electron one. XTLS instead keeps the electron count
    fixed by moving the photoelectron into a continuum orbital that is left out
    of the configuration list; dropping the electron outright is equivalent and
    keeps the final basis smaller.
    """
    sparse = _require_scipy_sparse()
    if final_basis.n_orbitals != initial_basis.n_orbitals:
        raise ValueError("bases must use the same orbital space")
    if final_basis.n_electrons != initial_basis.n_electrons - 1:
        raise ValueError("annihilation must lower the electron count by exactly one")
    if not 0 <= int(orbital) < initial_basis.n_orbitals:
        raise ValueError("orbital is outside the orbital space")

    lookup = _StateLookup(final_basis)
    states = np.asarray(initial_basis.states, dtype=np.int64)
    columns = np.arange(len(initial_basis))
    moved, source, signs = _apply_chain(states, columns, (("a", int(orbital)),))
    row, keep = lookup.index_of(moved)
    shape = (len(final_basis), len(initial_basis))
    if row.size == 0:
        return sparse.csr_matrix(shape, dtype=complex)
    return _assemble(
        sparse,
        [row],
        [source[keep]],
        [signs[keep].astype(complex)],
        shape,
    )


def two_body_sparse(basis: FockBasis, tensor, prefactor: float = 0.5, threshold: float = 0.0):
    """Sparse matrix for prefactor * sum_ijkl V[i,j,k,l] c_i^dag c_j^dag c_l c_k.

    `tensor` may span fewer orbitals than the basis. Orbital indices are shared,
    so a 2p+3d tensor can be applied directly to a 2p+3d+ligand basis without
    being padded up to the full orbital count first.
    """
    sparse = _require_scipy_sparse()
    tensor = np.asarray(tensor)
    if tensor.ndim != 4 or len(set(tensor.shape)) != 1:
        raise ValueError("tensor must have shape (n, n, n, n)")
    if tensor.shape[0] > basis.n_orbitals:
        raise ValueError("tensor spans more orbitals than the basis")
    if threshold > 0.0:
        nonzero = np.argwhere(np.abs(tensor) > threshold)
    else:
        nonzero = np.argwhere(tensor != 0)

    lookup = _StateLookup(basis)
    states = lookup.states
    columns = np.arange(len(basis))

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data: list[np.ndarray] = []
    for i, j, k, l in nonzero:
        moved, source, signs = _apply_chain(
            states,
            columns,
            (("a", k), ("a", l), ("c", j), ("c", i)),
        )
        row, keep = lookup.index_of(moved)
        if row.size:
            rows.append(row)
            cols.append(source[keep])
            data.append(prefactor * tensor[i, j, k, l] * signs[keep])
    return _assemble(sparse, rows, cols, data, (len(basis), len(basis)))


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


class _StateLookup:
    """Vectorized `bit pattern -> basis row` lookup."""

    def __init__(self, basis: FockBasis) -> None:
        self.states = np.asarray(basis.states, dtype=np.int64)
        self._order = np.argsort(self.states)
        self._sorted = self.states[self._order]

    def index_of(self, states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return `(rows, keep)` for the states that exist in the basis."""
        if states.size == 0:
            return np.empty(0, dtype=np.int64), np.zeros(0, dtype=bool)
        position = np.searchsorted(self._sorted, states)
        np.clip(position, 0, self._sorted.size - 1, out=position)
        keep = self._sorted[position] == states
        return self._order[position[keep]], keep


def _apply_chain(
    states: np.ndarray,
    columns: np.ndarray,
    operations: tuple[tuple[str, int], ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply a chain of creation/annihilation operators to many states at once.

    `operations` is applied left to right, each entry being `("a", orbital)` for
    an annihilation or `("c", orbital)` for a creation. States for which an
    operator vanishes are dropped. Returns the surviving states, the column
    index they came from, and the accumulated fermion sign.
    """
    signs = np.ones(states.size, dtype=np.int8)
    for kind, orbital in operations:
        occupied = (states >> int(orbital)) & 1
        keep = occupied == (1 if kind == "a" else 0)
        if not keep.all():
            states = states[keep]
            columns = columns[keep]
            signs = signs[keep]
            if states.size == 0:
                break
        parity = _popcount(states & ((np.int64(1) << int(orbital)) - 1)) & 1
        signs = signs * np.where(parity == 1, np.int8(-1), np.int8(1))
        states = states ^ (np.int64(1) << int(orbital))
    return states, columns, signs


_POPCOUNT_BYTE = np.array([bin(value).count("1") for value in range(256)], dtype=np.int64)


def _popcount(values: np.ndarray) -> np.ndarray:
    """Count set bits. Uses numpy's builtin when available."""
    builtin = getattr(np, "bitwise_count", None)
    if builtin is not None:
        return builtin(values)
    counts = np.zeros(values.shape, dtype=np.int64)
    remaining = values
    while True:
        counts += _POPCOUNT_BYTE[remaining & 0xFF]
        remaining = remaining >> 8
        if not remaining.any():
            return counts


def _assemble(sparse, rows, cols, data, shape):
    if not rows:
        return sparse.csr_matrix(shape, dtype=complex)
    return sparse.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=shape,
    ).tocsr()


def _require_scipy_sparse():
    try:
        import scipy.sparse as sparse
    except ImportError as exc:
        raise ImportError("scipy is required for sparse matrix support") from exc
    return sparse

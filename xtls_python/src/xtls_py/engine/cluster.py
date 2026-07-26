from __future__ import annotations

import numpy as np

from .basis import FockBasis
from .configuration import (
    basis_from_configurations,
    charge_transfer_configurations,
    xas_final_configurations,
    xas_initial_configurations,
)
from .operators import one_body_matrix


def ct_spin_orbital_labels() -> tuple[str, ...]:
    """Return labels for a 3d + ligand shell model.

    The first 10 spin-orbitals are 3d-like. The next 10 are ligand-like
    symmetry-adapted orbitals with matching d labels.
    """
    d_orbitals = ("xy", "yz", "zx", "x2-y2", "3z2-r2")
    spins = ("down", "up")
    d_labels = tuple(f"d_{orb}_{spin}" for orb in d_orbitals for spin in spins)
    ligand_labels = tuple(f"L_{orb}_{spin}" for orb in d_orbitals for spin in spins)
    return d_labels + ligand_labels


def p_ct_spin_orbital_labels() -> tuple[str, ...]:
    p_orbitals = ("x", "y", "z")
    d_orbitals = ("xy", "yz", "zx", "x2-y2", "3z2-r2")
    spins = ("down", "up")
    p_labels = tuple(f"p_{orb}_{spin}" for orb in p_orbitals for spin in spins)
    d_labels = tuple(f"d_{orb}_{spin}" for orb in d_orbitals for spin in spins)
    ligand_labels = tuple(f"L_{orb}_{spin}" for orb in d_orbitals for spin in spins)
    return p_labels + d_labels + ligand_labels


def charge_transfer_basis(
    n_d_electrons: int,
    max_ligand_holes: int = 2,
    n_d_spin_orbitals: int = 10,
) -> FockBasis:
    """Build a `3d^n + 3d^(n+1)L + ...` charge-transfer basis.

    Ligand orbitals are assumed full in the ionic `3d^n` configuration. A
    ligand hole count `h` corresponds to `3d^(n+h)L^h`.
    """
    if n_d_spin_orbitals <= 0:
        raise ValueError("n_d_spin_orbitals must be positive")
    if not 0 <= n_d_electrons <= n_d_spin_orbitals:
        raise ValueError("n_d_electrons must fit the d shell")
    if max_ligand_holes < 0:
        raise ValueError("max_ligand_holes must be non-negative")

    if n_d_spin_orbitals != 10:
        raise ValueError("the general configuration engine currently uses a 10-orbital d shell")
    return basis_from_configurations(
        charge_transfer_configurations(
            n_d_electrons,
            max_ligand_holes=max_ligand_holes,
        )
    )


def p_ct_xas_initial_basis(
    n_d_electrons: int,
    max_ligand_holes: int = 1,
) -> FockBasis:
    """Initial `2p^6(3d^n + 3d^(n+1)L + ...)` basis."""
    return basis_from_configurations(
        xas_initial_configurations(
            n_d_electrons,
            max_ligand_holes=max_ligand_holes,
        )
    )


def p_ct_xas_final_basis(
    n_d_electrons: int,
    max_ligand_holes: int = 1,
) -> FockBasis:
    """Final `2p^5(3d^(n+1) + 3d^(n+2)L + ...)` basis."""
    return basis_from_configurations(
        xas_final_configurations(
            n_d_electrons,
            max_ligand_holes=max_ligand_holes,
        )
    )


def ligand_hole_count(state: int, n_d_spin_orbitals: int = 10) -> int:
    ligand_mask = ((1 << n_d_spin_orbitals) - 1) << n_d_spin_orbitals
    ligand_occupied = ((state & ligand_mask) >> n_d_spin_orbitals).bit_count()
    return n_d_spin_orbitals - ligand_occupied


def p_ct_ligand_hole_count(state: int) -> int:
    return ligand_hole_count(state >> 6)


def p_ct_core_hole_count(state: int) -> int:
    p_occupied = (state & ((1 << 6) - 1)).bit_count()
    return 6 - p_occupied


def charge_transfer_energy_matrix(
    basis: FockBasis,
    delta: float,
    n_d_spin_orbitals: int = 10,
    u_charge_transfer: float = 0.0,
) -> np.ndarray:
    """Diagonal configuration energy for ligand-hole sectors.

    `delta` is the `3d^(n+1)L - 3d^n` charge-transfer energy. The optional
    `u_charge_transfer` adds a simple quadratic penalty for multiple ligand
    holes: `h*delta + h*(h-1)/2*u_charge_transfer`.
    """
    matrix = np.zeros((len(basis), len(basis)), dtype=float)
    for idx, state in enumerate(basis.states):
        holes = ligand_hole_count(state, n_d_spin_orbitals)
        matrix[idx, idx] = holes * delta + holes * (holes - 1) * u_charge_transfer / 2.0
    return matrix


def p_ct_charge_transfer_energy_values(
    basis: FockBasis,
    delta: float,
    u_charge_transfer: float = 0.0,
    core_hole_potential: float = 0.0,
) -> np.ndarray:
    values = np.zeros(len(basis), dtype=float)
    for idx, state in enumerate(basis.states):
        holes = p_ct_ligand_hole_count(state)
        core_holes = p_ct_core_hole_count(state)
        if core_holes:
            values[idx] = holes * delta + holes * (holes + 1) * u_charge_transfer / 2.0
            values[idx] -= core_holes * holes * core_hole_potential
        else:
            values[idx] = holes * delta + holes * (holes - 1) * u_charge_transfer / 2.0
    return values


def d_ligand_hybridization_matrix(
    n_d_spin_orbitals: int = 10,
    hopping: float | np.ndarray = 1.0,
) -> np.ndarray:
    """One-body `d-L` hybridization matrix for a symmetry-adapted ligand shell."""
    n_orbitals = 2 * n_d_spin_orbitals
    h = np.zeros((n_orbitals, n_orbitals), dtype=np.result_type(hopping, complex))
    hopping_values = np.asarray(hopping)
    if hopping_values.ndim == 0:
        hopping_values = np.full(n_d_spin_orbitals, hopping_values.item())
    if hopping_values.shape != (n_d_spin_orbitals,):
        raise ValueError("hopping must be scalar or have length n_d_spin_orbitals")
    for idx, value in enumerate(hopping_values):
        d_idx = idx
        ligand_idx = n_d_spin_orbitals + idx
        h[d_idx, ligand_idx] = value
        h[ligand_idx, d_idx] = np.conjugate(value)
    return h


def d_ligand_hybridization_from_orbital_matrix(h_orbital: np.ndarray) -> np.ndarray:
    """Build a `d-L` spin-orbital hybridization matrix from a 5x5 d-orbital block."""

    from .shells import spin_expand_orbital_matrix

    h_spin = spin_expand_orbital_matrix(np.asarray(h_orbital))
    h = np.zeros((20, 20), dtype=np.result_type(h_spin, complex))
    h[0:10, 10:20] = h_spin
    h[10:20, 0:10] = h_spin.conjugate().T
    return h


def p_ct_one_body_matrix(h_p=None, h_d=None, h_ligand=None, h_hybridization=None) -> np.ndarray:
    """Pad p, d, ligand, and d-ligand one-body blocks into a 26x26 matrix."""
    h = np.zeros((26, 26), dtype=complex)
    if h_p is not None:
        h_p = np.asarray(h_p)
        if h_p.shape != (6, 6):
            raise ValueError("h_p must have shape (6,6)")
        h[0:6, 0:6] = h_p
    if h_d is not None:
        h_d = np.asarray(h_d)
        if h_d.shape != (10, 10):
            raise ValueError("h_d must have shape (10,10)")
        h[6:16, 6:16] = h_d
    if h_ligand is not None:
        h_ligand = np.asarray(h_ligand)
        if h_ligand.shape != (10, 10):
            raise ValueError("h_ligand must have shape (10,10)")
        h[16:26, 16:26] = h_ligand
    if h_hybridization is not None:
        h_hybridization = np.asarray(h_hybridization)
        if h_hybridization.shape != (20, 20):
            raise ValueError("h_hybridization must have shape (20,20)")
        h[6:26, 6:26] += h_hybridization
    return h


def p_ct_dipole_matrix(polarization: str = "z") -> np.ndarray:
    """Pad the p -> d L-edge dipole matrix into the p+d+ligand space."""
    from .shells import pd_l_edge_dipole_matrix

    pd = pd_l_edge_dipole_matrix(polarization)
    h = np.zeros((26, 26), dtype=np.result_type(pd, complex))
    h[0:16, 0:16] = pd
    return h


def pad_d_one_body_to_ct(h_d: np.ndarray, n_d_spin_orbitals: int = 10) -> np.ndarray:
    h_d = np.asarray(h_d)
    if h_d.shape != (n_d_spin_orbitals, n_d_spin_orbitals):
        raise ValueError("h_d must have shape (n_d_spin_orbitals, n_d_spin_orbitals)")
    h = np.zeros((2 * n_d_spin_orbitals, 2 * n_d_spin_orbitals), dtype=np.result_type(h_d, complex))
    h[:n_d_spin_orbitals, :n_d_spin_orbitals] = h_d
    return h


def charge_transfer_hamiltonian(
    basis: FockBasis,
    h_d: np.ndarray,
    delta: float,
    hopping: float | np.ndarray,
    u_charge_transfer: float = 0.0,
) -> np.ndarray:
    """Convenience Hamiltonian for the first charge-transfer cluster model."""
    h_one = pad_d_one_body_to_ct(h_d) + d_ligand_hybridization_matrix(hopping=hopping)
    return (
        one_body_matrix(basis, h_one)
        + charge_transfer_energy_matrix(basis, delta, u_charge_transfer=u_charge_transfer)
    )

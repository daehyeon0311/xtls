from __future__ import annotations

import numpy as np

from .basis import FockBasis
from .configuration import (
    basis_from_configurations,
    core_xps_final_configurations,
    sector_energy,
    valence_xps_final_configurations,
    xas_final_configurations,
    xas_initial_configurations,
)


# Orbital layout of the 26 spin-orbitals: 2p, then 3d, then the ligand shell.
P_SPIN_ORBITALS = tuple(range(0, 6))
D_SPIN_ORBITALS = tuple(range(6, 16))
LIGAND_SPIN_ORBITALS = tuple(range(16, 26))

PHOTOEMISSION_SHELLS = {
    "2p": P_SPIN_ORBITALS,
    "3d": D_SPIN_ORBITALS,
    "ligand": LIGAND_SPIN_ORBITALS,
}


def p_ct_spin_orbital_labels() -> tuple[str, ...]:
    """Return labels for the 26 spin-orbitals of the `2p + 3d + ligand` model.

    The ligand orbitals are the symmetry-adapted combinations that pair with
    each d orbital, so they carry matching labels.
    """
    from .shells import D_ORBITALS, P_ORBITALS, SPINS

    p_labels = tuple(f"p_{orb}_{spin}" for orb in P_ORBITALS for spin in SPINS)
    d_labels = tuple(f"d_{orb}_{spin}" for orb in D_ORBITALS for spin in SPINS)
    ligand_labels = tuple(f"L_{orb}_{spin}" for orb in D_ORBITALS for spin in SPINS)
    return p_labels + d_labels + ligand_labels


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


def p_ct_core_xps_final_basis(
    n_d_electrons: int,
    max_ligand_holes: int = 1,
) -> FockBasis:
    """Final `2p^5(3d^n + 3d^(n+1)L + ...)` basis for core-level XPS."""
    return basis_from_configurations(
        core_xps_final_configurations(
            n_d_electrons,
            max_ligand_holes=max_ligand_holes,
        )
    )


def p_ct_valence_xps_final_basis(
    n_d_electrons: int,
    max_ligand_holes: int = 1,
) -> FockBasis:
    """Final `2p^6(3d^(n-1) + 3d^n L + ...)` basis for valence-band XPS."""
    return basis_from_configurations(
        valence_xps_final_configurations(
            n_d_electrons,
            max_ligand_holes=max_ligand_holes,
        )
    )


def p_ct_photoemission_orbitals(shell: str = "2p") -> tuple[int, ...]:
    """Spin-orbital indices the photoelectron can be removed from.

    Photoemission intensity sums these channels incoherently, and that sum is
    invariant under a unitary change of basis inside the shell, so the cubic
    harmonics used here give the same spectrum as spherical ones would.
    """
    try:
        return PHOTOEMISSION_SHELLS[shell]
    except KeyError:
        raise ValueError(
            f"shell must be one of {sorted(PHOTOEMISSION_SHELLS)}, got {shell!r}"
        ) from None


def ligand_hole_count(state: int, n_d_spin_orbitals: int = 10) -> int:
    ligand_mask = ((1 << n_d_spin_orbitals) - 1) << n_d_spin_orbitals
    ligand_occupied = ((state & ligand_mask) >> n_d_spin_orbitals).bit_count()
    return n_d_spin_orbitals - ligand_occupied


def p_ct_ligand_hole_count(state: int) -> int:
    return ligand_hole_count(state >> 6)


def p_ct_core_hole_count(state: int) -> int:
    p_occupied = (state & ((1 << 6) - 1)).bit_count()
    return 6 - p_occupied


def p_ct_charge_transfer_energy_values(
    basis: FockBasis,
    delta: float,
    u_charge_transfer: float = 0.0,
    core_hole_potential: float = 0.0,
    d_electron_offset: int = 0,
) -> np.ndarray:
    """Diagonal sector energy for every state of a `2p + 3d + ligand` basis.

    The per-sector formula lives in `configuration.sector_energy`; this only
    classifies each state by its ligand-hole and core-hole count.
    `d_electron_offset` selects the spectroscopy: `0` for the initial state and
    for core-level XPS, `1` for XAS, `-1` for valence-band XPS.
    """
    values = np.zeros(len(basis), dtype=float)
    for idx, state in enumerate(basis.states):
        values[idx] = sector_energy(
            ligand_holes=p_ct_ligand_hole_count(state),
            core_holes=p_ct_core_hole_count(state),
            delta=delta,
            u_charge_transfer=u_charge_transfer,
            core_hole_potential=core_hole_potential,
            d_electron_offset=d_electron_offset,
        )
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



from __future__ import annotations

from math import factorial, sqrt

import numpy as np

from .basis import FockBasis


D_ORBITALS = ("xy", "yz", "zx", "x2-y2", "3z2-r2")
P_ORBITALS = ("x", "y", "z")
SPINS = ("down", "up")


def d_spin_orbital_labels() -> tuple[str, ...]:
    """Return labels for the 10 spin-orbitals of a d shell."""
    return tuple(f"{orb}_{spin}" for orb in D_ORBITALS for spin in SPINS)


def p_spin_orbital_labels() -> tuple[str, ...]:
    """Return labels for the 6 spin-orbitals of a p shell."""
    return tuple(f"{orb}_{spin}" for orb in P_ORBITALS for spin in SPINS)


def pd_spin_orbital_labels() -> tuple[str, ...]:
    """Return labels for a combined 2p + 3d shell orbital space."""
    return p_spin_orbital_labels() + d_spin_orbital_labels()


def d_shell_basis(n_d_electrons: int) -> FockBasis:
    """Return fixed-N basis for a spinful d shell."""
    return FockBasis.fixed_n(n_orbitals=10, n_electrons=n_d_electrons)


def pd_xas_initial_basis(n_d_electrons: int) -> FockBasis:
    """Return the `2p^6 3d^n` basis in a combined p+d orbital space."""
    if not 0 <= n_d_electrons <= 10:
        raise ValueError("n_d_electrons must satisfy 0 <= n <= 10")
    states = []
    p_full = (1 << 6) - 1
    for d_state in FockBasis.fixed_n(10, n_d_electrons).states:
        states.append(p_full | (d_state << 6))
    return FockBasis.from_states(16, 6 + n_d_electrons, states)


def pd_xas_final_basis(n_d_electrons: int) -> FockBasis:
    """Return the `2p^5 3d^(n+1)` basis in a combined p+d orbital space."""
    if not 0 <= n_d_electrons <= 9:
        raise ValueError("n_d_electrons must satisfy 0 <= n <= 9 for XAS final states")
    states = []
    for p_state in FockBasis.fixed_n(6, 5).states:
        for d_state in FockBasis.fixed_n(10, n_d_electrons + 1).states:
            states.append(p_state | (d_state << 6))
    return FockBasis.from_states(16, 6 + n_d_electrons, states)


def spin_expand_orbital_matrix(orbital_matrix: np.ndarray) -> np.ndarray:
    """Expand a 5x5 orbital matrix to a 10x10 spin-orbital matrix.

    The spin-orbital order is:
    `xy_down, xy_up, yz_down, yz_up, ...`.
    """
    orbital_matrix = np.asarray(orbital_matrix)
    if orbital_matrix.shape != (5, 5):
        raise ValueError("orbital_matrix must have shape (5, 5)")
    spin_matrix = np.zeros((10, 10), dtype=np.result_type(orbital_matrix, complex))
    for orbital_i in range(5):
        for orbital_j in range(5):
            for spin in range(2):
                spin_matrix[2 * orbital_i + spin, 2 * orbital_j + spin] = orbital_matrix[
                    orbital_i, orbital_j
                ]
    return _real_if_close(spin_matrix)


def d_angular_momentum_matrices() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Lx, Ly, Lz matrices in the cubic d-orbital basis.

    The cubic basis order is `xy, yz, zx, x2-y2, 3z2-r2`. The spherical basis
    order used internally is m = -2, -1, 0, 1, 2.
    """
    l_quantum = 2
    m_values = np.arange(-l_quantum, l_quantum + 1)
    lz_sph = np.diag(m_values.astype(float)).astype(complex)
    lp_sph = np.zeros((5, 5), dtype=complex)
    for col, m in enumerate(m_values[:-1]):
        row = col + 1
        lp_sph[row, col] = np.sqrt(l_quantum * (l_quantum + 1) - m * (m + 1))
    lm_sph = lp_sph.conjugate().T
    lx_sph = (lp_sph + lm_sph) / 2
    ly_sph = (lp_sph - lm_sph) / (2j)

    transform = _spherical_to_cubic_transform()
    lx = transform @ lx_sph @ transform.conjugate().T
    ly = transform @ ly_sph @ transform.conjugate().T
    lz = transform @ lz_sph @ transform.conjugate().T
    return lx, ly, lz


def p_angular_momentum_matrices() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Lx, Ly, Lz matrices in the cubic p-orbital basis `x,y,z`."""
    l_quantum = 1
    m_values = np.arange(-l_quantum, l_quantum + 1)
    lz_sph = np.diag(m_values.astype(float)).astype(complex)
    lp_sph = np.zeros((3, 3), dtype=complex)
    for col, m in enumerate(m_values[:-1]):
        row = col + 1
        lp_sph[row, col] = np.sqrt(l_quantum * (l_quantum + 1) - m * (m + 1))
    lm_sph = lp_sph.conjugate().T
    lx_sph = (lp_sph + lm_sph) / 2
    ly_sph = (lp_sph - lm_sph) / (2j)

    transform = _spherical_to_p_cubic_transform()
    lx = transform @ lx_sph @ transform.conjugate().T
    ly = transform @ ly_sph @ transform.conjugate().T
    lz = transform @ lz_sph @ transform.conjugate().T
    return lx, ly, lz


def spin_matrices() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Sx, Sy, Sz for spin order `down, up`."""
    sx = np.array([[0, 0.5], [0.5, 0]], dtype=complex)
    sy = np.array([[0, 0.5j], [-0.5j, 0]], dtype=complex)
    sz = np.array([[-0.5, 0], [0, 0.5]], dtype=complex)
    return sx, sy, sz


def d_spin_orbit_matrix(zeta: float = 1.0) -> np.ndarray:
    """Return the 10x10 one-body matrix for zeta * L dot S."""
    lx, ly, lz = d_angular_momentum_matrices()
    sx, sy, sz = spin_matrices()
    matrix = zeta * (np.kron(lx, sx) + np.kron(ly, sy) + np.kron(lz, sz))
    return _real_if_close(matrix)


def p_spin_orbit_matrix(zeta: float = 1.0) -> np.ndarray:
    """Return the 6x6 one-body matrix for zeta * L dot S in a p shell."""
    lx, ly, lz = p_angular_momentum_matrices()
    sx, sy, sz = spin_matrices()
    matrix = zeta * (np.kron(lx, sx) + np.kron(ly, sy) + np.kron(lz, sz))
    return _real_if_close(matrix)


def pd_l_edge_dipole_matrix(polarization: str = "z") -> np.ndarray:
    """Return a 16x16 one-body matrix for the `2p -> 3d` dipole operator.

    Rows are created orbitals and columns are annihilated orbitals. The combined
    p+d order is `p_x,p_y,p_z,d_xy,d_yz,d_zx,d_x2-y2,d_3z2-r2`, each with
    `down, up` spin channels.
    """
    orbital = np.zeros((8, 8), dtype=complex)
    p_to_d = _p_to_d_orbital_dipole_matrix(polarization)
    orbital[3:8, 0:3] = p_to_d

    spin_matrix = np.zeros((16, 16), dtype=complex)
    for orbital_i in range(8):
        for orbital_j in range(8):
            for spin in range(2):
                spin_matrix[2 * orbital_i + spin, 2 * orbital_j + spin] = orbital[
                    orbital_i,
                    orbital_j,
                ]
    return _real_if_close(spin_matrix)


def d_shell_kanamori_tensor(u: float, j_hund: float, n_orbitals: int = 5) -> np.ndarray:
    """Return a spin-orbital Coulomb tensor for a d-shell Kanamori model.

    The spin-orbital order is the same as `d_spin_orbital_labels`: each orbital
    has `down, up` spin channels. The tensor can be passed directly to
    `two_body_matrix`.

    Parameters
    ----------
    u:
        Intra-orbital Coulomb repulsion.
    j_hund:
        Hund exchange. The inter-orbital direct term is `u - 2*j_hund`.
    n_orbitals:
        Number of spatial orbitals. The d shell uses 5.
    """
    if n_orbitals <= 0:
        raise ValueError("n_orbitals must be positive")
    u_prime = u - 2.0 * j_hund
    orbital_tensor = np.zeros((n_orbitals, n_orbitals, n_orbitals, n_orbitals), dtype=float)

    for a in range(n_orbitals):
        orbital_tensor[a, a, a, a] = u
        for b in range(n_orbitals):
            if a == b:
                continue
            orbital_tensor[a, b, a, b] = u_prime
            orbital_tensor[a, b, b, a] = j_hund

    n_spin_orbitals = 2 * n_orbitals
    tensor = np.zeros(
        (n_spin_orbitals, n_spin_orbitals, n_spin_orbitals, n_spin_orbitals),
        dtype=float,
    )
    for p in range(n_spin_orbitals):
        a, spin_p = divmod(p, 2)
        for q in range(n_spin_orbitals):
            b, spin_q = divmod(q, 2)
            for r in range(n_spin_orbitals):
                c, spin_r = divmod(r, 2)
                if spin_p != spin_r:
                    continue
                for s in range(n_spin_orbitals):
                    d, spin_s = divmod(s, 2)
                    if spin_q == spin_s:
                        tensor[p, q, r, s] = orbital_tensor[a, b, c, d]
    return tensor


def d_shell_slater_to_kanamori(f0: float, f2: float, f4: float) -> tuple[float, float]:
    """Convert d-shell Slater integrals to average Kanamori U and J.

    This is the common spherical-average reduction:
    `U = F0 + 4/49 * (F2 + F4)` and `J = (F2 + F4) / 14`.
    It is useful for a compact interacting solver while the full anisotropic
    Slater tensor is being built.
    """
    u = f0 + 4.0 * (f2 + f4) / 49.0
    j_hund = (f2 + f4) / 14.0
    return u, j_hund


def d_shell_slater_kanamori_tensor(f0: float, f2: float, f4: float) -> np.ndarray:
    """Return a d-shell Kanamori tensor parameterized by Slater integrals."""
    u, j_hund = d_shell_slater_to_kanamori(f0, f2, f4)
    return d_shell_kanamori_tensor(u, j_hund)


def slater_spherical_orbital_tensor(
    l_quantum: int,
    slater_integrals: dict[int, float],
) -> np.ndarray:
    """Return the anisotropic Coulomb tensor in a spherical-harmonic shell.

    The output tensor is `<m1,m2|V|m3,m4>` with m ordered from `-l` to `+l`.
    `slater_integrals` maps even k values to Slater integrals, for example
    `{0: F0, 2: F2, 4: F4}` for a d shell.
    """
    if l_quantum < 0:
        raise ValueError("l_quantum must be non-negative")
    n_orbitals = 2 * l_quantum + 1
    m_values = list(range(-l_quantum, l_quantum + 1))
    tensor = np.zeros((n_orbitals, n_orbitals, n_orbitals, n_orbitals), dtype=complex)

    for k, f_k in slater_integrals.items():
        if k < 0 or k > 2 * l_quantum or k % 2:
            raise ValueError("Slater k values must be even and satisfy 0 <= k <= 2*l")
        if f_k == 0:
            continue
        for a, m1 in enumerate(m_values):
            for b, m2 in enumerate(m_values):
                for c, m3 in enumerate(m_values):
                    for d, m4 in enumerate(m_values):
                        value = 0.0
                        for q in range(-k, k + 1):
                            left = _reduced_spherical_harmonic_matrix(l_quantum, k, q, m1, m3)
                            right = ((-1) ** q) * _reduced_spherical_harmonic_matrix(
                                l_quantum,
                                k,
                                -q,
                                m2,
                                m4,
                            )
                            value += left * right
                        tensor[a, b, c, d] += f_k * value
    return _real_if_close(tensor)


def orbital_tensor_to_spin_orbital(orbital_tensor: np.ndarray) -> np.ndarray:
    """Expand an orbital Coulomb tensor into spin-orbital form."""
    orbital_tensor = np.asarray(orbital_tensor)
    if orbital_tensor.ndim != 4 or len(set(orbital_tensor.shape)) != 1:
        raise ValueError("orbital_tensor must have shape (n,n,n,n)")
    n_orbitals = orbital_tensor.shape[0]
    n_spin_orbitals = 2 * n_orbitals
    tensor = np.zeros(
        (n_spin_orbitals, n_spin_orbitals, n_spin_orbitals, n_spin_orbitals),
        dtype=np.result_type(orbital_tensor, complex),
    )
    for p in range(n_spin_orbitals):
        a, spin_p = divmod(p, 2)
        for q in range(n_spin_orbitals):
            b, spin_q = divmod(q, 2)
            for r in range(n_spin_orbitals):
                c, spin_r = divmod(r, 2)
                if spin_p != spin_r:
                    continue
                for s in range(n_spin_orbitals):
                    d, spin_s = divmod(s, 2)
                    if spin_q == spin_s:
                        tensor[p, q, r, s] = orbital_tensor[a, b, c, d]
    return _real_if_close(tensor)


def transform_orbital_tensor(orbital_tensor: np.ndarray, transform: np.ndarray) -> np.ndarray:
    """Transform `<ab|V|cd>` with `|new> = transform |old>`."""
    orbital_tensor = np.asarray(orbital_tensor)
    transform = np.asarray(transform)
    if orbital_tensor.ndim != 4 or len(set(orbital_tensor.shape)) != 1:
        raise ValueError("orbital_tensor must have shape (n,n,n,n)")
    n_orbitals = orbital_tensor.shape[0]
    if transform.shape != (n_orbitals, n_orbitals):
        raise ValueError("transform must have shape (n,n)")
    transformed = np.einsum(
        "am,bn,co,dp,mnop->abcd",
        transform.conjugate(),
        transform.conjugate(),
        transform,
        transform,
        orbital_tensor,
        optimize=True,
    )
    return _real_if_close(transformed)


def d_shell_slater_tensor(f0: float, f2: float, f4: float, cubic: bool = True) -> np.ndarray:
    """Return the full anisotropic d-shell spin-orbital Slater tensor."""
    orbital_tensor = slater_spherical_orbital_tensor(2, {0: f0, 2: f2, 4: f4})
    if cubic:
        orbital_tensor = transform_orbital_tensor(orbital_tensor, _spherical_to_cubic_transform())
    return orbital_tensor_to_spin_orbital(orbital_tensor)


def pd_shell_slater_tensor(
    fdd2: float,
    fdd4: float,
    fpd2: float,
    gpd1: float,
    gpd3: float,
    fdd0: float = 0.0,
    fpd0: float = 0.0,
) -> np.ndarray:
    """Return a combined 2p+3d spin-orbital Coulomb tensor.

    The d-d block uses `Fdd0/Fdd2/Fdd4`. The p-d direct block uses
    `Fpd0/Fpd2`; the p-d exchange block uses `Gpd1/Gpd3`. Appendix A provides
    the multipole values `Fdd2`, `Fdd4`, `Fpd2`, `Gpd1`, and `Gpd3`. The
    monopole values are kept as explicit optional parameters because they are
    strongly screened and usually fitted separately.
    """
    spherical = np.zeros((8, 8, 8, 8), dtype=complex)
    spherical[3:8, 3:8, 3:8, 3:8] = slater_spherical_orbital_tensor(
        2,
        {0: fdd0, 2: fdd2, 4: fdd4},
    )
    _add_pd_direct_exchange_spherical(spherical, fpd0, fpd2, gpd1, gpd3)

    transform = np.zeros((8, 8), dtype=complex)
    transform[0:3, 0:3] = _spherical_to_p_cubic_transform()
    transform[3:8, 3:8] = _spherical_to_cubic_transform()
    cubic = transform_orbital_tensor(spherical, transform)
    return orbital_tensor_to_spin_orbital(cubic)


def _add_pd_direct_exchange_spherical(
    tensor: np.ndarray,
    fpd0: float,
    fpd2: float,
    gpd1: float,
    gpd3: float,
) -> None:
    p_m = list(range(-1, 2))
    d_m = list(range(-2, 3))
    direct = _mixed_direct_tensor(2, 1, {0: fpd0, 2: fpd2})
    exchange = _mixed_exchange_tensor(2, 1, {1: gpd1, 3: gpd3})

    for a, _md1 in enumerate(d_m):
        d1 = 3 + a
        for b, _mp1 in enumerate(p_m):
            p1 = b
            for c, _md2 in enumerate(d_m):
                d2 = 3 + c
                for e, _mp2 in enumerate(p_m):
                    p2 = e
                    direct_value = direct[a, b, c, e]
                    exchange_value = exchange[a, b, e, c]
                    tensor[d1, p1, d2, p2] += direct_value
                    tensor[p1, d1, p2, d2] += direct_value
                    tensor[d1, p1, p2, d2] += exchange_value
                    tensor[p1, d1, d2, p2] += np.conjugate(exchange_value)


def _mixed_direct_tensor(
    l_left: int,
    l_right: int,
    slater_integrals: dict[int, float],
) -> np.ndarray:
    left_m = list(range(-l_left, l_left + 1))
    right_m = list(range(-l_right, l_right + 1))
    tensor = np.zeros((len(left_m), len(right_m), len(left_m), len(right_m)), dtype=complex)
    for k, value_k in slater_integrals.items():
        if value_k == 0:
            continue
        for a, m1 in enumerate(left_m):
            for b, m2 in enumerate(right_m):
                for c, m3 in enumerate(left_m):
                    for d, m4 in enumerate(right_m):
                        value = 0.0
                        for q in range(-k, k + 1):
                            value += _reduced_spherical_harmonic_between(
                                l_left,
                                l_left,
                                k,
                                q,
                                m1,
                                m3,
                            ) * ((-1) ** q) * _reduced_spherical_harmonic_between(
                                l_right,
                                l_right,
                                k,
                                -q,
                                m2,
                                m4,
                            )
                        tensor[a, b, c, d] += value_k * value
    return tensor


def _mixed_exchange_tensor(
    l_left: int,
    l_right: int,
    exchange_integrals: dict[int, float],
) -> np.ndarray:
    left_m = list(range(-l_left, l_left + 1))
    right_m = list(range(-l_right, l_right + 1))
    tensor = np.zeros((len(left_m), len(right_m), len(right_m), len(left_m)), dtype=complex)
    for k, value_k in exchange_integrals.items():
        if value_k == 0:
            continue
        for a, m1 in enumerate(left_m):
            for b, m2 in enumerate(right_m):
                for c, m3 in enumerate(right_m):
                    for d, m4 in enumerate(left_m):
                        value = 0.0
                        for q in range(-k, k + 1):
                            value += _reduced_spherical_harmonic_between(
                                l_left,
                                l_right,
                                k,
                                q,
                                m1,
                                m3,
                            ) * ((-1) ** q) * _reduced_spherical_harmonic_between(
                                l_right,
                                l_left,
                                k,
                                -q,
                                m2,
                                m4,
                            )
                        tensor[a, b, c, d] += value_k * value
    return tensor


def _reduced_spherical_harmonic_matrix(
    l_quantum: int,
    k: int,
    q: int,
    m_bra: int,
    m_ket: int,
) -> float:
    return _reduced_spherical_harmonic_between(
        l_quantum,
        l_quantum,
        k,
        q,
        m_bra,
        m_ket,
    )


def _reduced_spherical_harmonic_between(
    l_bra: int,
    l_ket: int,
    k: int,
    q: int,
    m_bra: int,
    m_ket: int,
) -> float:
    return (
        ((-1) ** m_bra)
        * sqrt((2 * l_bra + 1) * (2 * l_ket + 1))
        * _wigner_3j(l_bra, k, l_ket, 0, 0, 0)
        * _wigner_3j(l_bra, k, l_ket, -m_bra, q, m_ket)
    )


def _p_to_d_orbital_dipole_matrix(polarization: str) -> np.ndarray:
    pol = polarization.lower()
    if pol not in {"x", "y", "z"}:
        raise ValueError("polarization must be 'x', 'y', or 'z'")

    d_m = list(range(-2, 3))
    p_m = list(range(-1, 2))
    spherical_q: dict[int, complex]
    if pol == "x":
        spherical_q = {-1: 1 / np.sqrt(2), 1: -1 / np.sqrt(2)}
    elif pol == "y":
        spherical_q = {-1: 1j / np.sqrt(2), 1: 1j / np.sqrt(2)}
    else:
        spherical_q = {0: 1.0}

    spherical = np.zeros((5, 3), dtype=complex)
    for row, md in enumerate(d_m):
        for col, mp in enumerate(p_m):
            for q, coefficient in spherical_q.items():
                spherical[row, col] += coefficient * _reduced_spherical_harmonic_between(
                    2,
                    1,
                    1,
                    q,
                    md,
                    mp,
                )

    d_transform = _spherical_to_cubic_transform()
    p_transform = _spherical_to_p_cubic_transform()
    cubic = d_transform @ spherical @ p_transform.conjugate().T
    return _real_if_close(cubic)


def _wigner_3j(
    j1: int,
    j2: int,
    j3: int,
    m1: int,
    m2: int,
    m3: int,
) -> float:
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0

    delta = (
        factorial(j1 + j2 - j3)
        * factorial(j1 - j2 + j3)
        * factorial(-j1 + j2 + j3)
        / factorial(j1 + j2 + j3 + 1)
    )
    norm = sqrt(
        delta
        * factorial(j1 + m1)
        * factorial(j1 - m1)
        * factorial(j2 + m2)
        * factorial(j2 - m2)
        * factorial(j3 + m3)
        * factorial(j3 - m3)
    )
    phase = (-1) ** (j1 - j2 - m3)

    z_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
    z_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    total = 0.0
    for z in range(z_min, z_max + 1):
        denom = (
            factorial(z)
            * factorial(j1 + j2 - j3 - z)
            * factorial(j1 - m1 - z)
            * factorial(j2 + m2 - z)
            * factorial(j3 - j2 + m1 + z)
            * factorial(j3 - j1 - m2 + z)
        )
        total += ((-1) ** z) / denom
    return phase * norm * total


def _spherical_to_cubic_transform() -> np.ndarray:
    return np.array(
        [
            [1j / np.sqrt(2), 0, 0, 0, -1j / np.sqrt(2)],
            [0, 1j / np.sqrt(2), 0, 1j / np.sqrt(2), 0],
            [0, 1 / np.sqrt(2), 0, -1 / np.sqrt(2), 0],
            [1 / np.sqrt(2), 0, 0, 0, 1 / np.sqrt(2)],
            [0, 0, 1, 0, 0],
        ],
        dtype=complex,
    )


def _spherical_to_p_cubic_transform() -> np.ndarray:
    return np.array(
        [
            [1 / np.sqrt(2), 0, -1 / np.sqrt(2)],
            [1j / np.sqrt(2), 0, 1j / np.sqrt(2)],
            [0, 1, 0],
        ],
        dtype=complex,
    )


def _real_if_close(matrix: np.ndarray) -> np.ndarray:
    if np.allclose(matrix.imag, 0):
        return matrix.real
    return matrix

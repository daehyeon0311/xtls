from __future__ import annotations

import numpy as np

from .engine.shells import spherical_to_cubic_transform


def crystal_field(
    positions_sph: np.ndarray,
    ten_dq: float,
    r2: float,
    r4: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return crystal-field matrices and Tanaka Qkm coefficients.

    `positions_sph` rows are `[radius, theta, phi]` in radians.
    """
    qkm = np.zeros((3, 5), dtype=complex)

    for radius, theta, phi in positions_sph:
        thkm = np.zeros((3, 5), dtype=float)
        thkm[0, 0] = 1 / np.sqrt(2)
        thkm[1, 2] = np.sqrt(15) / 4 * np.sin(theta) ** 2
        thkm[1, 1] = -np.sqrt(15) / 2 * np.cos(theta) * np.sin(theta)
        thkm[1, 0] = np.sqrt(5) / (np.sqrt(2) * 2) * (
            2 * np.cos(theta) ** 2 - np.sin(theta) ** 2
        )
        thkm[2, 4] = 3 * np.sqrt(35) / 16 * np.sin(theta) ** 4
        thkm[2, 3] = (
            -3 * np.sqrt(35) / (4 * np.sqrt(2)) * np.cos(theta) * np.sin(theta) ** 3
        )
        thkm[2, 2] = (
            3
            * np.sqrt(5)
            / 8
            * np.sin(theta) ** 2
            * (6 * np.cos(theta) ** 2 - np.sin(theta) ** 2)
        )
        thkm[2, 1] = (
            -(3 * np.sqrt(5))
            / (4 * np.sqrt(2))
            * np.cos(theta)
            * np.sin(theta)
            * (4 * np.cos(theta) ** 2 - 3 * np.sin(theta) ** 2)
        )
        thkm[2, 0] = (
            3
            / (8 * np.sqrt(2))
            * (
                8 * np.cos(theta) ** 4
                - 24 * np.cos(theta) ** 2 * np.sin(theta) ** 2
                + 3 * np.sin(theta) ** 4
            )
        )

        ykm = np.zeros((3, 5), dtype=complex)
        for m in range(5):
            ykm[:, m] = (
                1
                / np.sqrt(2 * np.pi)
                * thkm[:, m]
                * np.exp(1j * m * phi)
            )

        qkm[0, :] += np.sqrt(4 * np.pi / 1) / radius * ykm[0, :]
        qkm[1, :] += np.sqrt(4 * np.pi / 5) / radius**3 * ykm[1, :]
        qkm[2, :] += np.sqrt(4 * np.pi / 9) / radius**5 * ykm[2, :]

    q0 = np.zeros((5, 5), dtype=complex)
    q2 = np.array(
        [
            [qkm[1, 0], -qkm[1, 1].conjugate(), qkm[1, 2].conjugate(), 0, 0],
            [qkm[1, 1], qkm[1, 0], -qkm[1, 1].conjugate(), qkm[1, 2].conjugate(), 0],
            [
                qkm[1, 2],
                qkm[1, 1],
                qkm[1, 0],
                -qkm[1, 1].conjugate(),
                qkm[1, 2].conjugate(),
            ],
            [0, qkm[1, 2], qkm[1, 1], qkm[1, 0], -qkm[1, 1].conjugate()],
            [0, 0, qkm[1, 2], qkm[1, 1], qkm[1, 0]],
        ],
        dtype=complex,
    )
    q4 = np.array(
        [
            [
                qkm[2, 0],
                -qkm[2, 1].conjugate(),
                qkm[2, 2].conjugate(),
                -qkm[2, 3].conjugate(),
                qkm[2, 4].conjugate(),
            ],
            [
                qkm[2, 1],
                qkm[2, 0],
                -qkm[2, 1].conjugate(),
                qkm[2, 2].conjugate(),
                -qkm[2, 3].conjugate(),
            ],
            [
                qkm[2, 2],
                qkm[2, 1],
                qkm[2, 0],
                -qkm[2, 1].conjugate(),
                qkm[2, 2].conjugate(),
            ],
            [qkm[2, 3], qkm[2, 2], qkm[2, 1], qkm[2, 0], -qkm[2, 1].conjugate()],
            [qkm[2, 4], qkm[2, 3], qkm[2, 2], qkm[2, 1], qkm[2, 0]],
        ],
        dtype=complex,
    )

    c0 = np.eye(5)
    c2 = np.array(
        [
            [-2 / 7, np.sqrt(6) / 7, -2 / 7, 0, 0],
            [-np.sqrt(6) / 7, 1 / 7, 1 / 7, -np.sqrt(6) / 7, 0],
            [-2 / 7, -1 / 7, 2 / 7, -1 / 7, -2 / 7],
            [0, -np.sqrt(6) / 7, 1 / 7, 1 / 7, -np.sqrt(6) / 7],
            [0, 0, -2 / 7, np.sqrt(6) / 7, -2 / 7],
        ],
        dtype=complex,
    )
    c4 = np.array(
        [
            [1, -np.sqrt(5), np.sqrt(15), -np.sqrt(35), np.sqrt(70)],
            [np.sqrt(5), -4, np.sqrt(30), -2 * np.sqrt(10), np.sqrt(35)],
            [np.sqrt(15), -np.sqrt(30), 6, -np.sqrt(30), np.sqrt(15)],
            [np.sqrt(35), -2 * np.sqrt(10), np.sqrt(30), -4, np.sqrt(5)],
            [np.sqrt(70), -np.sqrt(35), np.sqrt(15), -np.sqrt(5), 1],
        ],
        dtype=complex,
    ) / 21

    h_cry = q0 * c0 + r2 * q2 * c2 + r4 * q4 * c4
    e_cry = np.sort(np.real(np.linalg.eigvals(h_cry)))
    positive_count = int(np.sum(e_cry > 0))
    e_plus = np.sum(e_cry[5 - positive_count : 5])
    e_minus = np.sum(e_cry[: 5 - positive_count])
    dq = e_plus / positive_count - e_minus / (5 - positive_count)
    norm = ten_dq / dq
    h_cry = norm * h_cry
    e_cry = norm * e_cry

    transform = spherical_to_cubic_transform()
    h_cry_xyz = transform @ h_cry @ transform.conjugate().T
    qkm_tanaka = norm * np.vstack([qkm[0, :] * 1, qkm[1, :] * r2, qkm[2, :] * r4])
    return h_cry_xyz, h_cry, e_cry, qkm_tanaka


def hybridization(
    positions_sph: np.ndarray,
    pd_sigma: float,
    ratio: float,
    d_ref: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ligand hybridization matrices and compact Tanaka Bkm coefficients."""
    positions_xyz = np.column_stack(
        [
            positions_sph[:, 0] * np.sin(positions_sph[:, 1]) * np.cos(positions_sph[:, 2]),
            positions_sph[:, 0] * np.sin(positions_sph[:, 1]) * np.sin(positions_sph[:, 2]),
            positions_sph[:, 0] * np.cos(positions_sph[:, 1]),
        ]
    )
    pd_pi = pd_sigma / ratio
    pds = pd_sigma * (positions_sph[:, 0] / d_ref) ** (-3.5)
    pdp = pd_pi * (positions_sph[:, 0] / d_ref) ** (-3.5)
    l = positions_xyz[:, 0] / positions_sph[:, 0]
    m = positions_xyz[:, 1] / positions_sph[:, 0]
    n = positions_xyz[:, 2] / positions_sph[:, 0]

    ligand_count = positions_sph.shape[0]
    h_hyb = np.zeros((3 * ligand_count, 5), dtype=float)
    for idx in range(ligand_count):
        row = 3 * idx
        li, mi, ni = l[idx], m[idx], n[idx]
        sig, pi = pds[idx], pdp[idx]

        h_hyb[row, 0] = np.sqrt(3) * li**2 * mi * sig + mi * (1 - 2 * li**2) * pi
        h_hyb[row, 1] = np.sqrt(3) * li * mi * ni * sig - 2 * li * mi * ni * pi
        h_hyb[row, 2] = np.sqrt(3) * li**2 * ni * sig + ni * (1 - 2 * li**2) * pi
        h_hyb[row, 3] = (
            0.5 * np.sqrt(3) * li * (li**2 - mi**2) * sig
            + li * (1 - li**2 + mi**2) * pi
        )
        h_hyb[row, 4] = (
            li * (ni**2 - 0.5 * (li**2 + mi**2)) * sig
            - np.sqrt(3) * li * ni**2 * pi
        )

        h_hyb[row + 1, 0] = np.sqrt(3) * mi**2 * li * sig + li * (1 - 2 * mi**2) * pi
        h_hyb[row + 1, 1] = np.sqrt(3) * mi**2 * ni * sig + ni * (1 - 2 * mi**2) * pi
        h_hyb[row + 1, 2] = np.sqrt(3) * mi * ni * li * sig - 2 * mi * ni * li * pi
        h_hyb[row + 1, 3] = (
            0.5 * np.sqrt(3) * mi * (li**2 - mi**2) * sig
            - mi * (1 + li**2 - mi**2) * pi
        )
        h_hyb[row + 1, 4] = (
            mi * (ni**2 - 0.5 * (li**2 + mi**2)) * sig
            - np.sqrt(3) * mi * ni**2 * pi
        )

        h_hyb[row + 2, 0] = np.sqrt(3) * ni * li * mi * sig - 2 * ni * li * mi * pi
        h_hyb[row + 2, 1] = np.sqrt(3) * ni**2 * mi * sig + mi * (1 - 2 * ni**2) * pi
        h_hyb[row + 2, 2] = np.sqrt(3) * ni**2 * li * sig + li * (1 - 2 * ni**2) * pi
        h_hyb[row + 2, 3] = (
            0.5 * np.sqrt(3) * ni * (li**2 - mi**2) * sig
            - ni * (li**2 - mi**2) * pi
        )
        h_hyb[row + 2, 4] = (
            ni * (ni**2 - 0.5 * (li**2 + mi**2)) * sig
            + np.sqrt(3) * ni * (li**2 + mi**2) * pi
        )

    h_hyb_ld = np.diag(np.sqrt(np.sum(h_hyb**2, axis=0)))
    e_hyb = np.linalg.eigvals(h_hyb_ld)
    bkm = _compact_bkm(h_hyb_ld)
    return h_hyb_ld, h_hyb, e_hyb, bkm


def _compact_bkm(h: np.ndarray) -> np.ndarray:
    b = np.zeros((3, 5), dtype=complex)
    b[0, 0] = (h[3, 3] + h[0, 0] + h[1, 1] + h[2, 2] + h[4, 4]) / 5
    b[1, 0] = 0.5 * (-2 * h[3, 3] - 2 * h[0, 0] + h[2, 2] + h[1, 1] + 2 * h[4, 4])
    b[1, 1] = 0.5 * (
        np.sqrt(6)
        * (1j * h[0, 2] - h[0, 1] - h[2, 3] - 1j * h[1, 3] + 1j * h[1, 4])
        - np.sqrt(2) * h[2, 4]
    )
    b[1, 2] = 0.25 * (
        4 * np.sqrt(2) * (1j * h[0, 4] - h[3, 4])
        + np.sqrt(6) * (h[2, 2] - h[1, 1] - 2j * h[1, 2])
    )
    b[2, 0] = 0.3 * (h[3, 3] + h[0, 0] - 4 * h[2, 2] - 4 * h[1, 1] + 6 * h[4, 4])
    b[2, 1] = 0.3 * np.sqrt(5) * (
        -1j * h[0, 2]
        + h[0, 1]
        + h[2, 3]
        - 2 * np.sqrt(3) * h[2, 4]
        + 1j * h[1, 3]
        + 2 * np.sqrt(3) * 1j * h[1, 4]
    )
    b[2, 2] = 0.3 * np.sqrt(10) * (
        np.sqrt(3) * h[3, 4]
        - 1j * np.sqrt(3) * h[0, 4]
        + h[2, 2]
        - h[1, 1]
        - 2j * h[1, 2]
    )
    b[2, 3] = 0.3 * np.sqrt(35) * (
        1j * h[0, 2] + h[0, 1] - h[2, 3] + 1j * h[1, 3]
    )
    b[2, 4] = 3 / 20 * np.sqrt(70) * (h[3, 3] - h[0, 0] - 2j * h[0, 3])
    return b

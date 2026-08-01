"""Ba2FeSi2O7 unit cell, and Quantum ESPRESSO inputs for the exchange coupling.

    python dft/build_structure.py

Structure from the Rietveld refinement of Jang et al., Phys. Rev. B 104,
214434 (2021), Table II (T = 1.7 K). Note their site labels for O2 and O3 are
swapped relative to the coordinate forms: (0.3651, 0.8651, z) has y = x + 1/2
and so sits on 4e, while (0.0769, 0.1984, z) is a general 8f position. The
coordinates themselves are consistent, and reproduce the FeO4 tetrahedron to
0.4% in bond length.

J is extracted from the energy difference between ferromagnetic and Neel
alignments of the two Fe sites in the cell. With H = J sum_<ij> S_i . S_j over
the four nearest-neighbour bonds per cell,

    E_FM - E_AFM = 8 J S^2   ->   J = (E_FM - E_AFM) / (8 S^2),  S = 2
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent

LATTICE_A = 8.3194
LATTICE_C = 5.3336

# Wyckoff representatives, Table II of the reference.
REPRESENTATIVES = {
    "Ba": (0.1644, 0.6644, 0.5098),
    "Fe": (0.0000, 0.0000, 0.0000),
    "Si": (0.3645, 0.8645, 0.9609),
    "O1": (0.0000, 0.5000, 0.1383),
    "O2": (0.3651, 0.8651, 0.2642),
    "O3": (0.0769, 0.1984, 0.1694),
}

MASSES = {"Ba": 137.327, "Fe": 55.845, "Si": 28.0855, "O": 15.999}


def general_positions(x: float, y: float, z: float) -> list[tuple[float, float, float]]:
    """The eight general positions of P-42_1m (space group 113)."""
    return [
        (x, y, z),
        (-x, -y, z),
        (y, -x, -z),
        (-y, x, -z),
        (-x + 0.5, y + 0.5, -z),
        (x + 0.5, -y + 0.5, -z),
        (-y + 0.5, -x + 0.5, z),
        (y + 0.5, x + 0.5, z),
    ]


def orbit(x: float, y: float, z: float) -> list[np.ndarray]:
    """Distinct sites generated from one representative."""
    found: list[np.ndarray] = []
    for position in general_positions(x, y, z):
        candidate = np.mod(np.asarray(position, dtype=float), 1.0)
        if not any(
            np.allclose(np.mod(candidate - existing + 0.5, 1.0) - 0.5, 0.0, atol=1e-4)
            for existing in found
        ):
            found.append(candidate)
    return found


def build() -> dict[str, list[np.ndarray]]:
    return {label: orbit(*rep) for label, rep in REPRESENTATIVES.items()}


def coordination(sites: dict[str, list[np.ndarray]], cutoff: float = 2.4):
    """Oxygen neighbours of the Fe at the origin, as (distance, label, vector)."""
    cell = np.diag([LATTICE_A, LATTICE_A, LATTICE_C])
    iron = sites["Fe"][0]
    neighbours = []
    for label in ("O1", "O2", "O3"):
        for site in sites[label]:
            for shift in itertools.product([-1, 0, 1], repeat=3):
                vector = (site + np.array(shift) - iron) @ cell
                distance = float(np.linalg.norm(vector))
                if distance < cutoff:
                    neighbours.append((distance, label, vector))
    neighbours.sort(key=lambda entry: entry[0])
    return neighbours


def write_pw_input(
    path: Path,
    sites: dict[str, list[np.ndarray]],
    *,
    antiferromagnetic: bool,
    hubbard_u: float = 4.0,
) -> None:
    """One scf input, with the two Fe sites either parallel or antiparallel."""
    tag = ("afm" if antiferromagnetic else "fm") + f"_u{hubbard_u:g}".replace(".", "p")
    lines = [
        "&CONTROL",
        "  calculation = 'scf'",
        f"  prefix = '{tag}'",
        "  outdir = './tmp'",
        "  pseudo_dir = './pseudo'",
        "  tprnfor = .true.",
        "/",
        "&SYSTEM",
        "  ibrav = 6",
        f"  celldm(1) = {LATTICE_A / 0.529177210903:.8f}",
        f"  celldm(3) = {LATTICE_C / LATTICE_A:.8f}",
        "  nat = 24",
        "  ntyp = %d" % (5 if antiferromagnetic else 4),
        "  ecutwfc = 60.0",
        "  ecutrho = 480.0",
        "  occupations = 'smearing'",
        "  smearing = 'gaussian'",
        "  degauss = 0.01",
        "  nspin = 2",
    ]
    if antiferromagnetic:
        lines += [
            "  starting_magnetization(2) =  0.5",
            "  starting_magnetization(3) = -0.5",
        ]
    else:
        lines += ["  starting_magnetization(2) = 0.5"]
    lines += [
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-8",
        "  mixing_beta = 0.3",
        "  electron_maxstep = 200",
        "/",
        "ATOMIC_SPECIES",
        "  Ba 137.327  Ba.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "  Fe1 55.845  Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
    ]
    if antiferromagnetic:
        lines.append("  Fe2 55.845  Fe.pbe-spn-kjpaw_psl.0.2.1.UPF")
    lines += [
        "  Si 28.0855  Si.pbe-n-kjpaw_psl.1.0.0.UPF",
        "  O  15.999   O.pbe-n-kjpaw_psl.1.0.0.UPF",
        "",
        "ATOMIC_POSITIONS crystal",
    ]
    for index, site in enumerate(sites["Fe"]):
        name = "Fe2" if (antiferromagnetic and index == 1) else "Fe1"
        lines.append("  %-3s %.6f %.6f %.6f" % (name, *site))
    for label in ("Ba", "Si"):
        for site in sites[label]:
            lines.append("  %-3s %.6f %.6f %.6f" % (label, *site))
    for label in ("O1", "O2", "O3"):
        for site in sites[label]:
            lines.append("  %-3s %.6f %.6f %.6f" % ("O", *site))
    lines += ["", "K_POINTS automatic", "  4 4 6 0 0 0", ""]
    # DFT+U card syntax, required from QE 7.1 onward (the old lda_plus_u /
    # Hubbard_U(i) namelist form is rejected).
    lines += ["HUBBARD (ortho-atomic)", f"U Fe1-3d {hubbard_u}"]
    if antiferromagnetic:
        lines.append(f"U Fe2-3d {hubbard_u}")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def distort(sites: dict[str, list[np.ndarray]], delta_theta_deg: float, bond_length: float):
    """Reshape the FeO4 tetrahedron to a given distortion, leaving the cell fixed.

    Only the four O3 sites coordinating Fe move, and each keeps its azimuth
    while its polar angle and radius are set to the requested values. This is
    the same deformation the cluster model applies, so D(distortion) and
    J(distortion) refer to identical structures and can be compared directly.
    Si-O bonds are left to distort, which is the same approximation the cluster
    model makes by treating FeO4 in isolation.
    """
    cell = np.diag([LATTICE_A, LATTICE_A, LATTICE_C])
    inverse = np.linalg.inv(cell)
    theta = np.radians(54.7356 + delta_theta_deg)
    moved = {label: [site.copy() for site in positions] for label, positions in sites.items()}

    for iron in sites["Fe"]:
        for index, site in enumerate(sites["O3"]):
            for shift in itertools.product([-1, 0, 1], repeat=3):
                vector = (site + np.array(shift) - iron) @ cell
                if np.linalg.norm(vector) > 2.4:
                    continue
                azimuth = np.arctan2(vector[1], vector[0])
                sign = 1.0 if vector[2] >= 0 else -1.0
                rebuilt = np.array([
                    bond_length * np.sin(theta) * np.cos(azimuth),
                    bond_length * np.sin(theta) * np.sin(azimuth),
                    sign * bond_length * np.cos(theta),
                ])
                moved["O3"][index] = np.mod(iron + rebuilt @ inverse - np.array(shift), 1.0)
    return moved


def write_relax_input(
    path: Path,
    sites: dict[str, list[np.ndarray]],
    *,
    strain_c: float,
    hubbard_u: float = 4.0,
) -> None:
    """Relax the internal coordinates at a fixed, strained cell.

    Straining c and letting the ions settle is what actually answers the
    paper's open question. It replaces the assumed Poisson ratio with a
    computed response, and each relaxed structure then gives both its own
    distortion angle and, through a pair of fixed-spin runs, its own J.
    """
    tag = f"relax_c{strain_c:+.3f}".replace(".", "p").replace("+", "p").replace("-", "m")
    scaled_c = LATTICE_C * (1.0 + strain_c)
    lines = [
        "&CONTROL",
        "  calculation = 'relax'",
        f"  prefix = '{tag}'",
        "  outdir = './tmp'",
        "  pseudo_dir = './pseudo'",
        "  forc_conv_thr = 1.0d-4",
        "  nstep = 100",
        "/",
        "&SYSTEM",
        "  ibrav = 6",
        f"  celldm(1) = {LATTICE_A / 0.529177210903:.8f}",
        f"  celldm(3) = {scaled_c / LATTICE_A:.8f}",
        "  nat = 24",
        "  ntyp = 4",
        "  ecutwfc = 60.0",
        "  ecutrho = 480.0",
        "  occupations = 'smearing'",
        "  smearing = 'gaussian'",
        "  degauss = 0.01",
        "  nspin = 2",
        "  starting_magnetization(2) = 0.5",
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-7",
        "  mixing_beta = 0.3",
        "  electron_maxstep = 200",
        "/",
        "&IONS",
        "  ion_dynamics = 'bfgs'",
        "/",
        "ATOMIC_SPECIES",
        "  Ba 137.327  Ba.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "  Fe1 55.845  Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
        "  Si 28.0855  Si.pbe-n-kjpaw_psl.1.0.0.UPF",
        "  O  15.999   O.pbe-n-kjpaw_psl.1.0.0.UPF",
        "",
        "ATOMIC_POSITIONS crystal",
    ]
    for site in sites["Fe"]:
        lines.append("  %-3s %.6f %.6f %.6f" % ("Fe1", *site))
    for label in ("Ba", "Si"):
        for site in sites[label]:
            lines.append("  %-3s %.6f %.6f %.6f" % (label, *site))
    for label in ("O1", "O2", "O3"):
        for site in sites[label]:
            lines.append("  %-3s %.6f %.6f %.6f" % ("O", *site))
    lines += ["", "K_POINTS automatic", "  3 3 4 0 0 0", ""]
    lines += ["HUBBARD (ortho-atomic)", f"U Fe1-3d {hubbard_u}", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    sites = build()
    counts = {label: len(positions) for label, positions in sites.items()}
    total = sum(counts.values())
    print("site multiplicities:", counts, "-> total", total)
    if total != 24:
        raise SystemExit(f"expected 24 atoms in Ba2FeSi2O7 (Z=2), built {total}")

    neighbours = coordination(sites)
    print("\nFeO4 coordination:")
    for distance, label, vector in neighbours[:4]:
        angle = float(np.degrees(np.arccos(abs(vector[2]) / distance)))
        print(f"   {label}  d = {distance:.4f} A   theta = {angle:6.2f} deg   delta = {angle - 54.7356:+.3f}")
    distances = [entry[0] for entry in neighbours[:4]]
    angles = [float(np.degrees(np.arccos(abs(entry[2][2]) / entry[0]))) for entry in neighbours[:4]]
    print(f"\n   mean Fe-O        = {np.mean(distances):.4f} A    (cluster model uses 1.99626)")
    print(f"   mean delta theta = {np.mean(angles) - 54.7356:.3f} deg  (cluster model uses 7.8353)")

    cell = np.diag([LATTICE_A, LATTICE_A, LATTICE_C])
    iron = sites["Fe"]
    separation = min(
        float(np.linalg.norm((iron[1] + np.array(shift) - iron[0]) @ cell))
        for shift in itertools.product([-1, 0, 1], repeat=3)
    )
    print(f"   Fe-Fe nearest    = {separation:.4f} A")

    # J is sensitive to U, so the scan is generated up front.
    (HERE / "pseudo").mkdir(exist_ok=True)
    for hubbard_u in (3.0, 4.0, 5.0, 6.0):
        suffix = f"u{hubbard_u:g}".replace(".", "p")
        write_pw_input(HERE / f"scf_fm_{suffix}.in", sites, antiferromagnetic=False, hubbard_u=hubbard_u)
        write_pw_input(HERE / f"scf_afm_{suffix}.in", sites, antiferromagnetic=True, hubbard_u=hubbard_u)

    # J(distortion): reshape FeO4 exactly as the cluster model does, then run
    # fixed-geometry scf pairs. U = 5 eV, which reproduces the measured J to
    # 2.6% once the S=1/S=2 convention is accounted for.
    # Same theta-R relation the cluster model uses: a straight line through the
    # two measured structures, pinned at BFSO.
    bfso_angle, bfso_bond = 7.8353, 1.99626
    sfso_angle, sfso_bond = 4.8600, 1.96900
    slope = (bfso_bond - sfso_bond) / (bfso_angle - sfso_angle)

    for angle in (4.0, 12.0):
        bond = bfso_bond + slope * (angle - bfso_angle)
        shifted = distort(sites, angle, bond)
        neighbours = coordination(shifted)
        realised = float(np.degrees(np.arccos(abs(neighbours[0][2][2]) / neighbours[0][0])))
        print(f"   distortion {angle:5.2f} deg -> realised {realised - 54.7356:.3f} deg, "
              f"Fe-O {neighbours[0][0]:.4f} A")
        suffix = f"th{angle:g}".replace(".", "p")
        write_pw_input(HERE / f"scf_fm_{suffix}.in", shifted, antiferromagnetic=False, hubbard_u=5.0)
        write_pw_input(HERE / f"scf_afm_{suffix}.in", shifted, antiferromagnetic=True, hubbard_u=5.0)
    payload = {
        "a": LATTICE_A,
        "c": LATTICE_C,
        "space_group": "P-42_1m (113)",
        "source": "Jang et al., Phys. Rev. B 104, 214434 (2021), Table II, T = 1.7 K",
        "sites": {label: [list(map(float, site)) for site in positions] for label, positions in sites.items()},
    }
    (HERE / "bfso_sites.json").write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print("\nwrote scf_fm.in, scf_afm.in, bfso_sites.json")


if __name__ == "__main__":
    main()

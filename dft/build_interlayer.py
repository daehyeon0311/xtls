"""Doubled cell along c, for the interlayer exchange J'.

    python dft/build_interlayer.py

The critical ratio is not a constant: quantum Monte Carlo gives alpha_c = 0.18
in two dimensions and 0.10 in three, and Do et al.'s 0.158 for Ba2FeSi2O7 sits
between them at J'/J = 0.1. Straining c changes the interlayer separation, so
alpha_c moves along with everything else. Computing J' pins it down.

The cell is doubled along c so that two FeSi2O7 layers are present. In-plane
order is held antiferromagnetic throughout -- that is the observed structure --
and only the relative sign between layers is switched:

    E(interlayer FM) - E(interlayer AFM) = 8 J' S(S+1)

with four interlayer bonds per doubled cell, matching the form used for J.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
if str(HERE.parent) not in sys.path:
    sys.path.insert(0, str(HERE.parent))

from dft.build_structure import (  # noqa: E402
    LATTICE_A,
    LATTICE_C,
    REPRESENTATIVES,
    build,
)


def doubled_sites(sites: dict[str, list[np.ndarray]]) -> dict[str, list[np.ndarray]]:
    """Stack two copies along c, with fractional z halved."""
    doubled: dict[str, list[np.ndarray]] = {}
    for label, positions in sites.items():
        stacked = []
        for site in positions:
            stacked.append(np.array([site[0], site[1], site[2] / 2.0]))
            stacked.append(np.array([site[0], site[1], site[2] / 2.0 + 0.5]))
        doubled[label] = stacked
    return doubled


def write_input(
    path: Path,
    sites: dict[str, list[np.ndarray]],
    *,
    interlayer_antiparallel: bool,
    hubbard_u: float = 5.0,
    c_scale: float = 1.0,
    label: str = "ref",
) -> None:
    """One scf input on the doubled cell.

    The four Fe sites are two per layer. In-plane neighbours always oppose each
    other; `interlayer_antiparallel` flips the second layer as a whole.
    """
    tag = ("inter_afm_" if interlayer_antiparallel else "inter_fm_") + label
    iron = sites["Fe"]
    # doubled_sites interleaves the two layers: even index is z/2, odd is z/2+1/2
    signs = []
    for index in range(len(iron)):
        layer = index % 2
        in_plane = index // 2
        sign = 1 if in_plane == 0 else -1
        if interlayer_antiparallel and layer == 1:
            sign = -sign
        signs.append(sign)

    species = sorted({("Fe1" if s > 0 else "Fe2") for s in signs})
    lines = [
        "&CONTROL",
        "  calculation = 'scf'",
        f"  prefix = '{tag}'",
        "  outdir = './tmp'",
        "  pseudo_dir = './pseudo'",
        "/",
        "&SYSTEM",
        "  ibrav = 6",
        f"  celldm(1) = {LATTICE_A / 0.529177210903:.8f}",
        f"  celldm(3) = {2.0 * LATTICE_C * c_scale / LATTICE_A:.8f}",
        "  nat = 48",
        f"  ntyp = {2 + len(species)}",
        "  ecutwfc = 60.0",
        "  ecutrho = 480.0",
        "  occupations = 'smearing'",
        "  smearing = 'gaussian'",
        "  degauss = 0.01",
        "  nspin = 2",
    ]
    for index, name in enumerate(species, start=2):
        magnetisation = 0.5 if name == "Fe1" else -0.5
        lines.append(f"  starting_magnetization({index}) = {magnetisation:+.1f}")
    lines += [
        "/",
        "&ELECTRONS",
        # J' is an order of magnitude below J, but 1e-7 Ry is still 1.4e-6 eV,
        # far below the splitting being resolved.
        "  conv_thr = 1.0d-7",
        "  mixing_beta = 0.3",
        "  electron_maxstep = 200",
        "/",
        "ATOMIC_SPECIES",
        "  Ba 137.327  Ba.pbe-spn-kjpaw_psl.1.0.0.UPF",
    ]
    for name in species:
        lines.append(f"  {name} 55.845  Fe.pbe-spn-kjpaw_psl.0.2.1.UPF")
    lines += [
        "  Si 28.0855  Si.pbe-n-kjpaw_psl.1.0.0.UPF",
        "  O  15.999   O.pbe-n-kjpaw_psl.1.0.0.UPF",
        "",
        "ATOMIC_POSITIONS crystal",
    ]
    for site, sign in zip(iron, signs):
        lines.append("  %-3s %.6f %.6f %.6f" % ("Fe1" if sign > 0 else "Fe2", *site))
    for label_name in ("Ba", "Si"):
        for site in sites[label_name]:
            lines.append("  %-3s %.6f %.6f %.6f" % (label_name, *site))
    for label_name in ("O1", "O2", "O3"):
        for site in sites[label_name]:
            lines.append("  %-3s %.6f %.6f %.6f" % ("O", *site))
    # Halved reciprocal spacing along c, since the cell doubled there.
    lines += ["", "K_POINTS automatic", "  4 4 3 0 0 0", ""]
    lines += ["HUBBARD (ortho-atomic)"]
    for name in species:
        lines.append(f"U {name}-3d {hubbard_u}")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    sites = doubled_sites(build())
    counts = {label: len(positions) for label, positions in sites.items()}
    total = sum(counts.values())
    print("doubled cell:", counts, "-> total", total)
    if total != 48:
        raise SystemExit(f"expected 48 atoms, built {total}")

    cell = np.diag([LATTICE_A, LATTICE_A, 2.0 * LATTICE_C])
    iron = sites["Fe"]
    print("\nFe sites and their layer:")
    for index, site in enumerate(iron):
        print(f"  {index}  ({site[0]:.4f}, {site[1]:.4f}, {site[2]:.4f})"
              f"   layer {'A' if site[2] < 0.25 or site[2] > 0.75 else 'B'}")
    separation = min(
        float(np.linalg.norm((iron[j] + np.array([0, 0, k]) - iron[i]) @ cell))
        for i in range(len(iron)) for j in range(len(iron)) for k in (-1, 0, 1)
        if not (i == j and k == 0)
    )
    print(f"\nnearest Fe-Fe overall: {separation:.4f} A")

    o3 = np.array(REPRESENTATIVES["O3"])
    in_plane = np.linalg.norm(o3[:2] * LATTICE_A)
    z_reference = o3[2] * LATTICE_C
    for target, name in ((8.225, "ref"), (6.0, "th6")):
        scale = (in_plane / np.tan(np.radians(54.7356 + target))) / z_reference
        for antiparallel in (False, True):
            stem = ("inter_afm_" if antiparallel else "inter_fm_") + name
            write_input(HERE / f"{stem}.in", sites, interlayer_antiparallel=antiparallel,
                        c_scale=scale, label=name)
        print(f"wrote inter_fm_{name}.in / inter_afm_{name}.in  (c x {scale:.4f})")


if __name__ == "__main__":
    main()

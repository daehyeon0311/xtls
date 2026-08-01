"""Single-ion anisotropy from the cluster model.

    python run_anisotropy.py [input file]

The S = 2 manifold of a compressed FeO4 tetrahedron splits into Sz = 0, +-1
and +-2. Fitting those five levels gives the single-ion anisotropy D, the
quartic correction A, and the tunnelling splitting of the Sz = +-2 doublet:

    E(Sz) = D Sz^2 + A Sz^4,     E(+2) - E(-2) = tunnelling splitting

Only the initial-state Hamiltonian is needed, so this runs in a fraction of a
second per geometry and a distortion scan is cheap. Nothing here computes a
spectrum.

Scanning the distortion angle also has to move the bond length: compressing a
tetrahedron along c stretches the in-plane directions. The two measured
structures fix that relation, and `bond_length_for` interpolates through them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
for path in (ROOT, SRC):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import run_xas as runner  # noqa: E402
from xtls_py.engine import (  # noqa: E402
    lowest_eigenpairs,
    one_body_sparse,
    p_ct_one_body_matrix,
    p_ct_xas_initial_basis,
    spin_expand_orbital_matrix,
    spin_matrices,
)
from xtls_py.geometry import crystal_field, reference_scale  # noqa: E402


IDEAL_TETRAHEDRAL_ANGLE_DEG = 54.7356

# Structure whose crystal field calibrates the whole scan. The point-charge
# model already fixes how the field varies with the ligand positions, so it is
# calibrated once here and then left alone; renormalizing at every angle would
# discard exactly the variation the scan is trying to measure.
CALIBRATION_STRUCTURE = "BFSO"

# The two crystallographic structures, which fix how the bond length tracks
# the distortion angle. Values from Phys. Rev. B 104, 214434 (2021) and the
# SFSO structure report, as used in the XLD paper.
REFERENCE_STRUCTURES = {
    "BFSO": {"delta_theta_deg": 7.8353, "bond_length_A": 1.99626},
    "SFSO": {"delta_theta_deg": 4.8600, "bond_length_A": 1.96900},
}


def bond_length_slope(references=REFERENCE_STRUCTURES) -> float:
    """dR/d(delta theta) in A/deg, from the two measured structures."""
    (first, second) = references.values()
    return (first["bond_length_A"] - second["bond_length_A"]) / (
        first["delta_theta_deg"] - second["delta_theta_deg"]
    )


def bond_length_for(
    delta_theta_deg: float,
    slope_factor: float = 1.0,
    anchor: str = CALIBRATION_STRUCTURE,
    references=REFERENCE_STRUCTURES,
) -> float:
    """Fe-O distance implied by a distortion angle.

    A straight line through the two measured structures, pinned at the anchor
    structure. This is the constraint the XLD paper expresses as a fixed
    Poisson ratio: compressing along c stretches the in-plane bonds, so theta
    and R cannot be varied independently.

    `slope_factor` rescales dR/dtheta to test how much that assumption matters.
    1.0 follows the measured structures; 0.0 holds the bond length fixed and
    varies the angle alone.
    """
    reference = references[anchor]
    slope = bond_length_slope(references) * float(slope_factor)
    return reference["bond_length_A"] + slope * (delta_theta_deg - reference["delta_theta_deg"])


def spin_z_operator(n_orbitals: int = 26) -> np.ndarray:
    """One-body Sz over every spin-orbital of the 2p + 3d + ligand space."""
    _sx, _sy, sz = spin_matrices()
    block = np.zeros((n_orbitals, n_orbitals), dtype=complex)
    for orbital in range(n_orbitals // 2):
        block[2 * orbital : 2 * orbital + 2, 2 * orbital : 2 * orbital + 2] = sz
    return block


def calibration_scale(ten_dq: float | None = None, structure: str = CALIBRATION_STRUCTURE) -> float:
    """Crystal-field scale set once at the reference structure."""
    reference = REFERENCE_STRUCTURES[structure]
    runner.__dict__["ligand_angle_offset_deg"] = reference["delta_theta_deg"]
    runner.__dict__["ligand_radius"] = reference["bond_length_A"]
    positions = runner._ligand_positions()
    value = runner.ten_dq if ten_dq is None else float(ten_dq)
    return reference_scale(positions, value, runner.r2, runner.r4)


def calibrate_ten_dq(target_D_meV: float, structure: str = CALIBRATION_STRUCTURE,
                     bracket: tuple[float, float] = (0.05, 0.55)) -> float:
    """Find the `ten_dq` that reproduces a target D at the reference structure.

    `ten_dq` is a fitted quantity, not an observable, so its value depends on
    the normalization convention it was fitted under. Changing the convention
    means refitting it against something measurable -- here, D.

    D is not monotonic in `ten_dq`: it falls to a minimum near 0.6 eV and rises
    again, and past roughly 1.5 eV the field is strong enough that the S = 2
    manifold stops being the five lowest states at all. The default bracket
    stays on the weak-field branch, which is the physically relevant one for a
    high-spin d6 tetrahedron.
    """
    from scipy.optimize import brentq

    reference = REFERENCE_STRUCTURES[structure]

    def residual(value: float) -> float:
        levels, _sz = multiplet_levels(
            reference["delta_theta_deg"], reference["bond_length_A"], ten_dq=value
        )
        return fit_anisotropy(levels)["D_meV"] - target_D_meV

    return float(brentq(residual, bracket[0], bracket[1], xtol=1e-7))


def multiplet_levels(
    delta_theta_deg: float,
    bond_length_A: float | None = None,
    n_levels: int = 8,
    ten_dq: float | None = None,
    scale: float | None = None,
):
    """Lowest eigenvalues of the initial-state Hamiltonian, in meV.

    Returns the level energies relative to the ground state together with
    their Sz expectation values. Passing `scale` holds the crystal-field
    calibration fixed while the geometry changes.
    """
    if bond_length_A is None:
        bond_length_A = bond_length_for(delta_theta_deg)

    runner.__dict__["ligand_angle_offset_deg"] = float(delta_theta_deg)
    runner.__dict__["ligand_radius"] = float(bond_length_A)

    holes = runner._effective_max_ligand_holes()
    positions = runner._ligand_positions()
    field_ten_dq = runner.ten_dq if ten_dq is None else float(ten_dq)
    h_crystal_orbital, *_ = crystal_field(
        positions, field_ten_dq, runner.r2, runner.r4, scale=scale
    )
    h_crystal = spin_expand_orbital_matrix(h_crystal_orbital)
    h_hybridization, _label = runner._build_hybridization_matrix(positions)

    basis = p_ct_xas_initial_basis(runner.n_d_electrons, max_ligand_holes=holes)
    slater = runner._configuration_slater_entries(
        p_electrons=6, first_d=runner.n_d_electrons, max_holes=holes
    )
    hamiltonian = runner._build_initial_hamiltonian(basis, h_crystal, h_hybridization, slater)

    energies, vectors = lowest_eigenpairs(hamiltonian, k=min(n_levels, len(basis) - 1))
    sz_matrix = one_body_sparse(basis, p_ct_one_body_matrix(
        h_p=spin_z_operator(6), h_d=spin_z_operator(10), h_ligand=spin_z_operator(10)
    ))
    sz_values = np.array(
        [float(np.real(np.vdot(vectors[:, i], sz_matrix @ vectors[:, i]))) for i in range(vectors.shape[1])]
    )
    return (energies - energies[0]) * 1000.0, sz_values


def fit_anisotropy(levels_meV: np.ndarray) -> dict[str, float]:
    """Extract D, the quartic term, and the tunnelling splitting.

    The five lowest levels are the S = 2 manifold, ordered Sz = 0, +-1, +-2.
    Solving `E(Sz) = D Sz^2 + A Sz^4` on the doublet centroids gives

        D = (16 E1 - E2) / 12,   A = (E2 - 4 E1) / 12
    """
    if len(levels_meV) < 5:
        raise ValueError("need at least five levels for the S = 2 manifold")
    manifold = np.asarray(levels_meV[:5], dtype=float)
    e1 = float(np.mean(manifold[1:3]))
    e2 = float(np.mean(manifold[3:5]))
    return {
        "D_meV": (16.0 * e1 - e2) / 12.0,
        "A_meV": (e2 - 4.0 * e1) / 12.0,
        "D_simple_meV": e1,
        "tunnelling_splitting_meV": float(manifold[4] - manifold[3]),
        "doublet1_splitting_meV": float(manifold[2] - manifold[1]),
        "E1_meV": e1,
        "E2_meV": e2,
        "manifold_width_meV": e2,
    }


def scan(
    angles_deg,
    scale: float | None = None,
    slope_factor: float = 1.0,
    verbose: bool = True,
) -> list[dict[str, float]]:
    rows = []
    for angle in angles_deg:
        radius = bond_length_for(angle, slope_factor=slope_factor)
        levels, sz = multiplet_levels(angle, bond_length_A=radius, scale=scale)
        row = {"delta_theta_deg": float(angle), "bond_length_A": radius}
        row.update(fit_anisotropy(levels))
        row["gap_to_next_meV"] = float(levels[5]) if len(levels) > 5 else float("nan")
        rows.append(row)
        if verbose:
            print(
                f"  dtheta={angle:6.3f} deg  R={radius:.5f} A   "
                f"D={row['D_meV']:.4f}  A={row['A_meV']:+.5f}  "
                f"tunnel={row['tunnelling_splitting_meV']:.5f} meV"
            )
    return rows


def main(argv: list[str] | None = None) -> None:
    arguments = list(sys.argv[1:] if argv is None else argv)
    input_file = Path(arguments[0]) if arguments else ROOT / "inputs" / "Fe_Ba2FeSi2O7.py"
    runner._load_input_file(input_file)
    runner.__dict__["n_analyzed_states"] = 8

    target = 1.45  # D_BFSO from the XLD paper, consistent with neutron 1.42 meV
    print()
    print(f"calibrating ten_dq against D({CALIBRATION_STRUCTURE}) = {target} meV")
    print(f"  input file value : ten_dq = {runner.ten_dq:.4f} eV")
    calibrated = calibrate_ten_dq(target)
    print(f"  calibrated       : ten_dq = {calibrated:.4f} eV")
    scale = calibration_scale(calibrated)
    print(f"  crystal-field scale held fixed at {scale:.6g}")

    print()
    print("reference structures (calibrated, scale held fixed)")
    for name, structure in REFERENCE_STRUCTURES.items():
        levels, sz = multiplet_levels(
            structure["delta_theta_deg"], structure["bond_length_A"], scale=scale
        )
        fit = fit_anisotropy(levels)
        print(f"  {name}:  dtheta={structure['delta_theta_deg']:.4f} deg  R={structure['bond_length_A']:.5f} A")
        print("    levels (meV):", "  ".join(f"{value:.4f}" for value in levels[:5]))
        print("    <Sz>        :", "  ".join(f"{value:+.3f}" for value in sz[:5]))
        print(
            f"    D = {fit['D_meV']:.4f} meV   quartic A = {fit['A_meV']:+.5f} meV   "
            f"tunnelling = {fit['tunnelling_splitting_meV']:.5f} meV"
        )
    print()
    print("  paper values: D_BFSO = 1.45 meV, D_SFSO = 1.22 meV")
    print("  SFSO is a prediction here -- only BFSO was used for the calibration.")

    print()
    print("distortion scan")
    angles = np.arange(4.0, 12.01, 0.25)
    rows = scan(angles, scale=scale, verbose=False)
    critical = crossing(rows, target=1.65)
    for row in rows[:: max(1, len(rows) // 12)]:
        print(
            f"  dtheta={row['delta_theta_deg']:6.2f}  R={row['bond_length_A']:.5f}  "
            f"D={row['D_meV']:.4f} meV"
        )
    if critical is not None:
        print(f"\n  D reaches the critical 1.65 meV at dtheta = {critical:.3f} deg")
    else:
        print("\n  D does not reach the critical 1.65 meV within 4-12 deg")

    output = ROOT / "outputs" / "anisotropy"
    output.mkdir(parents=True, exist_ok=True)
    save_table(output / "anisotropy_scan.txt", rows)
    plot(output / "anisotropy_scan.png", rows, scale=scale)

    print()
    print("Poisson-ratio sensitivity")
    print(f"  measured slope dR/dtheta = {bond_length_slope():.6f} A/deg")
    print("  factor   R(12 deg)   D(4 deg)  D(12 deg)   change   D_c at")
    sensitivity = {}
    for factor in (0.0, 0.5, 1.0, 1.5, 2.0):
        sensitivity_rows = scan(angles, scale=scale, slope_factor=factor, verbose=False)
        sensitivity[factor] = sensitivity_rows
        first = sensitivity_rows[0]["D_meV"]
        last = sensitivity_rows[-1]["D_meV"]
        reached = crossing(sensitivity_rows, target=1.65)
        print(
            f"  {factor:5.1f}   {bond_length_for(12.0, factor):8.5f}   "
            f"{first:7.4f}   {last:7.4f}   {100 * (last / first - 1):+6.1f}%   "
            + (f"{reached:6.2f} deg" if reached is not None else "  never")
        )
    save_table(output / "poisson_sensitivity.txt", [
        {"slope_factor": factor, **row} for factor, rows_ in sensitivity.items() for row in rows_
    ])
    plot_sensitivity(output / "poisson_sensitivity.png", sensitivity)
    print(f"\nsaved: {output / 'anisotropy_scan.txt'}")
    print(f"saved: {output / 'poisson_sensitivity.txt'}")


def crossing(rows, target: float) -> float | None:
    """Distortion angle where D crosses a target value, linearly interpolated."""
    angles = np.array([row["delta_theta_deg"] for row in rows])
    values = np.array([row["D_meV"] for row in rows])
    above = np.where(values >= target)[0]
    if not above.size or above[0] == 0:
        return None
    index = above[0]
    span = values[index] - values[index - 1]
    if abs(span) < 1e-12:
        return float(angles[index])
    weight = (target - values[index - 1]) / span
    return float(angles[index - 1] + weight * (angles[index] - angles[index - 1]))


def plot_sensitivity(path: Path, sensitivity: dict[float, list[dict[str, float]]]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    figure, axis = plt.subplots(figsize=(5.6, 4.0))
    for factor, rows in sorted(sensitivity.items()):
        angles = [row["delta_theta_deg"] for row in rows]
        values = [row["D_meV"] for row in rows]
        label = "bond length fixed" if factor == 0.0 else f"$dR/d\\theta$ x {factor:g}"
        axis.plot(angles, values, lw=1.6 if factor == 1.0 else 1.1,
                  ls="-" if factor == 1.0 else "--", label=label)
    axis.axhline(1.65, color="0.4", ls=":", lw=1.0)
    axis.text(4.1, 1.66, "critical $D_c$", fontsize=8, color="0.3")
    axis.set_xlabel(r"distortion angle $\Delta\theta$ (deg)")
    axis.set_ylabel("$D$ (meV)")
    axis.legend(frameon=False, fontsize=8)
    axis.tick_params(direction="in", top=True, right=True)
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)


def save_table(path: Path, rows) -> None:
    keys = list(rows[0])
    lines = [" ".join(keys)]
    for row in rows:
        lines.append(" ".join(f"{row[key]:.8g}" for key in keys))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot(path: Path, rows, scale: float | None = None) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    angles = [row["delta_theta_deg"] for row in rows]
    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))

    axes[0].plot(angles, [row["D_meV"] for row in rows], "-", color="crimson", lw=1.8)
    axes[0].axhline(1.65, color="0.4", ls="--", lw=1.0)
    axes[0].text(angles[0], 1.66, "critical $D_c$", fontsize=8, color="0.3")
    for name, structure in REFERENCE_STRUCTURES.items():
        levels, _sz = multiplet_levels(
            structure["delta_theta_deg"], structure["bond_length_A"], scale=scale
        )
        value = fit_anisotropy(levels)["D_meV"]
        axes[0].plot(structure["delta_theta_deg"], value, "o", color="black", ms=5)
        axes[0].annotate(name, (structure["delta_theta_deg"], value),
                         textcoords="offset points", xytext=(6, -10), fontsize=8)
    axes[0].set_xlabel(r"distortion angle $\Delta\theta$ (deg)")
    axes[0].set_ylabel("$D$ (meV)")

    axes[1].plot(angles, [row["A_meV"] for row in rows], "-", color="royalblue", lw=1.5, label="quartic $A$")
    axes[1].plot(angles, [row["tunnelling_splitting_meV"] for row in rows], "-",
                 color="seagreen", lw=1.5, label=r"$S_z=\pm2$ splitting")
    axes[1].axhline(0.0, color="0.7", lw=0.6)
    axes[1].set_xlabel(r"distortion angle $\Delta\theta$ (deg)")
    axes[1].set_ylabel("energy (meV)")
    axes[1].legend(frameon=False, fontsize=9)

    for axis in axes:
        axis.tick_params(direction="in", top=True, right=True)
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)


if __name__ == "__main__":
    main()

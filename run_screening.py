"""Zero-field splitting across the A2MX2O7 melilite family.

    python run_screening.py

Holds the FeO4 geometry of Ba2FeSi2O7 and swaps the transition-metal ion,
scanning the compressive distortion for each. The question is which d count
puts the single-ion anisotropy closest to the D/J critical ratio that
separates the quantum-paramagnetic and antiferromagnetic ground states.

Only the initial-state Hamiltonian is needed, so the whole family scan runs in
a couple of minutes. The charge-transfer parameters are deliberately held
common across ions: the point is the systematic trend with d count, not a fit
to any one compound.
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

import run_anisotropy as anisotropy  # noqa: E402
import run_xas as runner  # noqa: E402
from xtls_py.engine import (  # noqa: E402
    get_appendix_a_3d,
    lowest_eigenpairs,
    one_body_sparse,
    p_ct_one_body_matrix,
    p_ct_xas_initial_basis,
    spin_expand_orbital_matrix,
)
from xtls_py.geometry import crystal_field  # noqa: E402


# High-spin d counts in a tetrahedral oxygen environment, where the ligand
# field is too weak to pair electrons. `term` is the free-ion ground term in Td.
IONS = [
    {"ion": "Mn2+", "element": "Mn", "d": 5, "spin": 2.5, "term": "6A1", "note": "orbital singlet, half filled"},
    {"ion": "Fe2+", "element": "Fe", "d": 6, "spin": 2.0, "term": "5E", "note": "BFSO/SFSO, orbital doublet"},
    {"ion": "Co2+", "element": "Co", "d": 7, "spin": 1.5, "term": "4A2", "note": "orbital singlet"},
    {"ion": "Ni2+", "element": "Ni", "d": 8, "spin": 1.0, "term": "3T1", "note": "orbital triplet"},
]

MAX_LIGAND_HOLES = 1


def sz_squared_ladder(spin: float) -> np.ndarray:
    """Sz^2 for the 2S+1 levels, ordered from smallest |Sz| upward.

    Easy-plane anisotropy (D > 0) puts small |Sz| lowest, and every level
    except Sz = 0 comes as a doublet.
    """
    values = []
    sz = 0.0 if float(spin).is_integer() else 0.5
    while sz <= spin + 1e-9:
        values.extend([sz**2] if sz == 0.0 else [sz**2, sz**2])
        sz += 1.0
    return np.array(values[: int(round(2 * spin + 1))])


def fit_zero_field_splitting(levels_meV: np.ndarray, spin: float) -> dict[str, float]:
    """Least-squares D for E(Sz) = D Sz^2 over the 2S+1 lowest levels."""
    count = int(round(2 * spin + 1))
    if len(levels_meV) < count + 1:
        raise ValueError("need the whole manifold plus one level above it")
    manifold = np.asarray(levels_meV[:count], dtype=float)
    ladder = sz_squared_ladder(spin)
    weights = ladder - ladder[0]
    denominator = float(np.dot(weights, weights))
    slope = float(np.dot(weights, manifold - manifold[0]) / denominator) if denominator else float("nan")
    predicted = manifold[0] + slope * weights
    residual = float(np.max(np.abs(manifold - predicted)))
    return {
        "D_meV": slope,
        "manifold_width_meV": float(manifold[-1] - manifold[0]),
        "fit_residual_meV": residual,
        "gap_above_meV": float(levels_meV[count] - manifold[-1]),
    }


def levels_for(entry: dict, delta_theta_deg: float, scale: float, n_levels: int) -> np.ndarray:
    """Lowest levels of one ion at one distortion, in meV."""
    runner.__dict__["element"] = entry["element"]
    runner.__dict__["n_d_electrons"] = entry["d"]
    runner.__dict__["ligand_angle_offset_deg"] = float(delta_theta_deg)
    runner.__dict__["ligand_radius"] = anisotropy.bond_length_for(delta_theta_deg)

    positions = runner._ligand_positions()
    h_crystal = spin_expand_orbital_matrix(
        crystal_field(positions, runner.ten_dq, runner.r2, runner.r4, scale=scale)[0]
    )
    h_hybridization, _label = runner._build_hybridization_matrix(positions)

    holes = min(MAX_LIGAND_HOLES, 10 - entry["d"])
    basis = p_ct_xas_initial_basis(entry["d"], max_ligand_holes=holes)
    slater = {
        hole: runner._scale_multipoles(get_appendix_a_3d(entry["element"], 6, entry["d"] + hole))
        for hole in range(holes + 1)
    }
    hamiltonian = runner._build_initial_hamiltonian(basis, h_crystal, h_hybridization, slater)
    energies, _vectors = lowest_eigenpairs(hamiltonian, k=min(n_levels, len(basis) - 1))
    return (energies - energies[0]) * 1000.0


def main() -> None:
    runner._load_input_file(ROOT / "inputs" / "Fe_Ba2FeSi2O7.py")
    scale = anisotropy.calibration_scale(runner.ten_dq)
    print(f"crystal-field scale fixed at the BFSO structure: {scale:.6g}")
    print(f"common charge transfer: Delta = {runner.delta} eV, U_dd = {runner.u_charge_transfer} eV")
    print(f"ligand holes: {MAX_LIGAND_HOLES}\n")

    reference_angle = anisotropy.REFERENCE_STRUCTURES["BFSO"]["delta_theta_deg"]

    print("at the BFSO distortion angle")
    print("  ion    d  S    term  D (meV)   manifold  gap above  fit resid")
    for entry in IONS:
        count = int(round(2 * entry["spin"] + 1))
        levels = levels_for(entry, reference_angle, scale, n_levels=count + 4)
        fit = fit_zero_field_splitting(levels, entry["spin"])
        print(
            f"  {entry['ion']:<6} {entry['d']}  {entry['spin']:<4.1f} {entry['term']:<5} "
            f"{fit['D_meV']:8.4f}  {fit['manifold_width_meV']:8.3f}  "
            f"{fit['gap_above_meV']:9.2f}  {fit['fit_residual_meV']:8.4f}"
        )
        print(f"         levels: {'  '.join(f'{value:.3f}' for value in levels[:count])}")

    print("\ndistortion scan")
    angles = np.arange(4.0, 12.01, 0.5)
    table: dict[str, list[dict[str, float]]] = {}
    for entry in IONS:
        count = int(round(2 * entry["spin"] + 1))
        rows = []
        for angle in angles:
            levels = levels_for(entry, angle, scale, n_levels=count + 4)
            fit = fit_zero_field_splitting(levels, entry["spin"])
            rows.append({"delta_theta_deg": float(angle), **fit})
        table[entry["ion"]] = rows
        first, last = rows[0]["D_meV"], rows[-1]["D_meV"]
        worst_residual = max(row["fit_residual_meV"] for row in rows)
        smallest_gap = min(row["gap_above_meV"] for row in rows)
        print(
            f"  {entry['ion']:<6} D(4 deg) = {first:8.4f}   D(12 deg) = {last:8.4f}   "
            f"change {100 * (last / first - 1) if first else float('nan'):+7.1f}%   "
            f"max resid {worst_residual:.4f}   min gap {smallest_gap:.1f} meV"
        )

    output = ROOT / "outputs" / "screening"
    output.mkdir(parents=True, exist_ok=True)
    save(output / "family_scan.txt", table)
    plot(output / "family_scan.png", table)
    print(f"\nsaved: {output / 'family_scan.txt'}")


def save(path: Path, table: dict[str, list[dict[str, float]]]) -> None:
    keys = ["delta_theta_deg", "D_meV", "manifold_width_meV", "fit_residual_meV", "gap_above_meV"]
    lines = ["ion " + " ".join(keys)]
    for ion, rows in table.items():
        for row in rows:
            lines.append(ion + " " + " ".join(f"{row[key]:.8g}" for key in keys))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot(path: Path, table: dict[str, list[dict[str, float]]]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for entry in IONS:
        rows = table[entry["ion"]]
        angles = [row["delta_theta_deg"] for row in rows]
        values = [row["D_meV"] for row in rows]
        label = f"{entry['ion']} (d$^{entry['d']}$, S={entry['spin']:g})"
        axes[0].plot(angles, values, lw=1.6, label=label)
        axes[1].plot(angles, np.abs(values), lw=1.6, label=label)
    axes[0].axhline(1.65, color="0.4", ls="--", lw=1.0)
    axes[0].text(4.1, 1.70, "$D_c$ for BFSO", fontsize=8, color="0.3")
    axes[0].set_ylabel("$D$ (meV)")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("$|D|$ (meV, log)")
    for axis in axes:
        axis.set_xlabel(r"distortion angle $\Delta\theta$ (deg)")
        axis.legend(frameon=False, fontsize=8)
        axis.tick_params(direction="in", top=True, right=True)
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)


if __name__ == "__main__":
    main()

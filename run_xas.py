"""Spyder-friendly charge-transfer L-edge XAS/LD runner.

Choose an input file below, open this file in Spyder, and press F5.
"""

from __future__ import annotations

import json
import runpy
import sys
from dataclasses import replace
from math import comb
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from xtls_py.engine import (  # noqa: E402
    available_appendix_a_3d,
    continued_fraction_spectrum,
    d_angular_momentum_matrices,
    d_ligand_hybridization_from_orbital_matrix,
    d_ligand_hybridization_matrix,
    d_shell_slater_tensor,
    d_spin_orbit_matrix,
    diagonal_sparse,
    get_appendix_a_3d,
    lowest_eigenpairs,
    one_body_between_bases_sparse,
    one_body_sparse,
    p_ct_charge_transfer_energy_values,
    p_ct_dipole_matrix,
    p_ct_ligand_hole_count,
    p_ct_one_body_matrix,
    p_ct_spin_orbital_labels,
    p_ct_xas_final_basis,
    p_ct_xas_initial_basis,
    p_spin_orbit_matrix,
    pd_shell_slater_tensor,
    sector_energy,
    spin_matrices,
    spin_expand_orbital_matrix,
    two_body_sparse,
)
from xtls_py.geometry import crystal_field, hybridization  # noqa: E402
from xtls_py.spectrum import (  # noqa: E402
    GAUSSIAN_HWHM_TO_SIGMA,
    interpolate_width,
    save_xy,
    variable_gaussian_broaden,
)


# ---------------------------------------------------------------------------
# Input file. Put cluster/calculation/plot parameters there.

INPUT_FILE = ROOT / "inputs" / "Fe_Ba2FeSi2O7.py"


def main() -> None:
    _load_input_file()
    _validate_parameters()
    effective_estimate_holes = _effective_max_ligand_holes()
    estimated_initial, estimated_final = _estimate_basis_sizes(effective_estimate_holes)
    print("parameter check: OK")
    print("element/configuration:", f"{element} d{n_d_electrons}", f"h={effective_estimate_holes}")
    print("geometry:", coordination_geometry, _ligand_radius_summary())
    if make_experimental_geometry_curves:
        print("LD geometry:", f"grazing={grazing_angle_deg:g} deg", f"inplane={inplane_curve}")
    print("estimated effective max_ligand_holes:", effective_estimate_holes)
    print("estimated basis size:", f"initial={estimated_initial}", f"final={estimated_final}")
    if overlay_xtls and not Path(xtls_path).exists():
        print("overlay warning: XTLS overlay file not found:", xtls_path)
    if overlay_xtls and plot_ld_stacked:
        print("overlay note: XTLS overlay is drawn only in the two-panel comparison plot, not in stacked LD mode.")
    if estimate_only:
        print("estimate_only=True, skipping Lanczos calculation.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    display_energy, curves, metadata = calculate_spectrum()
    effective_max_ligand_holes = int(metadata["effective_max_ligand_holes"])

    stem = (
        f"{element}_d{n_d_electrons}_CT_L_lanczos_h{effective_max_ligand_holes}"
        f"_{_scale_tag()}"
    ).replace(".", "p")
    out_txt = output_dir / f"{stem}_spectrum.txt"
    out_png = output_dir / f"{stem}_spectrum.png"
    out_params = output_dir / f"{stem}_parameters.json"
    out_config = output_dir / f"{stem}_configuration_energies.txt"
    out_analysis = output_dir / f"{stem}_state_analysis.txt"

    if save_txt:
        save_xy(out_txt, display_energy, curves)
    if save_parameters:
        _save_parameters(out_params, metadata)
    if save_configuration_energies:
        _save_configuration_energies(out_config, metadata["configuration_energies"])
    if save_state_analysis:
        _save_state_analysis(out_analysis, metadata["state_analysis"])

    print("element:", element)
    if ion is not None:
        print("ion:", ion)
    print("n_d_electrons:", n_d_electrons)
    print("requested max_ligand_holes:", max_ligand_holes)
    print("effective max_ligand_holes:", metadata["effective_max_ligand_holes"])
    if metadata["effective_max_ligand_holes"] < max_ligand_holes:
        print("ligand hole truncation:", metadata["ligand_hole_truncation_reason"])
    print("ten_dq:", ten_dq)
    print("delta:", delta)
    print("pd_sigma:", pd_sigma)
    print("pd_ratio:", pd_ratio)
    print("u_charge_transfer:", u_charge_transfer)
    print("core_hole_potential:", core_hole_potential)
    print("fdd2_scale:", fdd2_scale)
    print("fdd4_scale:", fdd4_scale)
    print("fpd2_scale:", fpd2_scale)
    print("gpd1_scale:", gpd1_scale)
    print("gpd3_scale:", gpd3_scale)
    print("hybridization_mode:", metadata["hybridization_mode"])
    print("spectrum_method:", metadata["spectrum_method"])
    print("hopping:", hopping)
    print("initial states:", metadata["initial_states"])
    print("analyzed states:", metadata["analyzed_states"])
    print("temperature K:", temperature_kelvin)
    print("initial basis size:", metadata["initial_basis_size"])
    print("final basis size:", metadata["final_basis_size"])
    print("initial ground energy:", round(metadata["initial_ground_energy"], 6))
    print("final ground energy:", round(metadata["final_ground_energy"], 6))
    print("onset energy:", round(metadata["onset"], 6))
    print("n_recursion:", n_recursion)
    print("lorentzian_hwhm:", _format_width_summary(metadata["broadening"]["lorentzian_hwhm"]))
    print("gaussian_sigma:", _format_width_summary(metadata["broadening"]["gaussian_sigma"]))
    if metadata["experimental_geometry"]["enabled"]:
        print(
            "grazing geometry:",
            f"angle={metadata['experimental_geometry']['grazing_angle_deg']} deg,",
            f"pi = {metadata['experimental_geometry']['pi_inplane_weight']:.3f}*inplane"
            f" + {metadata['experimental_geometry']['pi_outofplane_weight']:.3f}*c",
        )
    print("energy range:", energy_min, energy_max, energy_step)
    print("configuration energies:")
    for row in metadata["configuration_energies"]:
        print(f"  {row['label']} {row['configuration']}: {row['energy_eV']:.6g} eV")
    print("initial low-state splitting:")
    for row in metadata["state_analysis"]["initial"][: min(5, len(metadata["state_analysis"]["initial"]))]:
        print(
            f"  state {row['index']}: {row['energy_relative_meV']:.3f} meV, "
            f"config={_compact_config_percent(row)}"
        )
    if save_txt:
        print("saved:", out_txt)
    if save_parameters:
        print("saved:", out_params)
    if save_configuration_energies:
        print("saved:", out_config)
    if save_state_analysis:
        print("saved:", out_analysis)

    _plot(display_energy, curves, out_png)


_INPUT_KEYS = {
    "element",
    "n_d_electrons",
    "max_ligand_holes",
    "ten_dq",
    "delta",
    "pd_sigma",
    "pd_ratio",
    "d_ref",
    "u_charge_transfer",
    "core_hole_potential",
    "coordination_geometry",
    "ligand_radius",
    "ligand_angle_offset_deg",
    "r2",
    "r4",
    "fdd2_scale",
    "fdd4_scale",
    "fpd2_scale",
    "gpd1_scale",
    "gpd3_scale",
    "so3d_scale",
    "so2p_scale",
    "hybridization_mode",
    "hopping",
    "n_initial_states",
    "temperature_kelvin",
    "spectrum_method",
    "n_recursion",
    "lorentzian_hwhm",
    "gaussian_sigma",
    "gaussian_hwhm",
    "lorentzian_hwhm_points",
    "gaussian_hwhm_points",
    "energy_min",
    "energy_max",
    "energy_step",
    "energy_shift",
    "normalize",
    "make_experimental_geometry_curves",
    "grazing_angle_deg",
    "inplane_curve",
    "save_txt",
    "save_png",
    "save_parameters",
    "save_configuration_energies",
    "save_state_analysis",
    "estimate_only",
    "n_analyzed_states",
    "show_plot",
    "show_cluster_inset",
    "plot_ld_stacked",
    "plot_iso_only",
    "plot_xas_offset",
    "plot_ld_offset",
    "plot_ld_scale",
    "plot_use_absolute_energy",
    "plot_absolute_energy_offset",
    "plot_relative_energy_min",
    "plot_relative_energy_max",
    "output_dir",
    "overlay_xtls",
    "xtls_path",
    "xtls_energy_shift",
    "xtls_scale",
    "xtls_iso_column",
}

_OPTIONAL_INPUT_DEFAULTS = {
    "case_name": None,
    "ion": None,
    "ligand_positions_xyz": None,
    "ligand_positions_spherical": None,
}

_ION_EDGE_L3_EV = {
    "K": 295.0,
    "Ca": 346.0,
    "Sc": 399.0,
    "Ti": 458.0,
    "V": 515.0,
    "Cr": 576.0,
    "Mn": 640.0,
    "Fe": 707.0,
    "Co": 779.0,
    "Ni": 852.0,
    "Cu": 931.0,
    "Zn": 1022.0,
}

_D_ELECTRON_GROUP = {
    "K": 1,
    "Ca": 2,
    "Sc": 3,
    "Ti": 4,
    "V": 5,
    "Cr": 6,
    "Mn": 7,
    "Fe": 8,
    "Co": 9,
    "Ni": 10,
    "Cu": 11,
    "Zn": 12,
}


def _load_input_file() -> None:
    path = Path(INPUT_FILE)
    if not path.is_absolute():
        path = ROOT / path
    if not path.exists():
        raise FileNotFoundError(f"input file not found: {path}")

    values = runpy.run_path(
        str(path),
        init_globals={
            "ROOT": ROOT,
            "Path": Path,
        },
    )
    missing = sorted(key for key in _INPUT_KEYS if key not in values)
    if missing:
        preview = ", ".join(missing[:12])
        if len(missing) > 12:
            preview += f", ... ({len(missing)} total)"
        raise ValueError(f"input file is missing required parameters: {preview}")

    loaded = []
    for key, default in _OPTIONAL_INPUT_DEFAULTS.items():
        globals()[key] = values.get(key, default)
        if key in values:
            loaded.append(key)
    for key in sorted(_INPUT_KEYS):
        value = values[key]
        if key in {"output_dir", "xtls_path"} and value is not None and value != "auto":
            value = _input_path(value, path.parent)
        globals()[key] = value
        loaded.append(key)

    _apply_ion_preset()
    print("input file:", path)
    print("loaded input parameters:", len(loaded))


def _input_path(value, base: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return base / path


def _apply_ion_preset() -> None:
    if ion is None:
        return
    preset_element, charge = _parse_ion(ion)
    preset_d_count = _d_count_from_ion(preset_element, charge)

    global element, n_d_electrons, plot_absolute_energy_offset, output_dir
    if element == "auto":
        element = preset_element
    elif str(element).capitalize() != preset_element:
        raise ValueError(f"input ion {ion!r} implies element {preset_element}, but element={element!r}")

    if n_d_electrons == "auto":
        supported = _supported_xas_d_counts(preset_element)
        if preset_d_count not in supported:
            supported_text = ", ".join(f"d{count}" for count in supported) or "none"
            raise ValueError(
                f"{ion} gives nominal d{preset_d_count}, but Appendix-A 3d XAS pairs "
                f"available for {preset_element} are: {supported_text}. "
                "Set n_d_electrons manually if you intentionally want one of those entries."
            )
        n_d_electrons = preset_d_count
    elif int(n_d_electrons) != preset_d_count:
        print(
            "ion warning:",
            f"{ion} nominal d count is d{preset_d_count}, but n_d_electrons={n_d_electrons}",
        )

    if plot_absolute_energy_offset == "auto":
        plot_absolute_energy_offset = _ION_EDGE_L3_EV[preset_element]

    if output_dir == "auto":
        name = case_name or f"{preset_element}_d{preset_d_count}"
        output_dir = ROOT / "outputs" / str(name)

    print(
        "ion preset:",
        f"{ion} -> element={element}, d{n_d_electrons},",
        f"L3 plot offset~{plot_absolute_energy_offset} eV",
    )


def _parse_ion(text: str) -> tuple[str, int]:
    import re

    match = re.fullmatch(r"\s*([A-Za-z]{1,2})\s*([0-9]+)\s*\+\s*", str(text))
    if not match:
        raise ValueError('ion must look like "Fe2+", "Co4+", "Ni2+", etc.')
    element_symbol = match.group(1).capitalize()
    charge = int(match.group(2))
    if element_symbol not in _D_ELECTRON_GROUP:
        raise ValueError(f"no automatic d-count preset for ion {text!r}")
    if element_symbol not in _ION_EDGE_L3_EV:
        raise ValueError(f"no automatic L-edge offset preset for ion {text!r}")
    return element_symbol, charge


def _d_count_from_ion(element_symbol: str, charge: int) -> int:
    d_count = _D_ELECTRON_GROUP[element_symbol] - charge
    if not 0 <= d_count <= 10:
        raise ValueError(f"nominal d count for {element_symbol}{charge}+ is d{d_count}, outside 0..10")
    return d_count


def _supported_xas_d_counts(element_symbol: str) -> tuple[int, ...]:
    available = set(available_appendix_a_3d(element_symbol))
    counts = []
    for _element, p_electrons, d_electrons in available:
        if p_electrons != 6:
            continue
        if (element_symbol, 5, d_electrons + 1) in available:
            counts.append(d_electrons)
    return tuple(sorted(counts))


def calculate_spectrum() -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, object]]:
    effective_max_ligand_holes = _effective_max_ligand_holes()
    ligand_positions = _ligand_positions()
    h_cry_xyz, _, _, _ = crystal_field(
        ligand_positions,
        ten_dq,
        r2,
        r4,
    )
    h_crystal = spin_expand_orbital_matrix(h_cry_xyz)
    h_hybridization, hybridization_label = _build_hybridization_matrix(ligand_positions)

    initial_basis = p_ct_xas_initial_basis(
        n_d_electrons,
        max_ligand_holes=effective_max_ligand_holes,
    )
    final_basis = p_ct_xas_final_basis(
        n_d_electrons,
        max_ligand_holes=effective_max_ligand_holes,
    )
    initial_slater_by_hole = _configuration_slater_entries(
        p_electrons=6,
        first_d=n_d_electrons,
        max_holes=effective_max_ligand_holes,
    )
    final_slater_by_hole = _configuration_slater_entries(
        p_electrons=5,
        first_d=n_d_electrons + 1,
        max_holes=effective_max_ligand_holes,
    )

    h_initial = _build_initial_hamiltonian(
        initial_basis,
        h_crystal,
        h_hybridization,
        initial_slater_by_hole,
    )
    h_final = _build_final_hamiltonian(
        final_basis,
        h_crystal,
        h_hybridization,
        final_slater_by_hole,
    )

    n_spectrum_states = max(1, n_initial_states)
    n_analysis_states = max(n_spectrum_states, n_analyzed_states)
    e_initial, v_initial = lowest_eigenpairs(h_initial, k=n_analysis_states)
    e_initial_spectrum = e_initial[:n_spectrum_states]
    v_initial_spectrum = v_initial[:, :n_spectrum_states]
    weights = _thermal_weights(e_initial_spectrum, temperature_kelvin)
    e_final_analysis, v_final_analysis = lowest_eigenpairs(h_final, k=max(1, n_analyzed_states))
    # `lowest_eigenpairs` returns ascending eigenvalues, so the first analyzed
    # final state is already the final-state ground state.
    onset = float(e_final_analysis[0] - e_initial[0])
    relative_energy = _energy_grid()
    transition_energy = onset + relative_energy
    display_energy = relative_energy + energy_shift
    lorentzian_profile = interpolate_width(
        relative_energy,
        lorentzian_hwhm,
        lorentzian_hwhm_points,
    )
    gaussian_sigma_profile = _gaussian_sigma_profile(relative_energy)
    exact_final = None
    if spectrum_method == "exact":
        print("spectrum_method=exact: full diagonalizing final Hamiltonian...")
        exact_final = _full_diagonalize(h_final)

    curves = {}
    for polarization in ("x", "y", "z"):
        dipole = one_body_between_bases_sparse(
            final_basis,
            initial_basis,
            p_ct_dipole_matrix(polarization),
        )
        curve = np.zeros_like(transition_energy)
        for state_idx, weight in enumerate(weights):
            if weight == 0.0:
                continue
            start_vector = dipole @ v_initial_spectrum[:, state_idx]
            if spectrum_method == "exact":
                curve += _exact_transition_spectrum(
                    exact_final,
                    start_vector,
                    transition_energy,
                    initial_energy=float(e_initial_spectrum[state_idx]),
                    broadening=lorentzian_profile,
                    weight=float(weight),
                )
            else:
                curve += weight * continued_fraction_spectrum(
                    h_final,
                    start_vector,
                    transition_energy,
                    reference_energy=float(e_initial_spectrum[state_idx]),
                    broadening=lorentzian_profile,
                    max_iter=n_recursion,
                    normalize=False,
                )
        curve = variable_gaussian_broaden(curve, display_energy, gaussian_sigma_profile)
        curves[polarization] = curve

    curves["iso"] = (curves["x"] + curves["y"] + curves["z"]) / 3.0
    curves["ld"] = (curves["x"] + curves["y"]) / 2.0 - curves["z"]
    if make_experimental_geometry_curves:
        curves.update(_experimental_geometry_curves(curves))
    curves = _normalize_columns(curves, display_energy, normalize)

    metadata = {
        "initial_basis_size": float(len(initial_basis)),
        "final_basis_size": float(len(final_basis)),
        "initial_ground_energy": float(e_initial[0]),
        "final_ground_energy": float(e_final_analysis[0]),
        "onset": float(onset),
        "requested_max_ligand_holes": int(max_ligand_holes),
        "effective_max_ligand_holes": int(effective_max_ligand_holes),
        "ligand_hole_truncation_reason": _ligand_hole_truncation_reason(effective_max_ligand_holes),
        "initial_states": int(n_spectrum_states),
        "analyzed_states": int(min(n_analyzed_states, len(e_initial))),
        "spectrum_method": spectrum_method,
        "hybridization_mode": hybridization_label,
        "broadening": _broadening_metadata(relative_energy, lorentzian_profile, gaussian_sigma_profile),
        "experimental_geometry": _experimental_geometry_metadata(),
        "configuration_energies": _configuration_energy_rows(effective_max_ligand_holes),
        "state_analysis": {
            "initial": _analyze_eigenstates(
                "initial",
                initial_basis,
                e_initial[:n_analyzed_states],
                v_initial[:, :n_analyzed_states],
                p_electrons=6,
                first_d_electrons=n_d_electrons,
                max_holes=effective_max_ligand_holes,
                thermal_weights=_thermal_weights(e_initial[:n_analyzed_states], temperature_kelvin),
            ),
            "final": _analyze_eigenstates(
                "final",
                final_basis,
                e_final_analysis[:n_analyzed_states],
                v_final_analysis[:, :n_analyzed_states],
                p_electrons=5,
                first_d_electrons=n_d_electrons + 1,
                max_holes=effective_max_ligand_holes,
                thermal_weights=None,
            ),
        },
        "initial_slater": {hole: _slater_as_dict(entry) for hole, entry in initial_slater_by_hole.items()},
        "final_slater": {hole: _slater_as_dict(entry) for hole, entry in final_slater_by_hole.items()},
    }
    return display_energy, curves, metadata


def _build_initial_hamiltonian(basis, h_crystal, h_hybridization, slater_by_hole):
    h_one = p_ct_one_body_matrix(h_d=h_crystal, h_hybridization=h_hybridization)
    hamiltonian = one_body_sparse(basis, h_one)
    hamiltonian = hamiltonian + diagonal_sparse(
        p_ct_charge_transfer_energy_values(
            basis,
            delta,
            u_charge_transfer,
            core_hole_potential=core_hole_potential,
        )
    )
    for holes, slater in slater_by_hole.items():
        sector_basis = _ligand_hole_sector_basis(basis, holes)
        h_so = p_ct_one_body_matrix(h_d=d_spin_orbit_matrix(slater.zeta_d))
        hamiltonian = hamiltonian + _embed_sector_matrix(
            basis,
            sector_basis,
            one_body_sparse(sector_basis, h_so),
        )
        v_pd = np.zeros((16, 16, 16, 16), dtype=complex)
        v_pd[6:16, 6:16, 6:16, 6:16] = d_shell_slater_tensor(
            0.0,
            slater.fdd2,
            slater.fdd4,
        )
        hamiltonian = hamiltonian + _embed_sector_matrix(
            basis,
            sector_basis,
            two_body_sparse(sector_basis, v_pd),
        )
    return hamiltonian


def _build_final_hamiltonian(basis, h_crystal, h_hybridization, slater_by_hole):
    h_one = p_ct_one_body_matrix(
        h_d=h_crystal,
        h_hybridization=h_hybridization,
    )
    hamiltonian = one_body_sparse(basis, h_one)
    hamiltonian = hamiltonian + diagonal_sparse(
        p_ct_charge_transfer_energy_values(
            basis,
            delta,
            u_charge_transfer,
            core_hole_potential=core_hole_potential,
        )
    )
    for holes, slater in slater_by_hole.items():
        sector_basis = _ligand_hole_sector_basis(basis, holes)
        h_so = p_ct_one_body_matrix(
            h_p=p_spin_orbit_matrix(slater.zeta_2p or 0.0),
            h_d=d_spin_orbit_matrix(slater.zeta_d),
        )
        hamiltonian = hamiltonian + _embed_sector_matrix(
            basis,
            sector_basis,
            one_body_sparse(sector_basis, h_so),
        )
        v_pd = pd_shell_slater_tensor(
            fdd2=slater.fdd2,
            fdd4=slater.fdd4,
            fpd2=slater.fpd2 or 0.0,
            gpd1=slater.gpd1 or 0.0,
            gpd3=slater.gpd3 or 0.0,
        )
        hamiltonian = hamiltonian + _embed_sector_matrix(
            basis,
            sector_basis,
            two_body_sparse(sector_basis, v_pd),
        )
    return hamiltonian


def _configuration_slater_entries(p_electrons: int, first_d: int, max_holes: int):
    entries = {}
    for holes in range(max_holes + 1):
        entries[holes] = _scale_multipoles(
            get_appendix_a_3d(element, p_electrons, first_d + holes)
        )
    return entries


def _validate_parameters() -> None:
    if element == "auto" or n_d_electrons == "auto":
        raise ValueError('element/n_d_electrons can be "auto" only when a valid ion preset is provided')
    if output_dir == "auto":
        raise ValueError('output_dir can be "auto" only when a valid ion preset is provided')
    if plot_absolute_energy_offset == "auto":
        raise ValueError('plot_absolute_energy_offset can be "auto" only when a valid ion preset is provided')
    if not element:
        raise ValueError("element is required")
    if not 0 <= n_d_electrons <= 10:
        raise ValueError("n_d_electrons must be between 0 and 10")
    if max_ligand_holes < 0:
        raise ValueError("max_ligand_holes must be non-negative")
    allowed_geometries = {
        "tetrahedral",
        "octahedral",
        "square_planar",
        "square_pyramidal",
        "trigonal_bipyramidal",
        "custom_xyz",
        "custom_spherical",
    }
    if coordination_geometry not in allowed_geometries:
        raise ValueError(f"coordination_geometry must be one of {sorted(allowed_geometries)}")
    if coordination_geometry == "custom_xyz" and not ligand_positions_xyz:
        raise ValueError("ligand_positions_xyz is required when coordination_geometry='custom_xyz'")
    if coordination_geometry == "custom_spherical" and not ligand_positions_spherical:
        raise ValueError(
            "ligand_positions_spherical is required when coordination_geometry='custom_spherical'"
        )
    if coordination_geometry not in {"custom_xyz", "custom_spherical"} and ligand_radius <= 0.0:
        raise ValueError("ligand_radius must be positive")
    if energy_step <= 0.0:
        raise ValueError("energy_step must be positive")
    if energy_min >= energy_max:
        raise ValueError("energy_min must be smaller than energy_max")
    if plot_relative_energy_min >= plot_relative_energy_max:
        raise ValueError("plot_relative_energy_min must be smaller than plot_relative_energy_max")
    if normalize not in {"max", "area", "none"}:
        raise ValueError('normalize must be "max", "area", or "none"')
    if spectrum_method not in {"lanczos", "exact"}:
        raise ValueError('spectrum_method must be "lanczos" or "exact"')
    if inplane_curve not in {"ab", "x", "y"}:
        raise ValueError('inplane_curve must be "ab", "x", or "y"')
    for name, value in (
        ("lorentzian_hwhm", lorentzian_hwhm),
        ("gaussian_sigma", gaussian_sigma),
        ("gaussian_hwhm", gaussian_hwhm),
    ):
        if value < 0.0:
            raise ValueError(f"{name} must be non-negative")
    for _energy, width in lorentzian_hwhm_points:
        if width < 0.0:
            raise ValueError("lorentzian_hwhm_points widths must be non-negative")
    for _energy, width in gaussian_hwhm_points:
        if width < 0.0:
            raise ValueError("gaussian_hwhm_points widths must be non-negative")
    get_appendix_a_3d(element, 6, n_d_electrons)
    get_appendix_a_3d(element, 5, n_d_electrons + 1)


def _estimate_basis_sizes(effective_max_holes: int) -> tuple[int, int]:
    initial_size = 0
    final_size = 0
    for holes in range(effective_max_holes + 1):
        initial_size += comb(10, n_d_electrons + holes) * comb(10, holes)
        final_size += 6 * comb(10, n_d_electrons + 1 + holes) * comb(10, holes)
    return initial_size, final_size


def _effective_max_ligand_holes() -> int:
    if max_ligand_holes < 0:
        raise ValueError("max_ligand_holes must be 0 or larger")
    effective = -1
    for holes in range(max_ligand_holes + 1):
        initial_d = n_d_electrons + holes
        final_d = n_d_electrons + 1 + holes
        if initial_d > 10 or final_d > 10:
            break
        try:
            get_appendix_a_3d(element, 6, initial_d)
            get_appendix_a_3d(element, 5, final_d)
        except KeyError:
            break
        effective = holes
    if effective < 0:
        get_appendix_a_3d(element, 6, n_d_electrons)
        get_appendix_a_3d(element, 5, n_d_electrons + 1)
    return effective


def _ligand_hole_truncation_reason(effective_max_holes: int) -> str:
    if effective_max_holes >= max_ligand_holes:
        return "none"
    next_hole = effective_max_holes + 1
    initial_d = n_d_electrons + next_hole
    final_d = n_d_electrons + 1 + next_hole
    if initial_d > 10 or final_d > 10:
        return f"d-shell electron limit at L{next_hole}"
    missing = []
    for label, p_electrons, d_electrons in (
        ("initial", 6, initial_d),
        ("final", 5, final_d),
    ):
        try:
            get_appendix_a_3d(element, p_electrons, d_electrons)
        except KeyError:
            missing.append(f"{label} 2p{p_electrons} 3d{d_electrons}")
    if missing:
        return "missing Appendix-A Slater entries for " + ", ".join(missing)
    return f"L{next_hole} is unavailable"


def _ligand_hole_sector_basis(basis, holes: int):
    states = [state for state in basis.states if p_ct_ligand_hole_count(state) == holes]
    return basis.__class__.from_states(basis.n_orbitals, basis.n_electrons, states)


def _embed_sector_matrix(global_basis, sector_basis, matrix):
    from scipy import sparse

    coo = matrix.tocoo()
    sector_indices = np.array([global_basis.index[state] for state in sector_basis.states])
    rows = sector_indices[coo.row]
    cols = sector_indices[coo.col]
    return sparse.coo_matrix(
        (coo.data, (rows, cols)),
        shape=(len(global_basis), len(global_basis)),
    ).tocsr()


def _build_hybridization_matrix(ligand_positions: np.ndarray):
    if hybridization_mode == "geometry":
        h_hyb_orbital, _, _, _ = hybridization(
            ligand_positions,
            pd_sigma,
            pd_ratio,
            d_ref,
        )
        return d_ligand_hybridization_from_orbital_matrix(h_hyb_orbital), "geometry"
    if hybridization_mode == "scalar":
        return d_ligand_hybridization_matrix(hopping=hopping), "scalar"
    raise ValueError('hybridization_mode must be "geometry" or "scalar"')


def _thermal_weights(energies: np.ndarray, temperature_k: float) -> np.ndarray:
    weights = np.zeros(len(energies), dtype=float)
    if len(energies) == 0:
        return weights
    if temperature_k <= 0.0 or len(energies) == 1:
        weights[0] = 1.0
        return weights
    boltzmann_ev_per_k = 8.617333262145e-5
    beta_energy = boltzmann_ev_per_k * temperature_k
    shifted = np.asarray(energies) - float(energies[0])
    weights = np.exp(-shifted / beta_energy)
    return weights / np.sum(weights)


def _energy_grid() -> np.ndarray:
    if energy_step <= 0.0:
        raise ValueError("energy_step must be positive")
    count = int(np.floor((energy_max - energy_min) / energy_step + 0.5)) + 1
    if count <= 1:
        raise ValueError("energy range must contain at least two points")
    return energy_min + energy_step * np.arange(count)


def _gaussian_sigma_profile(energy: np.ndarray) -> float | np.ndarray:
    if len(gaussian_hwhm_points) > 0:
        return np.asarray(interpolate_width(energy, 0.0, gaussian_hwhm_points)) * GAUSSIAN_HWHM_TO_SIGMA
    if gaussian_sigma > 0.0:
        return float(gaussian_sigma)
    if gaussian_hwhm > 0.0:
        return float(gaussian_hwhm * GAUSSIAN_HWHM_TO_SIGMA)
    return 0.0


def _broadening_metadata(
    energy: np.ndarray,
    lorentzian_profile: float | np.ndarray,
    gaussian_sigma_profile: float | np.ndarray,
) -> dict[str, object]:
    return {
        "lorentzian_hwhm": _width_summary(lorentzian_profile),
        "gaussian_sigma": _width_summary(gaussian_sigma_profile),
        "gaussian_hwhm_equivalent": _width_summary(
            np.asarray(gaussian_sigma_profile) / GAUSSIAN_HWHM_TO_SIGMA
        ),
        "lorentzian_hwhm_points": [list(point) for point in lorentzian_hwhm_points],
        "gaussian_hwhm_points": [list(point) for point in gaussian_hwhm_points],
        "energy_min": float(np.min(energy)),
        "energy_max": float(np.max(energy)),
    }


def _width_summary(width: float | np.ndarray) -> dict[str, float | bool]:
    values = np.asarray(width, dtype=float)
    if values.ndim == 0:
        value = float(values)
        return {"energy_dependent": False, "min": value, "max": value, "mean": value}
    return {
        "energy_dependent": True,
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
    }


def _format_width_summary(summary: object) -> str:
    assert isinstance(summary, dict)
    if summary.get("energy_dependent"):
        return f"{summary['min']:.4g}-{summary['max']:.4g} eV"
    return f"{summary['mean']:.4g} eV"


def _full_diagonalize(matrix):
    dense = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
    return np.linalg.eigh(dense)


def _exact_transition_spectrum(
    exact_final,
    start_vector: np.ndarray,
    transition_energy: np.ndarray,
    *,
    initial_energy: float,
    broadening: float | np.ndarray,
    weight: float,
    chunk_size: int = 512,
) -> np.ndarray:
    final_energies, final_vectors = exact_final
    amplitudes = final_vectors.conjugate().T @ start_vector
    strengths = weight * np.abs(amplitudes) ** 2
    spectrum = np.zeros_like(transition_energy, dtype=float)

    eta = np.asarray(broadening, dtype=float)
    if eta.ndim == 0:
        eta_value = float(eta)
        if eta_value <= 0.0:
            raise ValueError("broadening must be positive for exact spectrum")
        for start in range(0, len(final_energies), chunk_size):
            stop = min(start + chunk_size, len(final_energies))
            delta = transition_energy[:, None] + initial_energy - final_energies[None, start:stop]
            spectrum += np.sum(
                strengths[None, start:stop] * eta_value / np.pi / (delta**2 + eta_value**2),
                axis=1,
            )
        return spectrum

    if eta.shape != transition_energy.shape:
        raise ValueError("array broadening must have the same shape as energy")
    if np.any(eta <= 0.0):
        raise ValueError("broadening must be positive for exact spectrum")
    eta_column = eta[:, None]
    for start in range(0, len(final_energies), chunk_size):
        stop = min(start + chunk_size, len(final_energies))
        delta = transition_energy[:, None] + initial_energy - final_energies[None, start:stop]
        spectrum += np.sum(
            strengths[None, start:stop] * eta_column / np.pi / (delta**2 + eta_column**2),
            axis=1,
        )
    return spectrum


def _scale_multipoles(entry):
    return replace(
        entry,
        zeta_d=entry.zeta_d * so3d_scale,
        fdd2=entry.fdd2 * fdd2_scale,
        fdd4=entry.fdd4 * fdd4_scale,
        zeta_2p=None if entry.zeta_2p is None else entry.zeta_2p * so2p_scale,
        fpd2=None if entry.fpd2 is None else entry.fpd2 * fpd2_scale,
        gpd1=None if entry.gpd1 is None else entry.gpd1 * gpd1_scale,
        gpd3=None if entry.gpd3 is None else entry.gpd3 * gpd3_scale,
    )


def _configuration_energy_rows(max_holes: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for holes in range(max_holes + 1):
        rows.append(
            {
                "label": f"#i{holes + 1}",
                "state": "initial",
                "configuration": _configuration_label(6, n_d_electrons + holes, holes),
                "energy_eV": _configuration_energy(holes, core_holes=0),
            }
        )
    for holes in range(max_holes + 1):
        rows.append(
            {
                "label": f"#f{holes + 1}",
                "state": "final",
                "configuration": _configuration_label(5, n_d_electrons + 1 + holes, holes),
                "energy_eV": _configuration_energy(holes, core_holes=1),
            }
        )
    return rows


def _configuration_energy(holes: int, core_holes: int) -> float:
    return sector_energy(
        ligand_holes=holes,
        core_holes=core_holes,
        delta=delta,
        u_charge_transfer=u_charge_transfer,
        core_hole_potential=core_hole_potential,
    )


def _configuration_label(p_electrons: int, d_electrons: int, holes: int) -> str:
    label = f"2p{p_electrons} 3d{d_electrons}"
    if holes == 1:
        return f"{label} L"
    if holes > 1:
        return f"{label} L{holes}"
    return label


def _slater_as_dict(entry) -> dict[str, float | None | str | int]:
    return {
        "element": entry.element,
        "p_electrons": entry.p_electrons,
        "d_electrons": entry.d_electrons,
        "r2": entry.r2,
        "r4": entry.r4,
        "zeta_d": entry.zeta_d,
        "fdd2": entry.fdd2,
        "fdd4": entry.fdd4,
        "zeta_2p": entry.zeta_2p,
        "fpd2": entry.fpd2,
        "gpd1": entry.gpd1,
        "gpd3": entry.gpd3,
    }


def _analyze_eigenstates(
    state_kind: str,
    basis,
    energies: np.ndarray,
    vectors: np.ndarray,
    *,
    p_electrons: int,
    first_d_electrons: int,
    max_holes: int,
    thermal_weights: np.ndarray | None,
) -> list[dict[str, object]]:
    if len(energies) == 0:
        return []
    labels = p_ct_spin_orbital_labels()
    moment_ops = _d_shell_moment_operators(basis)
    ground = float(energies[0])
    rows: list[dict[str, object]] = []
    n_states = min(len(energies), vectors.shape[1])
    for state_index in range(n_states):
        probabilities = np.abs(np.asarray(vectors[:, state_index])) ** 2
        total = float(np.sum(probabilities))
        if total > 0.0:
            probabilities = probabilities / total
        occupations = _orbital_occupations(basis, probabilities)
        rows.append(
            {
                "kind": state_kind,
                "index": state_index + 1,
                "energy_eV": float(energies[state_index]),
                "energy_relative_eV": float(energies[state_index] - ground),
                "energy_relative_meV": float(1000.0 * (energies[state_index] - ground)),
                "thermal_weight": (
                    None
                    if thermal_weights is None or state_index >= len(thermal_weights)
                    else float(thermal_weights[state_index])
                ),
                "configuration_weights_percent": _configuration_weights_percent(
                    basis,
                    probabilities,
                    p_electrons=p_electrons,
                    first_d_electrons=first_d_electrons,
                    max_holes=max_holes,
                ),
                "shell_occupations": _shell_occupations(occupations),
                "moments_3d": _moment_expectations(vectors[:, state_index], moment_ops),
                "orbital_occupations": [
                    {"orbital": labels[orbital], "n": float(value)}
                    for orbital, value in enumerate(occupations)
                ],
            }
        )
    return rows


def _orbital_occupations(basis, probabilities: np.ndarray) -> np.ndarray:
    occupations = np.zeros(basis.n_orbitals, dtype=float)
    for basis_index, state in enumerate(basis.states):
        weight = float(probabilities[basis_index])
        if weight == 0.0:
            continue
        for orbital in range(basis.n_orbitals):
            if (state >> orbital) & 1:
                occupations[orbital] += weight
    return occupations


def _shell_occupations(occupations: np.ndarray) -> dict[str, float]:
    p_total = float(np.sum(occupations[0:6]))
    d_total = float(np.sum(occupations[6:16]))
    ligand_total = float(np.sum(occupations[16:26]))
    return {
        "p_electrons": p_total,
        "p_holes": 6.0 - p_total,
        "d_electrons": d_total,
        "ligand_electrons": ligand_total,
        "ligand_holes": 10.0 - ligand_total,
    }


def _configuration_weights_percent(
    basis,
    probabilities: np.ndarray,
    *,
    p_electrons: int,
    first_d_electrons: int,
    max_holes: int,
) -> dict[str, float]:
    weights = {
        _configuration_label(p_electrons, first_d_electrons + holes, holes): 0.0
        for holes in range(max_holes + 1)
    }
    for basis_index, state in enumerate(basis.states):
        holes = p_ct_ligand_hole_count(state)
        label = _configuration_label(p_electrons, first_d_electrons + holes, holes)
        weights[label] = weights.get(label, 0.0) + 100.0 * float(probabilities[basis_index])
    return weights


def _compact_config_percent(row: dict[str, object]) -> str:
    weights = row["configuration_weights_percent"]
    assert isinstance(weights, dict)
    pieces = []
    for label, percent in weights.items():
        ligand_part = "L0"
        if " L2" in label:
            ligand_part = "L2"
        elif " L" in label:
            ligand_part = "L1"
        pieces.append(f"{ligand_part}:{float(percent):.1f}%")
    return ", ".join(pieces)


def _d_shell_moment_operators(basis) -> dict[str, object]:
    identity_orbital = np.eye(5, dtype=complex)
    identity_spin = np.eye(2, dtype=complex)
    lx_orb, ly_orb, lz_orb = d_angular_momentum_matrices()
    sx_spin, sy_spin, sz_spin = spin_matrices()

    d_ops = {
        "Sx": np.kron(identity_orbital, sx_spin),
        "Sy": np.kron(identity_orbital, sy_spin),
        "Sz": np.kron(identity_orbital, sz_spin),
        "Lx": np.kron(lx_orb, identity_spin),
        "Ly": np.kron(ly_orb, identity_spin),
        "Lz": np.kron(lz_orb, identity_spin),
    }
    d_ops["Jx"] = d_ops["Lx"] + d_ops["Sx"]
    d_ops["Jy"] = d_ops["Ly"] + d_ops["Sy"]
    d_ops["Jz"] = d_ops["Lz"] + d_ops["Sz"]
    d_ops["Mx"] = -(2.0 * d_ops["Sx"] + d_ops["Lx"])
    d_ops["My"] = -(2.0 * d_ops["Sy"] + d_ops["Ly"])
    d_ops["Mz"] = -(2.0 * d_ops["Sz"] + d_ops["Lz"])

    return {name: one_body_sparse(basis, _embed_d_operator(matrix)) for name, matrix in d_ops.items()}


def _embed_d_operator(d_matrix: np.ndarray) -> np.ndarray:
    matrix = np.zeros((26, 26), dtype=complex)
    matrix[6:16, 6:16] = d_matrix
    return matrix


def _moment_expectations(vector: np.ndarray, ops: dict[str, object]) -> dict[str, object]:
    axes = ("x", "y", "z")
    moments: dict[str, object] = {}
    for prefix in ("S", "L", "J", "M"):
        components = {
            axis: _expect_sparse(vector, ops[f"{prefix}{axis}"])
            for axis in axes
        }
        moments[prefix] = components

    for prefix in ("S", "L", "J"):
        ox = ops[f"{prefix}x"]
        oy = ops[f"{prefix}y"]
        oz = ops[f"{prefix}z"]
        xx = _expect_square_sparse(vector, ox)
        yy = _expect_square_sparse(vector, oy)
        zz = _expect_square_sparse(vector, oz)
        total = xx + yy + zz
        moments[f"{prefix}2"] = total
        moments[f"{prefix}_quadrupole"] = {
            "3x2_minus_total": 3.0 * xx - total,
            "3y2_minus_total": 3.0 * yy - total,
            "3z2_minus_total": 3.0 * zz - total,
            "y2_minus_z2": yy - zz,
            "z2_minus_x2": zz - xx,
            "x2_minus_y2": xx - yy,
            "yz_plus_zy": _expect_sym_product_sparse(vector, oy, oz),
            "zx_plus_xz": _expect_sym_product_sparse(vector, oz, ox),
            "xy_plus_yx": _expect_sym_product_sparse(vector, ox, oy),
        }
    return moments


def _expect_sparse(vector: np.ndarray, operator) -> float:
    value = np.vdot(vector, operator @ vector)
    return float(np.real_if_close(value))


def _expect_square_sparse(vector: np.ndarray, operator) -> float:
    operated = operator @ vector
    value = np.vdot(operated, operated)
    return float(np.real_if_close(value))


def _expect_sym_product_sparse(vector: np.ndarray, left, right) -> float:
    value = np.vdot(vector, left @ (right @ vector)) + np.vdot(vector, right @ (left @ vector))
    return float(np.real_if_close(value))


def _save_parameters(path: Path, metadata: dict[str, object]) -> None:
    path.write_text(
        json.dumps(_parameter_dict(metadata), indent=2, sort_keys=True),
        encoding="utf-8",
    )


def _save_configuration_energies(path: Path, rows: list[dict[str, object]]) -> None:
    lines = ["label state configuration energy_eV"]
    for row in rows:
        lines.append(
            f"{row['label']} {row['state']} \"{row['configuration']}\" {row['energy_eV']:.12g}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _save_state_analysis(path: Path, analysis: dict[str, list[dict[str, object]]]) -> None:
    lines = [
        "Python CT-Lanczos XTLS-like state analysis",
        "  Energies are absolute eigenvalues in eV.",
        "  dE is relative to the lowest state in each sector.",
        "  Occupations are expectation values of spin-orbital number operators.",
        "",
    ]
    for kind, xtls_tag in (("initial", "#i"), ("final", "#f")):
        rows = analysis.get(kind, [])
        lines.append(f" ====== {kind} sector =========")
        if not rows:
            lines.append("no states")
            lines.append("")
            continue
        for row in rows:
            lines.extend(_format_xtls_like_state(row, xtls_tag))
        lines.append("")
    lines.extend(_format_pyxtls_footer(analysis))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _format_xtls_like_state(row: dict[str, object], xtls_tag: str) -> list[str]:
    state_index = int(row["index"])
    energy = float(row["energy_eV"])
    splitting = float(row["energy_relative_meV"])
    line = (
        f"   # {state_index:11d}   E= {energy:20.12f}"
        f"     dE= {splitting:12.6f} meV"
    )
    if row["thermal_weight"] is not None:
        line += f"     thermal weight= {float(row['thermal_weight']):12.7g}"
    lines = [
        f" ====== state {xtls_tag} =========",
        line,
        " -----------------------------------------------------------",
        *_format_moment_block(row),
        " -----------------------------------------------------------",
        " ---------- configuration mixing ---------",
        f"   # {state_index:11d}   E= {energy:20.12f}",
    ]
    config_weights = row["configuration_weights_percent"]
    assert isinstance(config_weights, dict)
    for config_index, (label, percent) in enumerate(config_weights.items(), start=1):
        lines.append(
            f"     CI={xtls_tag}{config_index:<3d} {label:<16s}"
            f"{float(percent):20.12f}       %"
        )

    shell = row["shell_occupations"]
    assert isinstance(shell, dict)
    lines.extend(
        [
            " -----------------------------------------------------------",
            " Shell occupation summary",
            f"  2p electrons = {float(shell['p_electrons']):10.6f}"
            f"    2p holes = {float(shell['p_holes']):10.6f}",
            f"  3d electrons = {float(shell['d_electrons']):10.6f}"
            f"    ligand holes = {float(shell['ligand_holes']):10.6f}",
        ]
    )

    orbital_rows = row["orbital_occupations"]
    assert isinstance(orbital_rows, list)
    occupations = {
        str(item["orbital"]): float(item["n"])
        for item in orbital_rows
        if isinstance(item, dict)
    }
    lines.extend(
        _format_orbit_occupation_block(
            "2p",
            ("x", "y", "z"),
            occupations,
            label_prefix="p",
        )
    )
    lines.extend(
        _format_orbit_occupation_block(
            "3d",
            ("xy", "yz", "zx", "x2-y2", "3z2-r2"),
            occupations,
            label_prefix="d",
        )
    )
    lines.extend(
        _format_orbit_occupation_block(
            "Ld",
            ("xy", "yz", "zx", "x2-y2", "3z2-r2"),
            occupations,
            label_prefix="L",
        )
    )
    lines.extend(_format_number_operator_block(occupations))
    return lines


def _format_moment_block(row: dict[str, object]) -> list[str]:
    moments = row.get("moments_3d")
    if not isinstance(moments, dict):
        return []
    lines = [
        "  Orbit ==>  3d",
        _format_vector_line("S", moments),
        _format_vector_line("L", moments),
        _format_vector_line("J", moments),
        _format_vector_line("M", moments),
        "",
    ]
    for prefix in ("S", "L", "J"):
        lines.extend(_format_square_block(prefix, moments))
    lines.extend(
        [
            " ----------------------------------------------------------",
            ' Total moments over the orbitals specified in "Mag" command.',
            (
                f"<S_2 S_1>= {_moment_scalar(moments, 'S2'):12.6f}"
                f"    <L_2 L_1>= {_moment_scalar(moments, 'L2'):12.6f}"
                f"    <J_2 J_1>= {_moment_scalar(moments, 'J2'):12.6f}"
            ),
        ]
    )
    return lines


def _format_vector_line(prefix: str, moments: dict[str, object]) -> str:
    vector = moments[prefix]
    assert isinstance(vector, dict)
    return (
        f"  {prefix}x= {float(vector['x']):12.6f}"
        f"      {prefix}y= {float(vector['y']):12.6f}"
        f"      {prefix}z= {float(vector['z']):12.6f}"
    )


def _format_square_block(prefix: str, moments: dict[str, object]) -> list[str]:
    quad = moments[f"{prefix}_quadrupole"]
    assert isinstance(quad, dict)
    total = _moment_scalar(moments, f"{prefix}2")
    return [
        f" {prefix}^2 = {total:12.6f}",
        (
            f"  {prefix}y{prefix}z+{prefix}z{prefix}y= {float(quad['yz_plus_zy']):12.6f}"
            f"      {prefix}z{prefix}x+{prefix}x{prefix}z= {float(quad['zx_plus_xz']):12.6f}"
            f"      {prefix}x{prefix}y+{prefix}y{prefix}x= {float(quad['xy_plus_yx']):12.6f}"
        ),
        (
            f"  3{prefix}x^2-{prefix}^2= {float(quad['3x2_minus_total']):12.6f}"
            f"      3{prefix}y^2-{prefix}^2= {float(quad['3y2_minus_total']):12.6f}"
            f"      3{prefix}z^2-{prefix}^2= {float(quad['3z2_minus_total']):12.6f}"
        ),
        (
            f"  {prefix}y^2-{prefix}z^2= {float(quad['y2_minus_z2']):12.6f}"
            f"      {prefix}z^2-{prefix}x^2= {float(quad['z2_minus_x2']):12.6f}"
            f"      {prefix}x^2-{prefix}y^2= {float(quad['x2_minus_y2']):12.6f}"
        ),
        "",
    ]


def _moment_scalar(moments: dict[str, object], key: str) -> float:
    return float(moments[key])


def _format_orbit_occupation_block(
    shell: str,
    orbitals: tuple[str, ...],
    occupations: dict[str, float],
    *,
    label_prefix: str | None = None,
) -> list[str]:
    prefix = shell if label_prefix is None else label_prefix
    lines = [
        f" Orbit => {shell}",
        "                ----- electron occupation ----",
        "   orbital          down        up       total",
    ]
    down_sum = 0.0
    up_sum = 0.0
    for orbital in orbitals:
        down = occupations.get(f"{prefix}_{orbital}_down", 0.0)
        up = occupations.get(f"{prefix}_{orbital}_up", 0.0)
        down_sum += down
        up_sum += up
        lines.append(f"   {orbital:<10s} {down:10.6f} {up:10.6f} {down + up:10.6f}")
    lines.append(f" sum          {down_sum:10.6f} {up_sum:10.6f}")
    lines.append(f" total        {down_sum + up_sum:10.6f}")
    return lines


def _format_number_operator_block(occupations: dict[str, float]) -> list[str]:
    labels = list(p_ct_spin_orbital_labels())
    values = [occupations.get(label, 0.0) for label in labels]
    lines = [
        " -- Expectation value of the number operator --",
        "    order: " + " ".join(labels),
    ]
    for start in range(0, len(values), 6):
        lines.append("    " + " ".join(f"{value:10.7f}" for value in values[start : start + 6]))
    return lines + [""]


def _format_pyxtls_footer(analysis: dict[str, list[dict[str, object]]]) -> list[str]:
    initial = analysis.get("initial", [])
    final = analysis.get("final", [])
    return [
        " additional information for PY-XTLS",
        " MODE",
        " xas",
        " DICRO1",
        "           3",
        " {",
        " x",
        " y",
        " z",
        " }",
        " ENERGY_INITIAL",
        f" {len(initial):11d}",
        *[f" {float(row['energy_eV']):20.12f}" for row in initial],
        " ENERGY_FINAL",
        f" {len(final):11d}",
        *[f" {float(row['energy_eV']):20.12f}" for row in final],
        " NOTE",
        " SPCD spectra are written in the companion *_spectrum.txt file.",
    ]


def _parameter_dict(metadata: dict[str, object]) -> dict[str, object]:
    return {
        "input_file": str(Path(INPUT_FILE).resolve() if Path(INPUT_FILE).is_absolute() else (ROOT / INPUT_FILE).resolve()),
        "case_name": case_name,
        "ion": ion,
        "element": element,
        "n_d_electrons": n_d_electrons,
        "max_ligand_holes": max_ligand_holes,
        "requested_max_ligand_holes": metadata["requested_max_ligand_holes"],
        "effective_max_ligand_holes": metadata["effective_max_ligand_holes"],
        "ligand_hole_truncation_reason": metadata["ligand_hole_truncation_reason"],
        "ten_dq": ten_dq,
        "delta": delta,
        "pd_sigma": pd_sigma,
        "pd_ratio": pd_ratio,
        "d_ref": d_ref,
        "u_charge_transfer": u_charge_transfer,
        "core_hole_potential": core_hole_potential,
        "coordination_geometry": coordination_geometry,
        "ligand_radius": ligand_radius,
        "ligand_angle_offset_deg": ligand_angle_offset_deg,
        "ligand_positions_xyz": globals().get("ligand_positions_xyz"),
        "ligand_positions_spherical": globals().get("ligand_positions_spherical"),
        "r2": r2,
        "r4": r4,
        "fdd2_scale": fdd2_scale,
        "fdd4_scale": fdd4_scale,
        "fpd2_scale": fpd2_scale,
        "gpd1_scale": gpd1_scale,
        "gpd3_scale": gpd3_scale,
        "so3d_scale": so3d_scale,
        "so2p_scale": so2p_scale,
        "hybridization_mode": hybridization_mode,
        "hopping": hopping,
        "n_initial_states": n_initial_states,
        "n_analyzed_states": n_analyzed_states,
        "temperature_kelvin": temperature_kelvin,
        "spectrum_method": spectrum_method,
        "n_recursion": n_recursion,
        "lorentzian_hwhm": lorentzian_hwhm,
        "gaussian_sigma": gaussian_sigma,
        "gaussian_hwhm": gaussian_hwhm,
        "lorentzian_hwhm_points": lorentzian_hwhm_points,
        "gaussian_hwhm_points": gaussian_hwhm_points,
        "energy_min": energy_min,
        "energy_max": energy_max,
        "energy_step": energy_step,
        "energy_shift": energy_shift,
        "normalize": normalize,
        "make_experimental_geometry_curves": make_experimental_geometry_curves,
        "grazing_angle_deg": grazing_angle_deg,
        "inplane_curve": inplane_curve,
        "save_state_analysis": save_state_analysis,
        "estimate_only": estimate_only,
        "plot_use_absolute_energy": plot_use_absolute_energy,
        "plot_absolute_energy_offset": plot_absolute_energy_offset,
        "plot_relative_energy_min": plot_relative_energy_min,
        "plot_relative_energy_max": plot_relative_energy_max,
        "overlay_xtls": overlay_xtls,
        "xtls_path": str(xtls_path),
        "xtls_energy_shift": xtls_energy_shift,
        "xtls_scale": xtls_scale,
        "xtls_iso_column": xtls_iso_column,
        "metadata": {
            "initial_basis_size": metadata["initial_basis_size"],
            "final_basis_size": metadata["final_basis_size"],
            "initial_ground_energy": metadata["initial_ground_energy"],
            "final_ground_energy": metadata["final_ground_energy"],
            "onset": metadata["onset"],
            "requested_max_ligand_holes": metadata["requested_max_ligand_holes"],
            "effective_max_ligand_holes": metadata["effective_max_ligand_holes"],
            "ligand_hole_truncation_reason": metadata["ligand_hole_truncation_reason"],
            "initial_states": metadata["initial_states"],
            "analyzed_states": metadata["analyzed_states"],
            "hybridization_mode": metadata["hybridization_mode"],
            "broadening": metadata["broadening"],
            "experimental_geometry": metadata["experimental_geometry"],
        },
        "configuration_energies": metadata["configuration_energies"],
        "initial_slater": metadata["initial_slater"],
        "final_slater": metadata["final_slater"],
    }


def _ligand_positions() -> np.ndarray:
    if coordination_geometry == "custom_xyz":
        return _xyz_rows_to_spherical(ligand_positions_xyz)
    if coordination_geometry == "custom_spherical":
        return _spherical_rows_to_array(ligand_positions_spherical)
    if coordination_geometry == "octahedral":
        return np.array(
            [
                [ligand_radius, np.pi / 2.0, 0.0],
                [ligand_radius, np.pi / 2.0, np.pi],
                [ligand_radius, np.pi / 2.0, np.pi / 2.0],
                [ligand_radius, np.pi / 2.0, 3.0 * np.pi / 2.0],
                [ligand_radius, 0.0, 0.0],
                [ligand_radius, np.pi, 0.0],
            ],
            dtype=float,
        )
    if coordination_geometry == "square_planar":
        return np.array(
            [
                [ligand_radius, np.pi / 2.0, 0.0],
                [ligand_radius, np.pi / 2.0, np.pi / 2.0],
                [ligand_radius, np.pi / 2.0, np.pi],
                [ligand_radius, np.pi / 2.0, 3.0 * np.pi / 2.0],
            ],
            dtype=float,
        )
    if coordination_geometry == "square_pyramidal":
        return np.array(
            [
                [ligand_radius, np.pi / 2.0, 0.0],
                [ligand_radius, np.pi / 2.0, np.pi / 2.0],
                [ligand_radius, np.pi / 2.0, np.pi],
                [ligand_radius, np.pi / 2.0, 3.0 * np.pi / 2.0],
                [ligand_radius, 0.0, 0.0],
            ],
            dtype=float,
        )
    if coordination_geometry == "trigonal_bipyramidal":
        return np.array(
            [
                [ligand_radius, np.pi / 2.0, 0.0],
                [ligand_radius, np.pi / 2.0, 2.0 * np.pi / 3.0],
                [ligand_radius, np.pi / 2.0, 4.0 * np.pi / 3.0],
                [ligand_radius, 0.0, 0.0],
                [ligand_radius, np.pi, 0.0],
            ],
            dtype=float,
        )
    if coordination_geometry != "tetrahedral":
        raise ValueError(f"unsupported coordination_geometry: {coordination_geometry}")
    angle = (54.7356 + ligand_angle_offset_deg) * np.pi / 180.0
    return np.array(
        [
            [ligand_radius, angle, -45.0 * np.pi / 180.0],
            [ligand_radius, angle, 135.0 * np.pi / 180.0],
            [ligand_radius, np.pi - angle, 45.0 * np.pi / 180.0],
            [ligand_radius, np.pi - angle, -135.0 * np.pi / 180.0],
        ],
        dtype=float,
    )


def _xyz_rows_to_spherical(rows) -> np.ndarray:
    if not rows:
        raise ValueError("ligand_positions_xyz must contain at least one ligand")
    xyz = []
    for row in rows:
        values = list(row.values()) if isinstance(row, dict) else list(row)
        if len(values) == 4:
            values = values[1:]
        if len(values) != 3:
            raise ValueError("each ligand_positions_xyz row must be (x, y, z) or (label, x, y, z)")
        xyz.append([float(values[0]), float(values[1]), float(values[2])])
    return _cartesian_to_spherical(np.asarray(xyz, dtype=float))


def _spherical_rows_to_array(rows) -> np.ndarray:
    if not rows:
        raise ValueError("ligand_positions_spherical must contain at least one ligand")
    spherical = []
    for row in rows:
        values = list(row.values()) if isinstance(row, dict) else list(row)
        if len(values) == 4:
            values = values[1:]
        if len(values) != 3:
            raise ValueError(
                "each ligand_positions_spherical row must be (r, theta_deg, phi_deg) "
                "or (label, r, theta_deg, phi_deg)"
            )
        radius = float(values[0])
        theta = np.deg2rad(float(values[1]))
        phi = np.deg2rad(float(values[2]))
        spherical.append([radius, theta, phi])
    return np.asarray(spherical, dtype=float)


def _cartesian_to_spherical(xyz: np.ndarray) -> np.ndarray:
    radius = np.linalg.norm(xyz, axis=1)
    if np.any(radius <= 0.0):
        raise ValueError("custom ligand positions must not be at the metal center")
    theta = np.arccos(np.clip(xyz[:, 2] / radius, -1.0, 1.0))
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    return np.column_stack([radius, theta, phi])


def _ligand_radius_summary() -> str:
    try:
        ligand_positions = _ligand_positions()
    except Exception:
        if coordination_geometry in {"custom_xyz", "custom_spherical"}:
            return "custom metal-ligand distances"
        return f"metal-ligand={ligand_radius:g} A"
    distances = ligand_positions[:, 0]
    if np.max(distances) - np.min(distances) < 1e-9:
        return f"metal-ligand={float(distances[0]):g} A"
    return f"metal-ligand={float(np.mean(distances)):.4g} A avg ({float(np.min(distances)):.4g}-{float(np.max(distances)):.4g})"


def _scale_tag() -> str:
    scales = (fdd2_scale, fdd4_scale, fpd2_scale, gpd1_scale, gpd3_scale)
    if all(abs(scale - scales[0]) < 1e-12 for scale in scales):
        return f"scale{scales[0]:.3f}"
    return "separate_scales"


def _experimental_geometry_curves(columns: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    inplane = _inplane_intensity(columns)
    angle = np.deg2rad(grazing_angle_deg)
    inplane_weight = float(np.sin(angle) ** 2)
    outofplane_weight = float(np.cos(angle) ** 2)
    c_grazing = inplane_weight * inplane + outofplane_weight * columns["z"]
    return {
        "ab": inplane,
        "c_grazing": c_grazing,
        "ld_grazing": inplane - c_grazing,
        "c_pure": columns["z"],
        "ld_pure": inplane - columns["z"],
    }


def _inplane_intensity(columns: dict[str, np.ndarray]) -> np.ndarray:
    if inplane_curve == "ab":
        return (columns["x"] + columns["y"]) / 2.0
    if inplane_curve == "x":
        return columns["x"]
    if inplane_curve == "y":
        return columns["y"]
    raise ValueError('inplane_curve must be "ab", "x", or "y"')


def _experimental_geometry_metadata() -> dict[str, object]:
    angle = np.deg2rad(grazing_angle_deg)
    return {
        "enabled": make_experimental_geometry_curves,
        "grazing_angle_deg": grazing_angle_deg,
        "grazing_angle_definition": "beam angle from sample ab plane",
        "inplane_curve": inplane_curve,
        "sigma_polarization": f"E//{inplane_curve}; pure in-plane channel",
        "pi_polarization": "E in incidence plane; mixed in-plane and c channel",
        "pi_inplane_weight": float(np.sin(angle) ** 2),
        "pi_outofplane_weight": float(np.cos(angle) ** 2),
        "formula": "I_pi_grazing = sin(angle)^2 * I_inplane + cos(angle)^2 * I_c",
        "ld_formula": "LD_grazing = I_inplane - I_pi_grazing",
    }


def _normalize_columns(columns: dict[str, np.ndarray], energy: np.ndarray, mode: str):
    if mode == "none":
        return columns
    reference = columns["iso"]
    if mode == "max":
        scale = max(float(np.max(np.abs(reference))), 1e-30)
    elif mode == "area":
        scale = max(abs(float(np.trapz(reference, energy))), 1e-30)
    else:
        raise ValueError('normalize must be "max", "area", or "none"')
    return {name: values / scale for name, values in columns.items()}


def _plot(energy: np.ndarray, curves: dict[str, np.ndarray], out_png: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; text output was still saved.")
        return

    plot_energy = energy + plot_absolute_energy_offset if plot_use_absolute_energy else energy
    plot_mask = _plot_energy_mask(energy)

    if plot_iso_only:
        fig, ax = plt.subplots(figsize=(6.0, 3.6))
        ax.plot(
            plot_energy[plot_mask],
            curves["iso"][plot_mask],
            label=f"{element} L2,3 XAS calc",
            color="black",
            linewidth=1.8,
        )
        ax.set_xlabel("Photon energy (eV)" if plot_use_absolute_energy else "Relative photon energy (eV)")
        ax.set_ylabel("Intensity (arb. unit)")
        ax.legend()
        fig.tight_layout()
        if save_png:
            fig.savefig(out_png, dpi=200)
            print("saved:", out_png)
        if show_plot:
            fig.canvas.draw_idle()
            plt.show()
        else:
            plt.close(fig)
        return

    if plot_ld_stacked and "ab" in curves and "c_grazing" in curves and "ld_grazing" in curves:
        fig, ax = plt.subplots(figsize=(5.2, 6.0))
        reference = max(float(np.max(np.abs(curves["ab"][plot_mask]))), 1e-30)
        ax.plot(
            plot_energy[plot_mask],
            (curves["ab"] / reference + plot_xas_offset)[plot_mask],
            label="E // ab",
            color="tab:blue",
            linewidth=1.8,
        )
        ax.plot(
            plot_energy[plot_mask],
            (curves["c_grazing"] / reference + plot_xas_offset)[plot_mask],
            label="E // c",
            color="red",
            linewidth=1.6,
        )
        ax.plot(
            plot_energy[plot_mask],
            (curves["ld_grazing"] / reference * plot_ld_scale + plot_ld_offset)[plot_mask],
            label="LD",
            color="limegreen",
            linewidth=1.7,
        )
        ax.axhline(plot_ld_offset, color="black", linewidth=1.0, linestyle=(0, (4, 4)))
        ax.set_xlabel("Photon energy (eV)" if plot_use_absolute_energy else "Relative photon energy (eV)")
        ax.set_ylabel("XAS Intensity (arb. unit)")
        ax.legend()
        ax.text(0.78, 0.72, "calc", transform=ax.transAxes, fontsize=13, color="black")
        fig.tight_layout()
        if show_cluster_inset:
            _plot_cluster_inset(fig, ax)
        if save_png:
            fig.savefig(out_png, dpi=200)
            print("saved:", out_png)
        if show_plot:
            fig.canvas.draw_idle()
            plt.show()
        else:
            plt.close(fig)
        return

    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
    axes[0].plot(plot_energy[plot_mask], curves["iso"][plot_mask], label="iso", color="black")
    axes[0].plot(plot_energy[plot_mask], curves["x"][plot_mask], label="X", alpha=0.75)
    axes[0].plot(plot_energy[plot_mask], curves["y"][plot_mask], label="Y", alpha=0.75)
    axes[0].plot(plot_energy[plot_mask], curves["z"][plot_mask], label="Z", alpha=0.75)
    if "ab" in curves and "c_grazing" in curves:
        axes[0].plot(plot_energy[plot_mask], curves["ab"][plot_mask], label="E//ab", color="tab:cyan", linewidth=1.2)
        axes[0].plot(
            plot_energy[plot_mask],
            curves["c_grazing"][plot_mask],
            label=f"E//c grazing {grazing_angle_deg:g} deg",
            color="tab:pink",
            linewidth=1.2,
        )
    axes[1].plot(plot_energy[plot_mask], curves["ld"][plot_mask], label="pure LD = in-plane - Z", color="tab:green")
    if "ld_grazing" in curves:
        axes[1].plot(
            plot_energy[plot_mask],
            curves["ld_grazing"][plot_mask],
            label="grazing LD = E//ab - E//c",
            color="tab:purple",
            linewidth=1.4,
        )
    if overlay_xtls and xtls_path is not None and Path(xtls_path).exists():
        _plot_xtls_overlay(axes)
    axes[1].axhline(0, color="0.6", linewidth=0.8)
    axes[0].set_ylabel("Intensity")
    axes[1].set_xlabel("Photon energy (eV)" if plot_use_absolute_energy else "Relative photon energy (eV)")
    axes[1].set_ylabel("LD")
    axes[0].legend()
    axes[1].legend()
    fig.tight_layout()
    if show_cluster_inset:
        _plot_cluster_inset(fig, axes[0])

    if save_png:
        fig.savefig(out_png, dpi=200)
        print("saved:", out_png)
    if show_plot:
        fig.canvas.draw_idle()
        plt.show()
    else:
        plt.close(fig)


def _plot_xtls_overlay(axes) -> None:
    table = _read_table(xtls_path)
    if xtls_iso_column in table:
        axes[0].plot(
            table["energy"] + xtls_energy_shift,
            _normalized(table[xtls_iso_column]) * xtls_scale,
            label="XTLS iso",
            color="tab:red",
            linewidth=1.8,
            alpha=0.75,
        )
    if all(name in table for name in ("x_broadened", "y_broadened", "z_broadened")):
        xtls_ld = (table["x_broadened"] + table["y_broadened"]) / 2.0 - table["z_broadened"]
        axes[1].plot(
            table["energy"] + xtls_energy_shift,
            _normalized(xtls_ld) * xtls_scale,
            label="XTLS LD",
            color="tab:orange",
            linewidth=1.6,
            alpha=0.75,
        )


def _read_table(path: str | Path) -> dict[str, np.ndarray]:
    path = Path(path)
    header = path.read_text(errors="replace").splitlines()[0].strip().split()
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return {name: data[:, idx] for idx, name in enumerate(header)}


def _normalized(values: np.ndarray) -> np.ndarray:
    scale = float(np.max(np.abs(values))) if values.size else 0.0
    if scale <= 0.0:
        return values.copy()
    return values / scale


def _plot_energy_mask(relative_energy: np.ndarray) -> np.ndarray:
    mask = np.ones_like(relative_energy, dtype=bool)
    if plot_relative_energy_min is not None:
        mask &= relative_energy >= plot_relative_energy_min
    if plot_relative_energy_max is not None:
        mask &= relative_energy <= plot_relative_energy_max
    return mask


def _plot_cluster_inset(fig, host_ax) -> None:
    ligand_positions = _ligand_positions()
    ligand_xyz = _spherical_to_cartesian(ligand_positions)
    bbox = host_ax.get_position()
    inset = fig.add_axes(
        [
            bbox.x0 + 0.025 * bbox.width,
            bbox.y0 + 0.46 * bbox.height,
            0.29 * bbox.width,
            0.48 * bbox.height,
        ],
        projection="3d",
    )
    inset.set_zorder(10)
    inset.patch.set_facecolor("white")
    inset.patch.set_alpha(0.86)

    metal = np.zeros(3)
    for ligand in ligand_xyz:
        inset.plot(
            [metal[0], ligand[0]],
            [metal[1], ligand[1]],
            [metal[2], ligand[2]],
            color="0.35",
            linewidth=1.0,
            alpha=0.9,
        )

    inset.scatter([0.0], [0.0], [0.0], s=58, color="tab:red", depthshade=True)
    inset.scatter(
        ligand_xyz[:, 0],
        ligand_xyz[:, 1],
        ligand_xyz[:, 2],
        s=36,
        color="tab:blue",
        depthshade=True,
    )
    inset.text(0.0, 0.0, 0.0, f" {element}", color="tab:red", fontsize=7)

    limit = max(float(np.max(np.abs(ligand_xyz))), ligand_radius) * 1.12
    inset.set_xlim(-limit, limit)
    inset.set_ylim(-limit, limit)
    inset.set_zlim(-limit, limit)
    try:
        inset.set_box_aspect((1, 1, 1))
    except AttributeError:
        pass
    inset.view_init(elev=18, azim=-42)
    inset.set_axis_off()
    inset.text2D(
        0.02,
        0.98,
        (
            f"M-L {np.mean(np.linalg.norm(ligand_xyz, axis=1)):.3f} A\n"
            f"L-M-L {_bond_angle_summary(ligand_xyz)} deg\n"
            f"10Dq {ten_dq:.3g} eV, pdSig {pd_sigma:.3g}"
        ),
        transform=inset.transAxes,
        va="top",
        ha="left",
        fontsize=7,
        color="0.1",
    )


def _spherical_to_cartesian(positions_sph: np.ndarray) -> np.ndarray:
    radius = positions_sph[:, 0]
    theta = positions_sph[:, 1]
    phi = positions_sph[:, 2]
    return np.column_stack(
        [
            radius * np.sin(theta) * np.cos(phi),
            radius * np.sin(theta) * np.sin(phi),
            radius * np.cos(theta),
        ]
    )


def _bond_angle_summary(ligand_xyz: np.ndarray) -> str:
    unit = ligand_xyz / np.linalg.norm(ligand_xyz, axis=1)[:, None]
    angles: list[float] = []
    for i in range(len(unit)):
        for j in range(i + 1, len(unit)):
            cosine = float(np.clip(np.dot(unit[i], unit[j]), -1.0, 1.0))
            angles.append(float(np.degrees(np.arccos(cosine))))
    rounded = sorted({round(angle, 1) for angle in angles})
    if len(rounded) <= 3:
        return "/".join(f"{angle:.1f}" for angle in rounded)
    return f"{min(rounded):.1f}-{max(rounded):.1f}"


if __name__ == "__main__":
    main()

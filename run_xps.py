"""Spyder-friendly charge-transfer cluster XPS runner.

Choose an input file below, open this file in Spyder, and press F5.

This covers the `Mode=XPS` / `Mode=PES` calculation of XTLS: a photoelectron is
removed from the cluster and the remaining `N-1` electron system is left in a
correlated multiplet final state. Three emission channels are available:

    photoemission_shell = "2p"      core-level XPS,      2p6 3d^n -> 2p5 3d^n
    photoemission_shell = "3d"      valence-band XPS,    2p6 3d^n -> 2p6 3d^(n-1)
    photoemission_shell = "ligand"  ligand-band XPS,     removes a ligand electron

Nothing here is specific to one compound: element, valence, coordination
geometry, charge-transfer parameters and hybridization are all read from the
input file.
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
    SlaterEntry,
    annihilation_between_bases_sparse,
    available_appendix_a_3d,
    continued_fraction_spectrum,
    d_ligand_hybridization_from_orbital_matrix,
    d_ligand_hybridization_matrix,
    d_shell_slater_tensor,
    d_spin_orbit_matrix,
    diagonal_sparse,
    get_appendix_a_3d,
    lowest_eigenpairs,
    one_body_sparse,
    p_ct_charge_transfer_energy_values,
    p_ct_core_xps_final_basis,
    p_ct_ligand_hole_count,
    p_ct_one_body_matrix,
    p_ct_photoemission_orbitals,
    p_ct_spin_orbital_labels,
    p_ct_valence_xps_final_basis,
    p_ct_xas_initial_basis,
    p_spin_orbit_matrix,
    pd_shell_slater_tensor,
    sector_energy,
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

INPUT_FILE = ROOT / "inputs" / "xps_NiO.py"


def main() -> None:
    _load_input_file()
    _validate_parameters()
    effective_estimate_holes = _effective_max_ligand_holes()
    estimated_initial, estimated_final = _estimate_basis_sizes(effective_estimate_holes)
    print("parameter check: OK")
    print("element/configuration:", f"{element} d{n_d_electrons}", f"h={effective_estimate_holes}")
    print("photoemission channel:", photoemission_shell, f"({_channel_description()})")
    print("geometry:", coordination_geometry, _ligand_radius_summary())
    print("estimated effective max_ligand_holes:", effective_estimate_holes)
    print("estimated basis size:", f"initial={estimated_initial}", f"final={estimated_final}")
    if estimate_only:
        print("estimate_only=True, skipping the spectrum calculation.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    display_energy, curves, metadata = calculate_spectrum()
    effective_max_ligand_holes = int(metadata["effective_max_ligand_holes"])

    stem = (
        f"{element}_d{n_d_electrons}_XPS_{photoemission_shell}"
        f"_{spectrum_method}_h{effective_max_ligand_holes}_{_scale_tag()}"
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
    print("ligand_ten_dq:", ligand_ten_dq)
    print("delta:", delta)
    print("u_charge_transfer:", u_charge_transfer)
    print("core_hole_potential:", core_hole_potential)
    print("hybridization_mode:", metadata["hybridization_mode"])
    print("spectrum_method:", spectrum_method)
    print("emission channels:", metadata["n_channels"])
    print("initial states:", metadata["initial_states"])
    print("temperature K:", temperature_kelvin)
    print("initial basis size:", metadata["initial_basis_size"])
    print("final basis size:", metadata["final_basis_size"])
    print("initial ground energy:", round(metadata["initial_ground_energy"], 6))
    print("final ground energy:", round(metadata["final_ground_energy"], 6))
    print("onset (lowest binding energy):", round(metadata["onset"], 6))
    print("total spectral weight:", round(metadata["sum_rule_weight"], 6))
    print("  (sum rule: expected", metadata["sum_rule_expected"], "for a filled shell)")
    print("lorentzian_hwhm:", _format_width_summary(metadata["broadening"]["lorentzian_hwhm"]))
    print("gaussian_sigma:", _format_width_summary(metadata["broadening"]["gaussian_sigma"]))
    print("energy range:", energy_min, energy_max, energy_step)
    print("configuration energies:")
    for row in metadata["configuration_energies"]:
        print(f"  {row['label']} {row['configuration']}: {row['energy_eV']:.6g} eV")
    print("initial state composition:")
    for row in metadata["state_analysis"]["initial"][:3]:
        print(f"  state {row['index']}: {row['energy_relative_meV']:.3f} meV, {row['composition']}")
    print("final state composition (lowest binding energy):")
    for row in metadata["state_analysis"]["final"][:3]:
        print(f"  state {row['index']}: {row['energy_relative_meV']:.3f} meV, {row['composition']}")
    if save_txt:
        print("saved:", out_txt)
    if save_parameters:
        print("saved:", out_params)
    if save_configuration_energies:
        print("saved:", out_config)
    if save_state_analysis:
        print("saved:", out_analysis)

    _plot(display_energy, curves, out_png)


# ---------------------------------------------------------------------------
# Input handling.

_INPUT_KEYS = {
    "element",
    "n_d_electrons",
    "max_ligand_holes",
    "photoemission_shell",
    "spin_resolved",
    "ten_dq",
    "ligand_ten_dq",
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
    "v_eg",
    "v_t2g",
    "n_initial_states",
    "temperature_kelvin",
    "spectrum_method",
    "n_recursion",
    "n_analyzed_states",
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
    "save_txt",
    "save_png",
    "save_parameters",
    "save_configuration_energies",
    "save_state_analysis",
    "estimate_only",
    "show_plot",
    "plot_binding_energy_axis",
    "plot_use_absolute_energy",
    "plot_absolute_energy_offset",
    "plot_relative_energy_min",
    "plot_relative_energy_max",
    "output_dir",
}

_OPTIONAL_INPUT_DEFAULTS = {
    "case_name": None,
    "ion": None,
    "ligand_positions_xyz": None,
    "ligand_positions_spherical": None,
}

# Approximate 2p3/2 binding energies, used only to put the plot on an absolute
# axis. They do not enter the calculation.
_ION_BINDING_2P32_EV = {
    "K": 294.0,
    "Ca": 347.0,
    "Sc": 399.0,
    "Ti": 459.0,
    "V": 516.0,
    "Cr": 577.0,
    "Mn": 641.0,
    "Fe": 711.0,
    "Co": 781.0,
    "Ni": 854.0,
    "Cu": 933.0,
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
        if key == "output_dir" and value is not None and value != "auto":
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
        supported = _supported_d_counts(preset_element)
        if preset_d_count not in supported:
            supported_text = ", ".join(f"d{count}" for count in supported) or "none"
            raise ValueError(
                f"{ion} gives nominal d{preset_d_count}, but the Appendix-A entries "
                f"available for {preset_element} in the {photoemission_shell} channel are: "
                f"{supported_text}. Set n_d_electrons manually to override."
            )
        n_d_electrons = preset_d_count
    elif int(n_d_electrons) != preset_d_count:
        print("ion warning:", f"{ion} nominal d count is d{preset_d_count}, but n_d_electrons={n_d_electrons}")

    if plot_absolute_energy_offset == "auto":
        plot_absolute_energy_offset = (
            _ION_BINDING_2P32_EV[preset_element] if photoemission_shell == "2p" else 0.0
        )

    if output_dir == "auto":
        name = case_name or f"{preset_element}_d{preset_d_count}_XPS_{photoemission_shell}"
        output_dir = ROOT / "outputs" / str(name)

    print("ion preset:", f"{ion} -> element={element}, d{n_d_electrons}")


def _parse_ion(text: str) -> tuple[str, int]:
    import re

    match = re.fullmatch(r"\s*([A-Za-z]{1,2})\s*([0-9]+)\s*\+\s*", str(text))
    if not match:
        raise ValueError('ion must look like "Fe2+", "Co4+", "Ni2+", etc.')
    element_symbol = match.group(1).capitalize()
    charge = int(match.group(2))
    if element_symbol not in _D_ELECTRON_GROUP:
        raise ValueError(f"no automatic d-count preset for ion {text!r}")
    return element_symbol, charge


def _d_count_from_ion(element_symbol: str, charge: int) -> int:
    d_count = _D_ELECTRON_GROUP[element_symbol] - charge
    if not 0 <= d_count <= 10:
        raise ValueError(f"nominal d count for {element_symbol}{charge}+ is d{d_count}, outside 0..10")
    return d_count


def _supported_d_counts(element_symbol: str) -> tuple[int, ...]:
    """d counts for which both the initial and the final sector have entries."""
    available = set(available_appendix_a_3d(element_symbol))
    counts = []
    for _element, p_electrons, d_electrons in available:
        if p_electrons != 6:
            continue
        if _final_slater_key(d_electrons, 0) in available:
            counts.append(d_electrons)
    return tuple(sorted(counts))


# ---------------------------------------------------------------------------
# Spectroscopy channel bookkeeping.
#
# Every channel is defined by three numbers: how the final sector's d count is
# shifted relative to the initial `3d^n`, how many core holes it carries, and
# which spin-orbitals the photoelectron can come from.


def _is_core_channel() -> bool:
    return photoemission_shell == "2p"


def _final_d_offset() -> int:
    """d-electron offset of the final sector relative to the initial one."""
    return 0 if _is_core_channel() else -1


def _final_core_holes() -> int:
    return 1 if _is_core_channel() else 0


def _final_p_electrons() -> int:
    return 5 if _is_core_channel() else 6


def _final_slater_key(initial_d: int, holes: int) -> tuple[str, int, int]:
    return (element, _final_p_electrons(), initial_d + _final_d_offset() + holes)


def _channel_description() -> str:
    if _is_core_channel():
        return "core-level XPS, 2p6 3d^n -> 2p5 3d^n"
    if photoemission_shell == "3d":
        return "valence-band XPS, 2p6 3d^n -> 2p6 3d^(n-1)"
    return "ligand-band XPS, removes a ligand electron"


def _final_basis(holes: int):
    if _is_core_channel():
        return p_ct_core_xps_final_basis(n_d_electrons, max_ligand_holes=holes)
    return p_ct_valence_xps_final_basis(n_d_electrons, max_ligand_holes=holes)


# ---------------------------------------------------------------------------
# The calculation.


def calculate_spectrum() -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, object]]:
    effective_max_ligand_holes = _effective_max_ligand_holes()
    ligand_positions = _ligand_positions()
    h_cry_xyz, _, _, _ = crystal_field(ligand_positions, ten_dq, r2, r4)
    h_crystal = spin_expand_orbital_matrix(h_cry_xyz)
    if ligand_ten_dq != 0.0:
        h_lig_xyz, _, _, _ = crystal_field(ligand_positions, ligand_ten_dq, r2, r4)
        h_ligand = spin_expand_orbital_matrix(h_lig_xyz)
    else:
        h_ligand = None
    h_hybridization, hybridization_label = _build_hybridization_matrix(ligand_positions)

    initial_basis = p_ct_xas_initial_basis(n_d_electrons, max_ligand_holes=effective_max_ligand_holes)
    final_basis = _final_basis(effective_max_ligand_holes)

    initial_slater_by_hole = _slater_entries(6, n_d_electrons, effective_max_ligand_holes)
    final_slater_by_hole = _slater_entries(
        _final_p_electrons(),
        n_d_electrons + _final_d_offset(),
        effective_max_ligand_holes,
    )

    h_initial = _build_hamiltonian(
        initial_basis,
        h_crystal,
        h_ligand,
        h_hybridization,
        initial_slater_by_hole,
        d_electron_offset=0,
        with_core_hole=False,
    )
    h_final = _build_hamiltonian(
        final_basis,
        h_crystal,
        h_ligand,
        h_hybridization,
        final_slater_by_hole,
        d_electron_offset=_final_d_offset(),
        with_core_hole=_is_core_channel(),
    )

    n_spectrum_states = max(1, n_initial_states)
    n_analysis_states = max(n_spectrum_states, n_analyzed_states)
    e_initial, v_initial = lowest_eigenpairs(h_initial, k=min(n_analysis_states, len(initial_basis) - 1))
    e_initial_spectrum = e_initial[:n_spectrum_states]
    v_initial_spectrum = v_initial[:, :n_spectrum_states]
    weights = _thermal_weights(e_initial_spectrum, temperature_kelvin)

    n_final_analysis = max(1, min(n_analyzed_states, len(final_basis) - 1))
    e_final_analysis, v_final_analysis = lowest_eigenpairs(h_final, k=n_final_analysis)
    onset = float(e_final_analysis[0] - e_initial[0])

    relative_energy = _energy_grid()
    transition_energy = onset + relative_energy
    display_energy = relative_energy + energy_shift
    lorentzian_profile = interpolate_width(relative_energy, lorentzian_hwhm, lorentzian_hwhm_points)
    gaussian_sigma_profile = _gaussian_sigma_profile(relative_energy)

    exact_final = None
    if spectrum_method == "exact":
        print("spectrum_method=exact: full diagonalizing the final Hamiltonian...")
        exact_final = _full_diagonalize(h_final)

    channels = p_ct_photoemission_orbitals(photoemission_shell)
    labels = p_ct_spin_orbital_labels()
    raw_curves: dict[str, np.ndarray] = {"total": np.zeros_like(transition_energy)}
    if spin_resolved:
        raw_curves["spin_up"] = np.zeros_like(transition_energy)
        raw_curves["spin_down"] = np.zeros_like(transition_energy)

    sum_rule_weight = 0.0
    for orbital in channels:
        operator = annihilation_between_bases_sparse(final_basis, initial_basis, orbital)
        curve = np.zeros_like(transition_energy)
        for state_idx, weight in enumerate(weights):
            if weight == 0.0:
                continue
            start_vector = operator @ v_initial_spectrum[:, state_idx]
            sum_rule_weight += float(weight) * float(np.vdot(start_vector, start_vector).real)
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
        raw_curves["total"] = raw_curves["total"] + curve
        if spin_resolved:
            key = "spin_up" if labels[orbital].endswith("_up") else "spin_down"
            raw_curves[key] = raw_curves[key] + curve

    curves = {
        name: variable_gaussian_broaden(values, display_energy, gaussian_sigma_profile)
        for name, values in raw_curves.items()
    }
    curves = _normalize_columns(curves, display_energy, normalize)

    metadata = {
        "photoemission_shell": photoemission_shell,
        "channel_description": _channel_description(),
        "n_channels": len(channels),
        "initial_basis_size": float(len(initial_basis)),
        "final_basis_size": float(len(final_basis)),
        "initial_ground_energy": float(e_initial[0]),
        "final_ground_energy": float(e_final_analysis[0]),
        "onset": onset,
        "sum_rule_weight": sum_rule_weight,
        "sum_rule_expected": len(channels) if _is_core_channel() else n_d_electrons,
        "requested_max_ligand_holes": int(max_ligand_holes),
        "effective_max_ligand_holes": int(effective_max_ligand_holes),
        "ligand_hole_truncation_reason": _ligand_hole_truncation_reason(effective_max_ligand_holes),
        "initial_states": int(n_spectrum_states),
        "spectrum_method": spectrum_method,
        "hybridization_mode": hybridization_label,
        "broadening": _broadening_metadata(relative_energy, lorentzian_profile, gaussian_sigma_profile),
        "configuration_energies": _configuration_energy_rows(effective_max_ligand_holes),
        "state_analysis": {
            "initial": _analyze_states(initial_basis, e_initial, v_initial, n_d_electrons, 6),
            "final": _analyze_states(
                final_basis,
                e_final_analysis,
                v_final_analysis,
                n_d_electrons + _final_d_offset(),
                _final_p_electrons(),
            ),
        },
        "initial_slater": {hole: _slater_as_dict(entry) for hole, entry in initial_slater_by_hole.items()},
        "final_slater": {hole: _slater_as_dict(entry) for hole, entry in final_slater_by_hole.items()},
    }
    return display_energy, curves, metadata


def _build_hamiltonian(
    basis,
    h_crystal,
    h_ligand,
    h_hybridization,
    slater_by_hole,
    *,
    d_electron_offset: int,
    with_core_hole: bool,
):
    """Assemble one sector's Hamiltonian.

    The same routine serves the initial state and every final channel; only the
    sector energy offset, the 2p spin-orbit term and the 2p-3d multipole terms
    depend on whether a core hole is present.
    """
    h_one = p_ct_one_body_matrix(
        h_d=h_crystal,
        h_ligand=h_ligand,
        h_hybridization=h_hybridization,
    )
    hamiltonian = one_body_sparse(basis, h_one)
    hamiltonian = hamiltonian + diagonal_sparse(
        p_ct_charge_transfer_energy_values(
            basis,
            delta,
            u_charge_transfer,
            core_hole_potential=core_hole_potential,
            d_electron_offset=d_electron_offset,
        )
    )
    for holes, slater in slater_by_hole.items():
        sector_basis = _ligand_hole_sector_basis(basis, holes)
        if len(sector_basis) == 0:
            continue
        h_so = p_ct_one_body_matrix(
            h_p=p_spin_orbit_matrix(slater.zeta_2p or 0.0) if with_core_hole else None,
            h_d=d_spin_orbit_matrix(slater.zeta_d),
        )
        hamiltonian = hamiltonian + _embed_sector_matrix(
            basis,
            sector_basis,
            one_body_sparse(sector_basis, h_so),
        )
        if with_core_hole:
            v_pd = pd_shell_slater_tensor(
                fdd2=slater.fdd2,
                fdd4=slater.fdd4,
                fpd2=slater.fpd2 or 0.0,
                gpd1=slater.gpd1 or 0.0,
                gpd3=slater.gpd3 or 0.0,
            )
        else:
            v_pd = np.zeros((16, 16, 16, 16), dtype=complex)
            v_pd[6:16, 6:16, 6:16, 6:16] = d_shell_slater_tensor(0.0, slater.fdd2, slater.fdd4)
        hamiltonian = hamiltonian + _embed_sector_matrix(
            basis,
            sector_basis,
            two_body_sparse(sector_basis, v_pd),
        )
    return hamiltonian


def _slater_entries(p_electrons: int, first_d: int, max_holes: int):
    entries = {}
    for holes in range(max_holes + 1):
        entries[holes] = _scale_multipoles(_slater_entry(p_electrons, first_d + holes))
    return entries


def _slater_entry(p_electrons: int, d_electrons: int):
    """Appendix-A multipoles, with closed d shells filled in.

    Charge-transfer sectors routinely reach `3d0` or `3d10`, which have no
    multiplet structure and therefore no tabulated Slater integrals. XTLS
    X-cards handle this by simply omitting the `Rk`/`Zeta` commands for those
    configurations, so zeros are the faithful equivalent. The 2p spin-orbit
    constant is the exception: it must survive, and it barely depends on the d
    count (for Ni it moves by 0.003 eV from 3d8 to 3d10), so the nearest
    tabulated value is used.
    """
    try:
        return get_appendix_a_3d(element, p_electrons, d_electrons)
    except KeyError:
        if d_electrons not in (0, 10):
            raise
    return SlaterEntry(
        element=element,
        p_electrons=p_electrons,
        d_electrons=d_electrons,
        r2=0.0,
        r4=0.0,
        zeta_d=0.0,
        fdd2=0.0,
        fdd4=0.0,
        zeta_2p=_nearest_zeta_2p(d_electrons) if p_electrons == 5 else None,
        fpd2=0.0 if p_electrons == 5 else None,
        gpd1=0.0 if p_electrons == 5 else None,
        gpd3=0.0 if p_electrons == 5 else None,
    )


def _nearest_zeta_2p(d_electrons: int) -> float:
    candidates = [
        (abs(entry_d - d_electrons), entry_d)
        for _element, p_electrons, entry_d in available_appendix_a_3d(element)
        if p_electrons == 5
    ]
    if not candidates:
        raise KeyError(f"no core-hole Appendix-A entry for {element}, cannot infer zeta_2p")
    _distance, nearest_d = min(candidates)
    zeta = get_appendix_a_3d(element, 5, nearest_d).zeta_2p
    if zeta is None:
        raise KeyError(f"Appendix-A entry {element} 2p5 3d{nearest_d} has no zeta_2p")
    return float(zeta)


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
        h_hyb_orbital, _, _, _ = hybridization(ligand_positions, pd_sigma, pd_ratio, d_ref)
        return d_ligand_hybridization_from_orbital_matrix(h_hyb_orbital), "geometry"
    if hybridization_mode == "symmetry":
        # XTLS `VOh(#sc1 #sc2 Ld 3d) = {V(eg), V(t2g)}`. D_ORBITALS order is
        # (xy, yz, zx, x2-y2, 3z2-r2), so t2g comes first and eg last.
        per_orbital = np.array([v_t2g, v_t2g, v_t2g, v_eg, v_eg], dtype=float)
        return d_ligand_hybridization_matrix(hopping=np.repeat(per_orbital, 2)), "symmetry"
    if hybridization_mode == "scalar":
        return d_ligand_hybridization_matrix(hopping=hopping), "scalar"
    raise ValueError('hybridization_mode must be "geometry", "symmetry", or "scalar"')


# ---------------------------------------------------------------------------
# Geometry, grids and broadening.


def _ligand_positions() -> np.ndarray:
    """Ligand sites as `[radius, theta, phi]` rows, angles in radians."""
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
    xyz = np.asarray(xyz, dtype=float)
    radius = np.linalg.norm(xyz, axis=1)
    if np.any(radius <= 0.0):
        raise ValueError("custom ligand positions must not be at the metal center")
    theta = np.arccos(np.clip(xyz[:, 2] / radius, -1.0, 1.0))
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    return np.column_stack([radius, theta, phi])


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
        spherical.append([float(values[0]), np.deg2rad(float(values[1])), np.deg2rad(float(values[2]))])
    return np.asarray(spherical, dtype=float)


def _ligand_radius_summary() -> str:
    try:
        positions = _ligand_positions()
    except Exception:
        if coordination_geometry in {"custom_xyz", "custom_spherical"}:
            return "custom metal-ligand distances"
        return f"metal-ligand={ligand_radius:g} A"
    distances = positions[:, 0]
    if np.max(distances) - np.min(distances) < 1e-9:
        return f"metal-ligand={float(distances[0]):g} A"
    return (
        f"metal-ligand={float(np.mean(distances)):.4g} A avg "
        f"({float(np.min(distances)):.4g}-{float(np.max(distances)):.4g})"
    )


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


def _broadening_metadata(energy, lorentzian_profile, gaussian_sigma_profile) -> dict[str, object]:
    return {
        "lorentzian_hwhm": _width_summary(lorentzian_profile),
        "gaussian_sigma": _width_summary(gaussian_sigma_profile),
        "energy_points": int(len(energy)),
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
    if np.any(eta <= 0.0):
        raise ValueError("broadening must be positive for the exact spectrum")
    eta_column = eta[:, None] if eta.ndim else float(eta)
    for start in range(0, len(final_energies), chunk_size):
        stop = min(start + chunk_size, len(final_energies))
        delta_e = transition_energy[:, None] + initial_energy - final_energies[None, start:stop]
        spectrum += np.sum(
            strengths[None, start:stop] * eta_column / np.pi / (delta_e**2 + np.asarray(eta_column) ** 2),
            axis=1,
        )
    return spectrum


def _normalize_columns(columns: dict[str, np.ndarray], energy: np.ndarray, mode: str):
    if mode == "none":
        return columns
    reference = columns["total"]
    if mode == "max":
        scale = max(float(np.max(np.abs(reference))), 1e-30)
    elif mode == "area":
        integrate = getattr(np, "trapezoid", np.trapz)
        scale = max(abs(float(integrate(reference, energy))), 1e-30)
    else:
        raise ValueError('normalize must be "max", "area", or "none"')
    return {name: values / scale for name, values in columns.items()}


# ---------------------------------------------------------------------------
# Sizing and validation.


def _estimate_basis_sizes(effective_max_holes: int) -> tuple[int, int]:
    initial_size = 0
    final_size = 0
    core_factor = 6 if _is_core_channel() else 1
    for holes in range(effective_max_holes + 1):
        initial_size += comb(10, n_d_electrons + holes) * comb(10, holes)
        final_d = n_d_electrons + _final_d_offset() + holes
        if 0 <= final_d <= 10:
            final_size += core_factor * comb(10, final_d) * comb(10, holes)
    return initial_size, final_size


def _effective_max_ligand_holes() -> int:
    if max_ligand_holes < 0:
        raise ValueError("max_ligand_holes must be 0 or larger")
    effective = -1
    for holes in range(max_ligand_holes + 1):
        initial_d = n_d_electrons + holes
        final_d = n_d_electrons + _final_d_offset() + holes
        if initial_d > 10 or not 0 <= final_d <= 10:
            break
        try:
            _slater_entry(6, initial_d)
            _slater_entry(_final_p_electrons(), final_d)
        except KeyError:
            break
        effective = holes
    if effective < 0:
        _slater_entry(6, n_d_electrons)
        _slater_entry(_final_p_electrons(), n_d_electrons + _final_d_offset())
    return effective


def _ligand_hole_truncation_reason(effective_max_holes: int) -> str:
    if effective_max_holes >= max_ligand_holes:
        return "none"
    next_hole = effective_max_holes + 1
    initial_d = n_d_electrons + next_hole
    final_d = n_d_electrons + _final_d_offset() + next_hole
    if initial_d > 10 or final_d > 10:
        return f"3d shell is full at h={next_hole}"
    return f"no Appendix-A entry for h={next_hole} (3d{initial_d} / 3d{final_d})"


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
    if photoemission_shell not in {"2p", "3d", "ligand"}:
        raise ValueError('photoemission_shell must be "2p", "3d", or "ligand"')
    if photoemission_shell == "3d" and n_d_electrons < 1:
        raise ValueError("valence-band XPS from the 3d shell needs at least one d electron")
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
        raise ValueError("ligand_positions_spherical is required when coordination_geometry='custom_spherical'")
    if coordination_geometry not in {"custom_xyz", "custom_spherical"} and ligand_radius <= 0.0:
        raise ValueError("ligand_radius must be positive")
    if hybridization_mode not in {"geometry", "symmetry", "scalar"}:
        raise ValueError('hybridization_mode must be "geometry", "symmetry", or "scalar"')
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
    for name, value in (
        ("lorentzian_hwhm", lorentzian_hwhm),
        ("gaussian_sigma", gaussian_sigma),
        ("gaussian_hwhm", gaussian_hwhm),
    ):
        if value < 0.0:
            raise ValueError(f"{name} must be non-negative")
    _slater_entry(6, n_d_electrons)
    _slater_entry(_final_p_electrons(), n_d_electrons + _final_d_offset())


# ---------------------------------------------------------------------------
# Reporting.


def _configuration_energy_rows(max_holes: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for holes in range(max_holes + 1):
        rows.append(
            {
                "label": f"#i{holes + 1}",
                "state": "initial",
                "configuration": _configuration_label(6, n_d_electrons + holes, holes),
                "energy_eV": sector_energy(
                    ligand_holes=holes,
                    core_holes=0,
                    delta=delta,
                    u_charge_transfer=u_charge_transfer,
                    core_hole_potential=core_hole_potential,
                    d_electron_offset=0,
                ),
            }
        )
    for holes in range(max_holes + 1):
        rows.append(
            {
                "label": f"#f{holes + 1}",
                "state": "final",
                "configuration": _configuration_label(
                    _final_p_electrons(),
                    n_d_electrons + _final_d_offset() + holes,
                    holes,
                ),
                "energy_eV": sector_energy(
                    ligand_holes=holes,
                    core_holes=_final_core_holes(),
                    delta=delta,
                    u_charge_transfer=u_charge_transfer,
                    core_hole_potential=core_hole_potential,
                    d_electron_offset=_final_d_offset(),
                ),
            }
        )
    return rows


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
        "zeta_d": entry.zeta_d,
        "fdd2": entry.fdd2,
        "fdd4": entry.fdd4,
        "zeta_2p": entry.zeta_2p,
        "fpd2": entry.fpd2,
        "gpd1": entry.gpd1,
        "gpd3": entry.gpd3,
    }


def _analyze_states(basis, energies, vectors, first_d: int, p_electrons: int) -> list[dict[str, object]]:
    """Decompose each eigenstate into its charge-transfer configurations."""
    hole_counts = np.array([p_ct_ligand_hole_count(state) for state in basis.states])
    rows: list[dict[str, object]] = []
    count = min(len(energies), max(1, n_analyzed_states))
    for idx in range(count):
        vector = vectors[:, idx]
        weights: dict[str, float] = {}
        for holes in np.unique(hole_counts):
            mask = hole_counts == holes
            weight = float(np.sum(np.abs(vector[mask]) ** 2))
            label = _configuration_label(p_electrons, first_d + int(holes), int(holes))
            weights[label] = weight
        rows.append(
            {
                "index": idx,
                "energy_eV": float(energies[idx]),
                "energy_relative_meV": float((energies[idx] - energies[0]) * 1000.0),
                "weights": weights,
                "composition": ", ".join(
                    f"{label} {100.0 * weight:.1f}%" for label, weight in weights.items() if weight > 5e-4
                ),
            }
        )
    return rows


def _parameter_dict(metadata: dict[str, object]) -> dict[str, object]:
    values = {key: globals()[key] for key in sorted(_INPUT_KEYS)}
    values.update({key: globals()[key] for key in _OPTIONAL_INPUT_DEFAULTS})
    values["output_dir"] = str(output_dir)
    values["input_file"] = str(INPUT_FILE)
    values["results"] = metadata
    return values


def _save_parameters(path: Path, metadata: dict[str, object]) -> None:
    path.write_text(
        json.dumps(_parameter_dict(metadata), indent=2, sort_keys=True, default=str),
        encoding="utf-8",
    )


def _save_configuration_energies(path: Path, rows: list[dict[str, object]]) -> None:
    lines = ["label state configuration energy_eV"]
    for row in rows:
        lines.append(f"{row['label']} {row['state']} \"{row['configuration']}\" {row['energy_eV']:.12g}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _save_state_analysis(path: Path, analysis: dict[str, list[dict[str, object]]]) -> None:
    lines = [
        "Python charge-transfer cluster XPS state analysis",
        f"  channel: {photoemission_shell} ({_channel_description()})",
        "  Energies are absolute eigenvalues in eV; dE is relative to the lowest state.",
        "  Weights are configuration projections of each eigenvector.",
        "",
    ]
    for kind, tag in (("initial", "#i"), ("final", "#f")):
        lines.append(f" ====== {kind} sector ({tag}) =========")
        rows = analysis.get(kind, [])
        if not rows:
            lines.append("no states")
            lines.append("")
            continue
        for row in rows:
            lines.append(
                f"state {row['index']:>3}  E = {row['energy_eV']:.9f} eV"
                f"  dE = {row['energy_relative_meV']:.3f} meV"
            )
            for label, weight in row["weights"].items():
                lines.append(f"    {label:<18} {100.0 * weight:8.4f} %")
        lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _plot(energy: np.ndarray, curves: dict[str, np.ndarray], out_png: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; text output was still saved.")
        return

    plot_energy = energy + plot_absolute_energy_offset if plot_use_absolute_energy else energy
    mask = (energy >= plot_relative_energy_min) & (energy <= plot_relative_energy_max)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(plot_energy[mask], curves["total"][mask], color="black", lw=1.4, label="total")
    if spin_resolved:
        ax.plot(plot_energy[mask], curves["spin_up"][mask], color="crimson", lw=1.0, ls="--", label="spin up")
        ax.plot(plot_energy[mask], curves["spin_down"][mask], color="royalblue", lw=1.0, ls="--", label="spin down")
        ax.legend(frameon=False, fontsize=9)

    ax.set_xlabel("Binding energy (eV)" if plot_use_absolute_energy else "Relative binding energy (eV)")
    ax.set_ylabel("Intensity (arb. units)")
    ax.set_title(f"{element} d{n_d_electrons} {photoemission_shell} XPS")
    if plot_binding_energy_axis:
        ax.invert_xaxis()
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()

    if save_png:
        fig.savefig(out_png, dpi=200)
        print("saved:", out_png)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)


def _scale_tag() -> str:
    scales = (fdd2_scale, fdd4_scale, fpd2_scale, gpd1_scale, gpd3_scale)
    if all(abs(scale - scales[0]) < 1e-12 for scale in scales):
        return f"scale{scales[0]:.3f}"
    return "separate_scales"


if __name__ == "__main__":
    main()

"""Input file for `run_xps.py` -- Ba2FeSi2O7 Fe 2p XPS.

Same cluster as `Fe_Ba2FeSi2O7.py` (the XAS/LD case): an Fe2+ d6 ion in a
flattened FeO4 tetrahedron. Only the spectroscopy differs, so the physical
parameters are deliberately kept identical -- fitting XAS and XPS with one
parameter set is the whole point of the charge-transfer cluster model, and
2p XPS satellites are far more sensitive to Delta than the XAS line shape is.

Fe d6 sectors for XPS:
    initial: 2p6 3d6, 2p6 3d7 L, 2p6 3d8 L2
    final:   2p5 3d6, 2p5 3d7 L, 2p5 3d8 L2      <- d count unchanged

Compare with the XAS final sectors (2p5 3d7, 2p5 3d8 L, 2p5 3d9 L2), which
carry one more d electron because the core electron is promoted rather than
emitted.
"""

# ---------------------------------------------------------------------------
# Cluster and electron configuration.

case_name = "Fe_Ba2FeSi2O7_2pXPS"
ion = "Fe2+"
element = "auto"
n_d_electrons = "auto"
max_ligand_holes = 2


# ---------------------------------------------------------------------------
# Photoemission channel.

photoemission_shell = "2p"  # "2p" (core), "3d" (valence), or "ligand"
spin_resolved = True


# ---------------------------------------------------------------------------
# XTLS-like cluster parameters. Kept identical to the XAS input.

ten_dq = 0.30
ligand_ten_dq = 0.00  # set to 2*Tpp to split the ligand molecular orbitals
delta = 5.5
u_charge_transfer = 4.25
core_hole_potential = 5.75

pd_sigma = -1.5
pd_ratio = -2.0
d_ref = 2.0


# ---------------------------------------------------------------------------
# Ligand geometry.

coordination_geometry = "tetrahedral"
ligand_radius = 1.99626
ligand_angle_offset_deg = 7.8353
r2 = 0.402
r4 = 0.380

ligand_positions_xyz = None
ligand_positions_spherical = None


# ---------------------------------------------------------------------------
# Hybridization.

hybridization_mode = "geometry"  # "symmetry", "geometry", or "scalar"
v_eg = 2.20  # used only when hybridization_mode = "symmetry"
v_t2g = 1.10
hopping = 0.8  # used only when hybridization_mode = "scalar"


# ---------------------------------------------------------------------------
# Slater integral and spin-orbit reduction factors.

fdd2_scale = 0.825
fdd4_scale = 0.825
fpd2_scale = 0.825
gpd1_scale = 0.825
gpd3_scale = 0.825
so3d_scale = 1.0
so2p_scale = 1.0


# ---------------------------------------------------------------------------
# Calculation settings.

n_initial_states = 10
temperature_kelvin = 370.0
spectrum_method = "lanczos"  # "lanczos" or "exact"
n_recursion = 500
n_analyzed_states = 10

energy_min = -5.0
energy_max = 35.0
energy_step = 0.01
energy_shift = 0.0
normalize = "max"  # "max", "area", or "none"


# ---------------------------------------------------------------------------
# Broadening.

lorentzian_hwhm = 0.35
gaussian_sigma = 0.00
gaussian_hwhm = 0.60
lorentzian_hwhm_points = []
gaussian_hwhm_points = []


# ---------------------------------------------------------------------------
# Output.

save_txt = True
save_png = True
save_parameters = True
save_configuration_energies = True
save_state_analysis = True
estimate_only = False
output_dir = "auto"


# ---------------------------------------------------------------------------
# Plot settings.

show_plot = True
plot_binding_energy_axis = True
plot_use_absolute_energy = True
plot_absolute_energy_offset = "auto"  # ~711 eV, the Fe 2p3/2 binding energy
plot_relative_energy_min = -3.0
plot_relative_energy_max = 32.0

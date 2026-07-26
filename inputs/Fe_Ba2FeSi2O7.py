"""Input file for `run_xas.py`.

Edit this file in Spyder to define the cluster, calculation, and plot settings.
Then run `run_xas.py`.
"""

# ---------------------------------------------------------------------------
# Cluster and electron configuration.

case_name = "Fe_Ba2FeSi2O7"
ion = "Fe2+"
element = "auto"
n_d_electrons = "auto"
max_ligand_holes = 2

# To start a different case, copy this input file and change for example:
# case_name = "Co_Ba2CoO4"; ion = "Co4+"
# case_name = "Ni_NiO"; ion = "Ni2+"
# case_name = "Mn_example"; ion = "Mn2+"
# The runner then fills element, nominal d count, output_dir, and L3 plot offset.

# Fe d6 charge-transfer sectors:
# initial: 2p6 3d6, 2p6 3d7 L, 2p6 3d8 L2
# final:   2p5 3d7, 2p5 3d8 L, 2p5 3d9 L2


# ---------------------------------------------------------------------------
# XTLS-like cluster parameters.

ten_dq = 0.30
delta = 5.5
pd_sigma = -1.5
pd_ratio = -2.0
d_ref = 2.0
u_charge_transfer = 4.25
core_hole_potential = 5.75


# ---------------------------------------------------------------------------
# Ligand geometry.

coordination_geometry = "tetrahedral"
# Available presets:
# "tetrahedral", "octahedral", "square_planar", "square_pyramidal",
# "trigonal_bipyramidal", "custom_xyz", "custom_spherical"
ligand_radius = 1.99626
ligand_angle_offset_deg = 7.8353
r2 = 0.402
r4 = 0.380

# For coordination_geometry = "custom_xyz", use Angstrom Cartesian coordinates
# relative to the metal atom at (0, 0, 0). The ligand label is optional.
ligand_positions_xyz = None
# ligand_positions_xyz = [
#     ("O1",  1.20,  1.20,  1.10),
#     ("O2", -1.20, -1.20,  1.05),
#     ("O3", -1.15,  1.25, -1.00),
#     ("O4",  1.18, -1.22, -1.08),
# ]

# For coordination_geometry = "custom_spherical", use
# (label, radius_A, theta_deg, phi_deg). The ligand label is optional.
ligand_positions_spherical = None
# ligand_positions_spherical = [
#     ("O1", 1.99,  62.5,  -45.0),
#     ("O2", 1.98,  62.0,  135.0),
#     ("O3", 2.01, 118.5,   45.0),
#     ("O4", 1.97, 119.0, -135.0),
# ]


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
# Hybridization.

hybridization_mode = "geometry"  # "geometry" or "scalar"
hopping = 0.8  # used only when hybridization_mode = "scalar"


# ---------------------------------------------------------------------------
# Calculation settings.

n_initial_states = 10
temperature_kelvin = 370.0
spectrum_method = "lanczos"  # "lanczos" or "exact"
n_recursion = 500
energy_min = -20.0
energy_max = 20.0
energy_step = 0.01
energy_shift = 0.0
normalize = "max"  # "max", "area", or "none"


# ---------------------------------------------------------------------------
# Broadening.

lorentzian_hwhm = 0.20
gaussian_sigma = 0.00
gaussian_hwhm = 0.22
lorentzian_hwhm_points = [(-3.0, 0.20), (10.0, 0.35)]
gaussian_hwhm_points = []


# ---------------------------------------------------------------------------
# Experimental LD geometry.

make_experimental_geometry_curves = True
grazing_angle_deg = 23.5
inplane_curve = "ab"  # "ab", "x", or "y"


# ---------------------------------------------------------------------------
# Output and state analysis.

save_txt = True
save_png = True
save_parameters = True
save_configuration_energies = True
save_state_analysis = True
estimate_only = False
n_analyzed_states = 10
output_dir = "auto"


# ---------------------------------------------------------------------------
# Plot settings.

show_plot = True
show_cluster_inset = False
plot_ld_stacked = True
plot_iso_only = False
plot_xas_offset = 0.70
plot_ld_offset = 0.00
plot_ld_scale = 0.60
plot_use_absolute_energy = True
plot_absolute_energy_offset = "auto"
plot_relative_energy_min = -5.0
plot_relative_energy_max = 25.0


# ---------------------------------------------------------------------------
# Optional XTLS overlay.

overlay_xtls = True
xtls_path = ROOT / "data" / "Fe_L_spectrum.txt"
xtls_energy_shift = 0.0
xtls_scale = 1.0
xtls_iso_column = "iso_broadened"

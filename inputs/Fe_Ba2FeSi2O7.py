"""Ba2FeSi2O7 -- Fe2+ (3d6) in a flattened FeO4 tetrahedron.

One file, both spectroscopies. Run it with:

    python run.py inputs/Fe_Ba2FeSi2O7.py xas
    python run.py inputs/Fe_Ba2FeSi2O7.py xps
    python run.py inputs/Fe_Ba2FeSi2O7.py both

Naming rule: settings written without a prefix are shared. Where the two
calculations need different values, write `xas_<name>` and `xps_<name>` -- a
prefixed setting overrides the bare one for that spectroscopy. Settings only
one of them knows about (LD geometry, photoemission channel) need no prefix,
since the other simply ignores them.

The cluster block below is therefore the single source of truth: change Delta
once and both spectra move together, which is the whole reason to fit them
jointly.

Charge-transfer sectors for Fe d6:
    initial:   2p6 3d6,  2p6 3d7 L,  2p6 3d8 L2
    XAS final: 2p5 3d7,  2p5 3d8 L,  2p5 3d9 L2    (core electron promoted)
    XPS final: 2p5 3d6,  2p5 3d7 L,  2p5 3d8 L2    (photoelectron emitted)
"""

# ===========================================================================
# SHARED -- the cluster itself
# ===========================================================================

case_name = "Fe_Ba2FeSi2O7"
ion = "Fe2+"
element = "auto"
n_d_electrons = "auto"
max_ligand_holes = 2

# To start a different case, copy this file and change for example:
# case_name = "Co_Ba2CoO4"; ion = "Co4+"
# case_name = "Mn_example";  ion = "Mn2+"
# The runner then fills element, nominal d count, output_dir and plot offset.


# --- Charge-transfer and crystal-field parameters --------------------------

ten_dq = 0.30
delta = 5.5
u_charge_transfer = 4.25
core_hole_potential = 5.75
ligand_ten_dq = 0.00  # ligand molecular-orbital splitting, XTLS `10Dq(Ld) = 2*Tpp`


# --- Ligand geometry -------------------------------------------------------

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


# --- Hybridization ---------------------------------------------------------

hybridization_mode = "geometry"  # "geometry", "symmetry", or "scalar"
pd_sigma = -1.5  # "geometry" mode
pd_ratio = -2.0
d_ref = 2.0
v_eg = 2.20  # "symmetry" mode, XTLS `VOh(Ld 3d) = {V(eg), V(t2g)}`
v_t2g = 1.10
hopping = 0.8  # "scalar" mode


# --- Slater integral and spin-orbit reduction factors ----------------------

fdd2_scale = 0.825
fdd4_scale = 0.825
fpd2_scale = 0.825
gpd1_scale = 0.825
gpd3_scale = 0.825
so3d_scale = 1.0
so2p_scale = 1.0


# --- Solver ----------------------------------------------------------------

n_initial_states = 10
temperature_kelvin = 370.0
spectrum_method = "lanczos"  # "lanczos" or "exact"
n_recursion = 500
n_analyzed_states = 10
energy_step = 0.01
energy_shift = 0.0
normalize = "max"  # "max", "area", or "none"


# --- Output ----------------------------------------------------------------

save_txt = True
save_png = True
save_parameters = True
save_configuration_energies = True
save_state_analysis = True
estimate_only = False
output_dir = "auto"

show_plot = True
plot_use_absolute_energy = True
plot_absolute_energy_offset = "auto"  # L3 edge for XAS, 2p3/2 binding energy for XPS


# ===========================================================================
# XAS -- absorption and linear dichroism
# ===========================================================================

xas_energy_min = -20.0
xas_energy_max = 20.0
xas_lorentzian_hwhm = 0.20
xas_gaussian_sigma = 0.00
xas_gaussian_hwhm = 0.22
xas_lorentzian_hwhm_points = [(-3.0, 0.20), (10.0, 0.35)]
xas_gaussian_hwhm_points = []
xas_plot_relative_energy_min = -5.0
xas_plot_relative_energy_max = 25.0

# XPS never reads these, so they need no prefix.
make_experimental_geometry_curves = True
grazing_angle_deg = 23.5
inplane_curve = "ab"  # "ab", "x", or "y"

show_cluster_inset = False
plot_ld_stacked = True
plot_iso_only = False
plot_xas_offset = 0.70
plot_ld_offset = 0.00
plot_ld_scale = 0.60

overlay_xtls = True
xtls_path = ROOT / "data" / "Fe_L_spectrum.txt"
xtls_energy_shift = 0.0
xtls_scale = 1.0
xtls_iso_column = "iso_broadened"


# ===========================================================================
# XPS -- photoemission
# ===========================================================================

xps_energy_min = -5.0
xps_energy_max = 35.0
xps_lorentzian_hwhm = 0.35
xps_gaussian_sigma = 0.00
xps_gaussian_hwhm = 0.60
xps_lorentzian_hwhm_points = []
xps_gaussian_hwhm_points = []
xps_plot_relative_energy_min = -3.0
xps_plot_relative_energy_max = 32.0

# XAS never reads these, so they need no prefix.
photoemission_shell = "2p"  # "2p" (core level), "3d" (valence band), or "ligand"
spin_resolved = True
plot_binding_energy_axis = True  # XPS convention: binding energy increases to the left

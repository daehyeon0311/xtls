"""NiO -- Ni2+ (3d8) in a NiO6 octahedron.

The XPS block reproduces the worked example in the XTLS 9.0 manual
(`Xtls900_man_Jpn.pdf`, pages 20-21) and is the reference case for validating
this implementation, since the manual prints every parameter it uses.

    XTLS X-card                       this file
    -------------------------------   -----------------------------------
    U3d3d = 7.3                       u_charge_transfer = 7.3
    U3d2p = 8.5                       core_hole_potential = 8.5
    Dlt   = 4.7                       delta = 4.7
    VEg   = 2.2, V(t2g) = VEg/2       v_eg = 2.2, v_t2g = 1.1
    TenDq = 0.7                       ten_dq = 0.7
    Tpp   = 0.7, 10Dq(Ld) = 2*Tpp     ligand_ten_dq = 1.4
    Mode  = XPS                       photoemission_shell = "2p"

The manual sets the Slater integrals directly to their Hartree-Fock values
rather than applying the usual `Red` reduction, so all scale factors are 1.0.
Those values agree with the Appendix-A table shipped here to three decimals.

CAVEAT: the manual gives no XAS card for NiO, so `xas` runs off the same
cluster but has not been checked against anything. Treat XPS as the validated
mode here.

Run with:
    python run.py inputs/NiO.py xps
"""

# ===========================================================================
# SHARED -- the cluster itself
# ===========================================================================

case_name = "NiO"
ion = "Ni2+"
element = "auto"
n_d_electrons = "auto"
max_ligand_holes = 2

# Ni d8 sectors, matching the manual's CNFG: block
#     initial:   2p6 3d8, 2p6 3d9 L, 2p6 3d10 L2   (#i1, #i2, #i3)
#     XPS final: 2p5 3d8, 2p5 3d9 L, 2p5 3d10 L2   (#f1, #f2, #f3)


# --- Charge-transfer and crystal-field parameters --------------------------

ten_dq = 0.70
delta = 4.70
u_charge_transfer = 7.30
core_hole_potential = 8.50
ligand_ten_dq = 1.40  # XTLS `10Dq(Ld) = 2*Tpp` with Tpp = 0.7


# --- Ligand geometry -------------------------------------------------------

coordination_geometry = "octahedral"
ligand_radius = 2.084  # rock-salt Ni-O distance
ligand_angle_offset_deg = 0.0
r2 = 0.402
r4 = 0.380

ligand_positions_xyz = None
ligand_positions_spherical = None


# --- Hybridization ---------------------------------------------------------
#
# "symmetry" is the XTLS `VOh` form and is how published cluster parameters
# are almost always quoted.

hybridization_mode = "symmetry"  # "geometry", "symmetry", or "scalar"
v_eg = 2.20
v_t2g = 1.10
pd_sigma = -1.5  # "geometry" mode
pd_ratio = -2.0
d_ref = 2.0
hopping = 0.8  # "scalar" mode


# --- Slater integral and spin-orbit reduction factors ----------------------
#
# The manual quotes bare Hartree-Fock values, so these stay at 1.0.
# XTLS otherwise defaults to Red = 0.8 for the Rk parameters.

fdd2_scale = 1.0
fdd4_scale = 1.0
fpd2_scale = 1.0
gpd1_scale = 1.0
gpd3_scale = 1.0
so3d_scale = 1.0
so2p_scale = 1.0


# --- Solver ----------------------------------------------------------------
#
# The Ni2+ 3A2g ground state is a spin triplet, so three initial states are
# averaged. XTLS calls this `Ninit`.

n_initial_states = 3
temperature_kelvin = 300.0
spectrum_method = "lanczos"  # "lanczos" or "exact"
n_recursion = 300
n_analyzed_states = 6
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
plot_absolute_energy_offset = "auto"


# ===========================================================================
# XAS -- absorption and linear dichroism   (not validated for this material)
# ===========================================================================

xas_energy_min = -15.0
xas_energy_max = 20.0
xas_lorentzian_hwhm = 0.25
xas_gaussian_sigma = 0.00
xas_gaussian_hwhm = 0.30
xas_lorentzian_hwhm_points = []
xas_gaussian_hwhm_points = []
xas_plot_relative_energy_min = -5.0
xas_plot_relative_energy_max = 20.0

# Octahedral NiO is cubic, so there is no linear dichroism to speak of.
make_experimental_geometry_curves = False
grazing_angle_deg = 23.5
inplane_curve = "ab"

show_cluster_inset = False
plot_ld_stacked = False
plot_iso_only = True
plot_xas_offset = 0.70
plot_ld_offset = 0.00
plot_ld_scale = 0.60

overlay_xtls = False
xtls_path = None
xtls_energy_shift = 0.0
xtls_scale = 1.0
xtls_iso_column = "iso_broadened"


# ===========================================================================
# XPS -- photoemission   (validated against the manual)
# ===========================================================================

xps_energy_min = -3.0
xps_energy_max = 30.0
xps_lorentzian_hwhm = 0.40
xps_gaussian_sigma = 0.00
xps_gaussian_hwhm = 0.80
xps_lorentzian_hwhm_points = []
xps_gaussian_hwhm_points = []
xps_plot_relative_energy_min = -2.0
xps_plot_relative_energy_max = 28.0

photoemission_shell = "2p"  # "2p" (core level), "3d" (valence band), or "ligand"
spin_resolved = False
plot_binding_energy_axis = True

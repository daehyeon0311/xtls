"""Input file for `run_xps.py` -- NiO Ni 2p XPS.

This reproduces the worked example in section 5 of the XTLS 9.0 manual
(`Xtls900_man_Jpn.pdf`, pages 20-21), which uses a NiO6 octahedral cluster.
It doubles as the reference case for the Python implementation, because the
manual prints every parameter it uses.

XTLS X-card                       this file
-------------------------------   ---------------------------------------
U3d3d = 7.3                       u_charge_transfer = 7.3
U3d2p = 8.5                       core_hole_potential = 8.5
Dlt   = 4.7                       delta = 4.7
VEg   = 2.2, V(t2g) = VEg/2       v_eg = 2.2, v_t2g = 1.1
TenDq = 0.7                       ten_dq = 0.7
Tpp   = 0.7, 10Dq(Ld) = 2*Tpp     ligand_ten_dq = 1.4
Mode  = XPS                       photoemission_shell = "2p"

The manual sets the Slater integrals directly to their Hartree-Fock values
(`Rk(#i1 3d 3d) = {12.234, 7.598}` and so on) rather than applying the usual
`Red` reduction, so all scale factors here are 1.0. Those values agree with
the Appendix-A table shipped with this package to three decimals.
"""

# ---------------------------------------------------------------------------
# Cluster and electron configuration.

case_name = "NiO_2pXPS"
ion = "Ni2+"
element = "auto"
n_d_electrons = "auto"
max_ligand_holes = 2

# Ni d8 charge-transfer sectors, matching the manual's CNFG: block:
# initial: 2p6 3d8, 2p6 3d9 L, 2p6 3d10 L2      (#i1, #i2, #i3)
# final:   2p5 3d8, 2p5 3d9 L, 2p5 3d10 L2      (#f1, #f2, #f3)


# ---------------------------------------------------------------------------
# Photoemission channel.

photoemission_shell = "2p"  # "2p" (core), "3d" (valence), or "ligand"
spin_resolved = False


# ---------------------------------------------------------------------------
# XTLS-like cluster parameters.

ten_dq = 0.70
ligand_ten_dq = 1.40  # XTLS `10Dq(Ld) = 2*Tpp`; splits the ligand molecular orbitals
delta = 4.70
u_charge_transfer = 7.30
core_hole_potential = 8.50

# Used only when hybridization_mode = "geometry".
pd_sigma = -1.5
pd_ratio = -2.0
d_ref = 2.0


# ---------------------------------------------------------------------------
# Ligand geometry.

coordination_geometry = "octahedral"
ligand_radius = 2.084  # NiO rock-salt Ni-O distance
ligand_angle_offset_deg = 0.0
r2 = 0.402
r4 = 0.380

ligand_positions_xyz = None
ligand_positions_spherical = None


# ---------------------------------------------------------------------------
# Hybridization.
#
# "symmetry" is the XTLS `VOh(#sc1 #sc2 Ld 3d) = {V(eg), V(t2g)}` form and is
# what published cluster parameters almost always quote.
# "geometry" instead builds the d-L hopping from the ligand positions using
# Harrison scaling; "scalar" gives every d orbital the same hopping.

hybridization_mode = "symmetry"  # "symmetry", "geometry", or "scalar"
v_eg = 2.20
v_t2g = 1.10
hopping = 0.8  # used only when hybridization_mode = "scalar"


# ---------------------------------------------------------------------------
# Slater integral and spin-orbit reduction factors.
#
# The manual's NiO example quotes bare Hartree-Fock values, so these stay at
# 1.0. XTLS otherwise defaults to Red = 0.8 for the Rk parameters.

fdd2_scale = 1.0
fdd4_scale = 1.0
fpd2_scale = 1.0
gpd1_scale = 1.0
gpd3_scale = 1.0
so3d_scale = 1.0
so2p_scale = 1.0


# ---------------------------------------------------------------------------
# Calculation settings.
#
# The Ni2+ 3A2g ground state is a spin triplet, so three initial states are
# averaged. XTLS calls this `Ninit`.

n_initial_states = 3
temperature_kelvin = 300.0
spectrum_method = "lanczos"  # "lanczos" or "exact"
n_recursion = 300
n_analyzed_states = 6

energy_min = -3.0
energy_max = 30.0
energy_step = 0.01
energy_shift = 0.0
normalize = "max"  # "max", "area", or "none"


# ---------------------------------------------------------------------------
# Broadening.

lorentzian_hwhm = 0.40
gaussian_sigma = 0.00
gaussian_hwhm = 0.80
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
plot_binding_energy_axis = True  # XPS convention: binding energy increases to the left
plot_use_absolute_energy = True
plot_absolute_energy_offset = "auto"  # ~854 eV, the Ni 2p3/2 binding energy
plot_relative_energy_min = -2.0
plot_relative_energy_max = 28.0

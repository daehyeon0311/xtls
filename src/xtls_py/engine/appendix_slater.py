from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SlaterEntry:
    """Appendix-A 3d/4d Slater data from Haverkort's thesis.

    Values are in eV. `r2` and `r4` are included for traceability but are not
    used directly by the current exact-diagonalization engine.
    """

    element: str
    p_electrons: int
    d_electrons: int
    r2: float | None
    r4: float | None
    zeta_d: float
    fdd2: float
    fdd4: float
    zeta_2p: float | None = None
    fpd2: float | None = None
    gpd1: float | None = None
    gpd3: float | None = None

# 3d table: K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
# 4d table: Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd

APPENDIX_A_3D_ROWS: tuple[tuple[str, int, int, float | None, ...], ...] = (
    # Rows from Appendix A, "Slater integrals for 3d elements".
    # `None` for r2/r4 corresponds to the starred entries in the appendix.
    ("K", 6, 1, None, None, 0.000, 0.000, 0.000, None, None, None, None),
    ("K", 5, 2, None, None, 0.002, 2.723, 1.622, 1.880, 1.045, 0.583, 0.329),
    ("K", 5, 1, 1.687, 6.860, 0.005, 0.000, 0.000, 1.880, 2.202, 1.318, 0.745),
    ("Ca", 6, 2, None, None, 0.001, 2.433, 1.452, None, None, None, None),
    ("Ca", 5, 3, None, None, 0.005, 4.205, 2.529, 2.405, 1.970, 1.204, 0.681),
    ("Ca", 6, 1, None, None, 0.004, 0.000, 0.000, None, None, None, None),
    ("Ca", 5, 2, 1.202, 3.577, 0.008, 6.258, 3.851, 2.405, 2.917, 1.851, 1.048),
    ("Ca", 5, 1, 0.781, 1.306, 0.011, 0.000, 0.000, 2.404, 3.793, 2.511, 1.423),
    ("Sc", 6, 3, None, None, 0.005, 3.855, 2.311, None, None, None, None),
    ("Sc", 5, 4, None, None, 0.010, 5.326, 3.224, 3.034, 2.696, 1.733, 0.982),
    ("Sc", 6, 2, 1.356, 4.559, 0.007, 5.851, 3.590, None, None, None, None),
    ("Sc", 5, 3, 0.943, 2.247, 0.013, 7.116, 4.385, 3.033, 3.532, 2.331, 1.321),
    ("Sc", 6, 1, 0.844, 1.535, 0.010, 0.000, 0.000, None, None, None, None),
    ("Sc", 5, 2, 0.658, 0.946, 0.017, 8.530, 5.321, 3.032, 4.332, 2.950, 1.674),
    ("Sc", 5, 1, 0.516, 0.535, 0.020, 0.000, 0.000, 3.032, 5.095, 3.582, 2.035),
    ("Ti", 6, 4, None, None, 0.009, 5.002, 3.019, None, None, None, None),
    ("Ti", 5, 5, 1.326, 6.323, 0.015, 6.261, 3.806, 3.778, 3.329, 2.216, 1.257),
    ("Ti", 6, 3, 1.029, 2.692, 0.012, 6.780, 4.167, None, None, None, None),
    ("Ti", 5, 4, 0.775, 1.546, 0.019, 7.888, 4.865, 3.777, 4.098, 2.783, 1.578),
    ("Ti", 6, 2, 0.699, 1.075, 0.016, 8.243, 5.132, None, None, None, None),
    ("Ti", 5, 3, 0.566, 0.713, 0.023, 9.213, 5.744, 3.776, 4.849, 3.376, 1.917),
    ("Ti", 6, 1, 0.538, 0.587, 0.019, 0.000, 0.000, None, None, None, None),
    ("Ti", 5, 2, 0.453, 0.421, 0.027, 10.342, 6.499, 3.776, 5.580, 3.989, 2.267),
    ("Ti", 5, 1, 0.381, 0.281, 0.032, 0.000, 0.000, 3.776, 6.301, 4.626, 2.632),
    ("V", 6, 5, 1.447, 7.507, 0.014, 5.970, 3.620, None, None, None, None),
    ("V", 5, 6, 1.044, 3.981, 0.022, 7.095, 4.324, 4.652, 3.911, 2.673, 1.517),
    ("V", 6, 4, 0.829, 1.786, 0.018, 7.599, 4.676, None, None, None, None),
    ("V", 5, 5, 0.655, 1.124, 0.026, 8.612, 5.314, 4.651, 4.634, 3.218, 1.827),
    ("V", 6, 3, 0.595, 0.793, 0.022, 8.961, 5.576, None, None, None, None),
    ("V", 5, 4, 0.494, 0.553, 0.031, 9.875, 6.152, 4.650, 5.351, 3.792, 2.154),
    ("V", 6, 2, 0.470, 0.456, 0.027, 10.126, 6.353, None, None, None, None),
    ("V", 5, 3, 0.403, 0.338, 0.036, 10.973, 6.887, 4.649, 6.056, 4.389, 2.495),
    ("V", 6, 1, 0.392, 0.299, 0.031, 0.000, 0.000, None, None, None, None),
    ("V", 5, 2, 0.342, 0.231, 0.041, 11.963, 7.554, 4.650, 6.757, 5.011, 2.852),
    ("Cr", 6, 6, 1.121, 4.612, 0.021, 6.833, 4.155, None, None, None, None),
    ("Cr", 5, 7, 0.856, 2.709, 0.030, 7.866, 4.803, 5.671, 4.460, 3.112, 1.767),
    ("Cr", 6, 5, 0.692, 1.267, 0.025, 8.355, 5.146, None, None, None, None),
    ("Cr", 5, 6, 0.565, 0.851, 0.035, 9.303, 5.742, 5.669, 5.151, 3.641, 2.068),
    ("Cr", 6, 4, 0.515, 0.605, 0.030, 9.648, 6.001, None, None, None, None),
    ("Cr", 5, 5, 0.437, 0.439, 0.041, 10.521, 6.551, 5.668, 5.840, 4.201, 2.387),
    ("Cr", 6, 3, 0.416, 0.362, 0.035, 10.776, 6.754, None, None, None, None),
    ("Cr", 5, 4, 0.361, 0.276, 0.047, 11.594, 7.270, 5.667, 6.524, 4.785, 2.721),
    ("Cr", 6, 2, 0.351, 0.243, 0.041, 11.793, 7.437, None, None, None, None),
    ("Cr", 5, 3, 0.310, 0.191, 0.053, 12.572, 7.928, 5.668, 7.209, 5.394, 3.070),
    ("Mn", 6, 7, None, None, 0.032, 5.312, 3.268, None, None, None, None),
    ("Mn", 5, 8, 0.721, 1.945, 0.040, 8.594, 5.255, 6.849, 4.987, 3.539, 2.010),
    ("Mn", 6, 6, 0.592, 0.942, 0.035, 9.072, 5.590, None, None, None, None),
    ("Mn", 5, 7, 0.495, 0.662, 0.046, 9.971, 6.156, 6.847, 5.652, 4.056, 2.304),
    ("Mn", 6, 5, 0.452, 0.475, 0.040, 10.315, 6.413, None, None, None, None),
    ("Mn", 5, 6, 0.390, 0.355, 0.053, 11.154, 6.942, 6.846, 6.320, 4.603, 2.617),
    ("Mn", 6, 4, 0.371, 0.293, 0.046, 11.414, 7.147, None, None, None, None),
    ("Mn", 5, 5, 0.326, 0.228, 0.059, 12.209, 7.648, 6.845, 6.987, 5.176, 2.944),
    ("Mn", 6, 3, 0.316, 0.201, 0.052, 12.414, 7.819, None, None, None, None),
    ("Mn", 5, 4, 0.282, 0.161, 0.066, 13.176, 8.299, 6.845, 7.657, 5.773, 3.287),
    ("Fe", 6, 8, None, None, 0.021, 2.087, 1.275, None, None, None, None),
    ("Fe", 5, 9, 0.620, 1.452, 0.053, 9.293, 5.688, 8.203, 5.498, 3.957, 2.248),
    ("Fe", 6, 7, 0.515, 0.723, 0.045, 9.761, 6.017, None, None, None, None),
    ("Fe", 5, 8, 0.438, 0.527, 0.059, 10.622, 6.559, 8.202, 6.142, 4.463, 2.537),
    ("Fe", 6, 6, 0.402, 0.380, 0.052, 10.965, 6.815, None, None, None, None),
    ("Fe", 5, 7, 0.351, 0.292, 0.067, 11.778, 7.327, 8.200, 6.792, 5.000, 2.843),
    ("Fe", 6, 5, 0.334, 0.240, 0.059, 12.042, 7.534, None, None, None, None),
    ("Fe", 5, 6, 0.296, 0.191, 0.074, 12.817, 8.023, 8.199, 7.445, 5.563, 3.165),
    ("Fe", 6, 4, 0.287, 0.167, 0.066, 13.029, 8.197, None, None, None, None),
    ("Fe", 5, 5, 0.258, 0.136, 0.082, 13.775, 8.667, 8.199, 8.101, 6.150, 3.502),
    ("Co", 6, 9, None, None, 0.027, 2.196, 1.348, None, None, None, None),
    ("Co", 5, 10, 0.541, 1.117, 0.067, 9.969, 6.108, 9.752, 5.996, 4.367, 2.482),
    ("Co", 6, 8, 0.453, 0.570, 0.059, 10.430, 6.431, None, None, None, None),
    ("Co", 5, 9, 0.392, 0.427, 0.075, 11.261, 6.954, 9.750, 6.624, 4.865, 2.766),
    ("Co", 6, 7, 0.360, 0.310, 0.066, 11.604, 7.209, None, None, None, None),
    ("Co", 5, 8, 0.318, 0.243, 0.083, 12.395, 7.707, 9.748, 7.259, 5.394, 3.068),
    ("Co", 6, 6, 0.302, 0.200, 0.074, 12.662, 7.916, None, None, None, None),
    ("Co", 5, 7, 0.270, 0.161, 0.092, 13.421, 8.394, 9.746, 7.899, 5.947, 3.384),
    ("Co", 6, 5, 0.262, 0.141, 0.082, 13.638, 8.572, None, None, None, None),
    ("Co", 5, 6, 0.237, 0.116, 0.101, 14.372, 9.034, 9.746, 8.544, 6.525, 3.716),
    ("Ni", 6, 9, 0.404, 0.458, 0.074, 11.084, 6.835, None, None, None, None),
    ("Ni", 5, 10, 0.353, 0.351, 0.093, 11.890, 7.343, 11.509, 7.098, 5.262, 2.993),
    ("Ni", 6, 8, 0.325, 0.256, 0.083, 12.233, 7.597, None, None, None, None),
    ("Ni", 5, 9, 0.289, 0.204, 0.102, 13.005, 8.084, 11.507, 7.720, 5.783, 3.290),
    ("Ni", 6, 7, 0.275, 0.168, 0.091, 13.276, 8.294, None, None, None, None),
    ("Ni", 5, 8, 0.248, 0.137, 0.112, 14.021, 8.763, 11.506, 8.349, 6.329, 3.602),
    ("Ni", 6, 6, 0.240, 0.120, 0.101, 14.244, 8.944, None, None, None, None),
    ("Ni", 5, 7, 0.218, 0.100, 0.122, 14.965, 9.399, 11.505, 8.984, 6.898, 3.929),
    ("Cu", 6, 9, 0.295, 0.214, 0.102, 12.854, 7.980, None, None, None, None),
    ("Cu", 5, 10, 0.265, 0.173, 0.124, 13.611, 8.457, 13.498, 8.177, 6.169, 3.510),
    ("Cu", 6, 8, 0.252, 0.143, 0.112, 13.885, 8.669, None, None, None, None),
    ("Cu", 5, 9, 0.228, 0.118, 0.135, 14.617, 9.130, 13.496, 8.796, 6.708, 3.818),
    ("Cu", 6, 7, 0.221, 0.103, 0.123, 14.845, 9.313, None, None, None, None),
    ("Cu", 5, 8, 0.202, 0.087, 0.147, 15.556, 9.762, 13.495, 9.421, 7.270, 4.141),
    ("Zn", 6, 9, 0.232, 0.122, 0.136, 14.489, 9.041, None, None, None, None),
    ("Zn", 5, 10, 0.211, 0.102, 0.162, 15.210, 9.495, 15.738, 9.240, 7.084, 4.033),
    ("Zn", 6, 8, 0.204, 0.089, 0.147, 15.443, 9.681, None, None, None, None),
    ("Zn", 5, 9, 0.187, 0.076, 0.175, 16.145, 10.124, 15.737, 9.857, 7.639, 4.352),
)


APPENDIX_A_3D: dict[tuple[str, int, int], SlaterEntry] = {
    (element, p_electrons, d_electrons): SlaterEntry(
        element,
        p_electrons,
        d_electrons,
        r2,
        r4,
        zeta_d,
        fdd2,
        fdd4,
        zeta_2p,
        fpd2,
        gpd1,
        gpd3,
    )
    for (
        element,
        p_electrons,
        d_electrons,
        r2,
        r4,
        zeta_d,
        fdd2,
        fdd4,
        zeta_2p,
        fpd2,
        gpd1,
        gpd3,
    ) in APPENDIX_A_3D_ROWS
}


APPENDIX_A_4D_ROWS: tuple[tuple[str, int, int, float | None, ...], ...] = (
    # Rows from Appendix A, "Slater integrals for 4d elements".
    ("Rb", 6, 1, None, None, 0.001, 0.000, 0.000, None, None, None, None),
    ("Rb", 5, 2, None, None, 0.009, 2.530, 1.553, 40.394, 0.232, 0.171, 0.099),
    ("Rb", 5, 1, None, None, 0.017, 0.000, 0.000, 40.393, 0.440, 0.333, 0.192),
    ("Sr", 6, 2, None, None, 0.008, 2.459, 1.510, None, None, None, None),
    ("Sr", 5, 3, None, None, 0.020, 3.585, 2.230, 45.518, 0.459, 0.350, 0.203),
    ("Sr", 6, 1, None, None, 0.016, 0.000, 0.000, None, None, None, None),
    ("Sr", 5, 2, 2.173, 9.473, 0.029, 4.784, 3.052, 45.518, 0.656, 0.508, 0.294),
    ("Sr", 5, 1, 1.579, 4.467, 0.038, 0.000, 0.000, 45.517, 0.857, 0.672, 0.389),
    ("Y", 6, 3, None, None, 0.019, 3.504, 2.177, None, None, None, None),
    ("Y", 5, 4, None, None, 0.033, 4.415, 2.774, 51.132, 0.666, 0.519, 0.301),
    ("Y", 6, 2, None, None, 0.028, 4.696, 2.992, None, None, None, None),
    ("Y", 5, 3, 1.689, 5.745, 0.042, 5.456, 3.497, 51.131, 0.855, 0.673, 0.390),
    ("Y", 6, 1, 1.620, 4.708, 0.038, 0.000, 0.000, None, None, None, None),
    ("Y", 5, 2, 1.306, 3.079, 0.052, 6.257, 4.067, 51.131, 1.047, 0.832, 0.482),
    ("Y", 5, 1, 1.089, 2.006, 0.062, 0.000, 0.000, 51.130, 1.238, 0.992, 0.575),
    ("Zr", 6, 4, None, None, 0.032, 4.338, 2.723, None, None, None, None),
    ("Zr", 5, 5, 1.949, 9.637, 0.047, 5.119, 3.238, 57.267, 0.863, 0.683, 0.396),
    ("Zr", 6, 3, 1.734, 6.062, 0.042, 5.379, 3.444, None, None, None, None),
    ("Zr", 5, 4, 1.381, 3.846, 0.058, 6.051, 3.892, 57.267, 1.046, 0.833, 0.483),
    ("Zr", 6, 2, 1.333, 3.213, 0.052, 6.188, 4.020, None, None, None, None),
    ("Zr", 5, 3, 1.111, 2.240, 0.068, 6.786, 4.418, 57.266, 1.232, 0.989, 0.574),
    ("Zr", 6, 1, 1.107, 2.076, 0.062, 0.000, 0.000, None, None, None, None),
    ("Zr", 5, 2, 0.947, 1.528, 0.079, 7.394, 4.859, 57.266, 1.419, 1.147, 0.666),
    ("Zr", 5, 1, 0.833, 1.135, 0.091, 0.000, 0.000, 57.266, 1.607, 1.309, 0.760),
    ("Nb", 6, 5, None, None, 0.047, 5.047, 3.190, None, None, None, None),
    ("Nb", 5, 6, 1.558, 6.078, 0.064, 5.745, 3.653, 63.959, 1.052, 0.842, 0.489),
    ("Nb", 6, 4, 1.411, 4.020, 0.057, 5.983, 3.844, None, None, None, None),
    ("Nb", 5, 5, 1.165, 2.738, 0.075, 6.600, 4.256, 63.958, 1.232, 0.991, 0.575),
    ("Nb", 6, 3, 1.130, 2.320, 0.068, 6.724, 4.375, None, None, None, None),
    ("Nb", 5, 4, 0.963, 1.691, 0.087, 7.288, 4.750, 63.958, 1.414, 1.144, 0.665),
    ("Nb", 6, 2, 0.960, 1.573, 0.079, 7.338, 4.820, None, None, None, None),
    ("Nb", 5, 3, 0.834, 1.194, 0.099, 7.866, 5.171, 63.958, 1.597, 1.302, 0.756),
    ("Nb", 6, 1, 0.843, 1.164, 0.091, 0.000, 0.000, None, None, None, None),
    ("Nb", 5, 2, 0.743, 0.907, 0.111, 8.370, 5.540, 63.958, 1.783, 1.462, 0.850),
    ("Mo", 6, 6, 0.000, 0.000, 0.000, 5.679, 3.608, None, None, None, None),
    ("Mo", 5, 7, 1.291, 4.126, 0.083, 6.320, 4.034, 71.243, 1.237, 0.998, 0.580),
    ("Mo", 6, 5, 1.186, 2.843, 0.075, 6.538, 4.213, None, None, None, None),
    ("Mo", 5, 6, 1.004, 2.034, 0.095, 7.118, 4.599, 71.242, 1.414, 1.147, 0.666),
    ("Mo", 6, 4, 0.977, 1.743, 0.087, 7.231, 4.711, None, None, None, None),
    ("Mo", 5, 5, 0.847, 1.313, 0.108, 7.769, 5.069, 71.242, 1.593, 1.299, 0.755),
    ("Mo", 6, 3, 0.845, 1.225, 0.099, 7.815, 5.135, None, None, None, None),
    ("Mo", 5, 4, 0.743, 0.953, 0.121, 8.325, 5.474, 71.241, 1.774, 1.455, 0.846),
    ("Mo", 6, 2, 0.750, 0.928, 0.112, 8.324, 5.508, None, None, None, None),
    ("Mo", 5, 3, 0.668, 0.738, 0.134, 8.814, 5.833, 71.242, 1.958, 1.615, 0.939),
    ("Tc", 6, 7, 1.316, 4.299, 0.083, 6.259, 3.992, None, None, None, None),
    ("Tc", 5, 8, 1.098, 2.951, 0.104, 6.860, 4.392, 79.157, 1.418, 1.153, 0.671),
    ("Tc", 6, 6, 1.019, 2.102, 0.095, 7.061, 4.560, None, None, None, None),
    ("Tc", 5, 7, 0.879, 1.560, 0.118, 7.613, 4.927, 79.156, 1.592, 1.300, 0.756),
    ("Tc", 6, 5, 0.858, 1.350, 0.108, 7.717, 5.033, None, None, None, None),
    ("Tc", 5, 6, 0.753, 1.044, 0.132, 8.236, 5.378, 79.155, 1.769, 1.452, 0.845),
    ("Tc", 6, 4, 0.751, 0.975, 0.121, 8.278, 5.441, None, None, None, None),
    ("Tc", 5, 5, 0.668, 0.774, 0.146, 8.773, 5.771, 79.155, 1.949, 1.607, 0.936),
    ("Tc", 6, 3, 0.674, 0.752, 0.135, 8.771, 5.803, None, None, None, None),
    ("Tc", 5, 4, 0.604, 0.608, 0.160, 9.250, 6.121, 79.156, 2.131, 1.767, 1.029),
    ("Ru", 6, 8, 1.116, 3.060, 0.104, 6.802, 4.352, None, None, None, None),
    ("Ru", 5, 9, 0.952, 2.193, 0.129, 7.373, 4.732, 87.739, 1.596, 1.306, 0.761),
    ("Ru", 6, 7, 0.890, 1.605, 0.118, 7.560, 4.890, None, None, None, None),
    ("Ru", 5, 8, 0.779, 1.226, 0.143, 8.090, 5.243, 87.738, 1.769, 1.453, 0.846),
    ("Ru", 6, 6, 0.762, 1.068, 0.132, 8.187, 5.344, None, None, None, None),
    ("Ru", 5, 7, 0.676, 0.844, 0.158, 8.690, 5.678, 87.737, 1.944, 1.604, 0.934),
    ("Ru", 6, 5, 0.674, 0.790, 0.146, 8.729, 5.739, None, None, None, None),
    ("Ru", 5, 6, 0.605, 0.637, 0.173, 9.213, 6.061, 87.737, 2.122, 1.759, 1.025),
    ("Ru", 6, 4, 0.609, 0.619, 0.161, 9.210, 6.093, None, None, None, None),
    ("Ru", 5, 5, 0.551, 0.508, 0.189, 9.680, 6.404, 87.738, 2.303, 1.918, 1.117),
    ("Rh", 6, 9, 0.965, 2.265, 0.129, 7.320, 4.695, None, None, None, None),
    ("Rh", 5, 10, 0.836, 1.678, 0.156, 7.866, 5.059, 97.031, 1.772, 1.458, 0.850),
    ("Rh", 6, 8, 0.787, 1.257, 0.144, 8.041, 5.209, None, None, None, None),
    ("Rh", 5, 9, 0.697, 0.982, 0.172, 8.554, 5.551, 97.030, 1.943, 1.605, 0.935),
    ("Rh", 6, 7, 0.683, 0.861, 0.159, 8.645, 5.647, None, None, None, None),
    ("Rh", 5, 8, 0.612, 0.692, 0.188, 9.135, 5.973, 97.029, 2.117, 1.755, 1.023),
    ("Rh", 6, 6, 0.610, 0.649, 0.175, 9.172, 6.031, None, None, None, None),
    ("Rh", 5, 7, 0.551, 0.531, 0.205, 9.645, 6.346, 97.029, 2.294, 1.910, 1.114),
    ("Rh", 6, 5, 0.554, 0.515, 0.191, 9.642, 6.377, None, None, None, None),
    ("Rh", 5, 6, 0.504, 0.428, 0.222, 10.104, 6.684, 97.030, 2.474, 2.068, 1.206),
    ("Pd", 6, 9, 0.704, 1.005, 0.173, 8.508, 5.518, None, None, None, None),
    ("Pd", 5, 10, 0.629, 0.800, 0.204, 9.007, 5.850, 107.075, 2.117, 1.756, 1.025),
    ("Pd", 6, 8, 0.617, 0.705, 0.189, 9.093, 5.943, None, None, None, None),
    ("Pd", 5, 9, 0.556, 0.575, 0.221, 9.572, 6.261, 107.074, 2.290, 1.907, 1.113),
    ("Pd", 6, 7, 0.555, 0.540, 0.206, 9.607, 6.319, None, None, None, None),
    ("Pd", 5, 8, 0.504, 0.447, 0.239, 10.072, 6.628, 107.074, 2.466, 2.061, 1.203),
    ("Pd", 6, 6, 0.507, 0.434, 0.223, 10.069, 6.659, None, None, None, None),
    ("Pd", 5, 7, 0.464, 0.364, 0.258, 10.522, 6.960, 107.075, 2.644, 2.219, 1.295),
    ("Ag", 6, 9, 0.561, 0.585, 0.223, 9.533, 6.233, None, None, None, None),
    ("Ag", 5, 10, 0.509, 0.483, 0.258, 10.002, 6.546, 117.915, 2.461, 2.057, 1.201),
    ("Ag", 6, 8, 0.508, 0.454, 0.241, 10.036, 6.602, None, None, None, None),
    ("Ag", 5, 9, 0.464, 0.380, 0.278, 10.493, 6.906, 117.915, 2.636, 2.210, 1.291),
    ("Ag", 6, 7, 0.467, 0.368, 0.260, 10.490, 6.936, None, None, None, None),
    ("Ag", 5, 8, 0.429, 0.312, 0.298, 10.938, 7.234, 117.916, 2.813, 2.368, 1.383),
    ("Cd", 6, 9, 0.467, 0.385, 0.280, 10.459, 6.881, None, None, None, None),
    ("Cd", 5, 10, 0.429, 0.326, 0.320, 10.910, 7.181, 129.599, 2.804, 2.359, 1.379),
    ("Cd", 6, 8, 0.431, 0.315, 0.300, 10.907, 7.211, None, None, None, None),
    ("Cd", 5, 9, 0.397, 0.269, 0.342, 11.350, 7.505, 129.600, 2.981, 2.516, 1.471),
)


APPENDIX_A_4D: dict[tuple[str, int, int], SlaterEntry] = {
    (element, p_electrons, d_electrons): SlaterEntry(
        element,
        p_electrons,
        d_electrons,
        r2,
        r4,
        zeta_d,
        fdd2,
        fdd4,
        zeta_2p,
        fpd2,
        gpd1,
        gpd3,
    )
    for (
        element,
        p_electrons,
        d_electrons,
        r2,
        r4,
        zeta_d,
        fdd2,
        fdd4,
        zeta_2p,
        fpd2,
        gpd1,
        gpd3,
    ) in APPENDIX_A_4D_ROWS
}


def register_appendix_a_3d(entry: SlaterEntry, *, overwrite: bool = False) -> None:
    """Register or extend an Appendix-A style 3d Slater table entry."""

    _register_entry(APPENDIX_A_3D, entry, overwrite=overwrite)


def register_appendix_a_4d(entry: SlaterEntry, *, overwrite: bool = False) -> None:
    """Register or extend an Appendix-A style 4d Slater table entry."""

    _register_entry(APPENDIX_A_4D, entry, overwrite=overwrite)


def _register_entry(
    table: dict[tuple[str, int, int], SlaterEntry],
    entry: SlaterEntry,
    *,
    overwrite: bool = False,
) -> None:
    key = _entry_key(entry.element, entry.p_electrons, entry.d_electrons)
    if not overwrite and key in table:
        raise KeyError(f"Appendix-A entry already exists for {key}")
    table[key] = SlaterEntry(
        element=key[0],
        p_electrons=entry.p_electrons,
        d_electrons=entry.d_electrons,
        r2=entry.r2,
        r4=entry.r4,
        zeta_d=entry.zeta_d,
        fdd2=entry.fdd2,
        fdd4=entry.fdd4,
        zeta_2p=entry.zeta_2p,
        fpd2=entry.fpd2,
        gpd1=entry.gpd1,
        gpd3=entry.gpd3,
    )


def available_appendix_a_3d(element: str | None = None) -> tuple[tuple[str, int, int], ...]:
    """Return available `(element, p_electrons, d_electrons)` table keys."""

    return _available_keys(APPENDIX_A_3D, element)


def available_appendix_a_4d(element: str | None = None) -> tuple[tuple[str, int, int], ...]:
    """Return available 4d `(element, p_electrons, d_electrons)` table keys."""

    return _available_keys(APPENDIX_A_4D, element)


def get_appendix_a_3d(element: str, p_electrons: int, d_electrons: int) -> SlaterEntry:
    return _get_entry(APPENDIX_A_3D, available_appendix_a_3d, "3d", element, p_electrons, d_electrons)


def get_appendix_a_4d(element: str, p_electrons: int, d_electrons: int) -> SlaterEntry:
    return _get_entry(APPENDIX_A_4D, available_appendix_a_4d, "4d", element, p_electrons, d_electrons)


def _get_entry(
    table: dict[tuple[str, int, int], SlaterEntry],
    available,
    label: str,
    element: str,
    p_electrons: int,
    d_electrons: int,
) -> SlaterEntry:
    key = _entry_key(element, p_electrons, d_electrons)
    try:
        return table[key]
    except KeyError as exc:
        available_text = ", ".join(f"{el} 2p{p} {label}{d}" for el, p, d in available(element))
        hint = f" Available for {key[0]}: {available_text}." if available_text else ""
        raise KeyError(f"No Appendix-A {label} entry for {key}.{hint}") from exc


def _entry_key(element: str, p_electrons: int, d_electrons: int) -> tuple[str, int, int]:
    return (element.capitalize(), p_electrons, d_electrons)


def _available_keys(
    table: dict[tuple[str, int, int], SlaterEntry],
    element: str | None = None,
) -> tuple[tuple[str, int, int], ...]:
    if element is None:
        keys = table.keys()
    else:
        normalized = element.capitalize()
        keys = (key for key in table if key[0] == normalized)
    return tuple(sorted(keys))

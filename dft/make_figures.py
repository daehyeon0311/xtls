"""Figures for the follow-up paper.

    python dft/make_figures.py

Two panels the argument rests on and neither exists yet: how D responds to the
two different deformations that get called "distortion", and how J responds to
the one the strain proposal actually means.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

OUTPUT = ROOT / "outputs" / "paper"

# Converged DFT, U = 5 eV, quantum convention
EXCHANGE = {6.000: 0.06620, 8.225: 0.09098, 10.000: 0.11110}
BOND_CAXIS = {6.000: 2.0292, 8.225: 1.9875, 10.000: 1.9575}
CRITICAL = 19.0
BFSO_ANGLE = 8.225


def load_scan(name: str):
    path = ROOT / "outputs" / "anisotropy" / name
    lines = path.read_text(encoding="utf-8").splitlines()
    keys = lines[0].split()
    rows = [dict(zip(keys, map(float, line.split()))) for line in lines[1:]]
    return (
        np.array([row["delta_theta_deg"] for row in rows]),
        np.array([row["D_meV"] for row in rows]),
        np.array([row["bond_length_A"] for row in rows]),
    )


def figure_geometry(plt) -> None:
    """D under the two deformations that both get called 'distortion'."""
    angles_p, d_poisson, bond_p = load_scan("anisotropy_scan.txt")
    angles_c, d_caxis, bond_c = load_scan("anisotropy_scan_caxis.txt")

    figure, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))

    axes[0].plot(angles_p, bond_p, lw=1.8, color="tab:orange",
                 label="A-site substitution\n(interpolating BFSO$\\leftrightarrow$SFSO)")
    axes[0].plot(angles_c, bond_c, lw=1.8, color="tab:blue",
                 label="uniaxial $c$ strain")
    axes[0].set_ylabel(r"Fe$-$O distance ($\mathrm{\AA}$)")
    axes[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    axes[0].set_title("the two deformations move Fe$-$O oppositely", fontsize=9)

    axes[1].plot(angles_p, d_poisson, lw=1.8, color="tab:orange")
    axes[1].plot(angles_c, d_caxis, lw=1.8, color="tab:blue")
    axes[1].set_ylabel("$D$ (meV)")
    axes[1].set_title("and so give different $D(\\Delta\\theta)$", fontsize=9)
    axes[1].annotate(f"{100 * (d_poisson[-1] / d_poisson[0] - 1):+.0f}%",
                     (angles_p[-1], d_poisson[-1]), textcoords="offset points",
                     xytext=(-34, 4), fontsize=8, color="tab:orange")
    axes[1].annotate(f"{100 * (d_caxis[-1] / d_caxis[0] - 1):+.0f}%",
                     (angles_c[-1], d_caxis[-1]), textcoords="offset points",
                     xytext=(-30, -12), fontsize=8, color="tab:blue")

    for axis in axes:
        axis.axvline(BFSO_ANGLE, color="0.6", ls=":", lw=1.0)
        axis.set_xlabel(r"distortion angle $\Delta\theta$ (deg)")
        axis.tick_params(direction="in", top=True, right=True)
    figure.tight_layout()
    figure.savefig(OUTPUT / "fig_geometry.png", dpi=300)
    plt.close(figure)


def figure_exchange(plt) -> None:
    """J(distortion) and the resulting phase boundary.

    Restricted to the range where J was actually computed. Extrapolating it
    would put a flat segment on the plot that reads as data.
    """
    known = np.array(sorted(EXCHANGE))
    values = np.array([EXCHANGE[a] for a in known])

    angles_all, d_all, _bond = load_scan("anisotropy_scan_caxis.txt")
    inside = (angles_all >= known[0]) & (angles_all <= known[-1])
    angles, d_values = angles_all[inside], d_all[inside]
    j_curve = np.interp(angles, known, values)

    figure, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))

    # Relative change makes the disparity legible; the absolute scales differ
    # by more than an order of magnitude.
    reference = float(np.interp(BFSO_ANGLE, angles, j_curve))
    d_reference = float(np.interp(BFSO_ANGLE, angles, d_values))
    axes[0].plot(angles, 100 * (j_curve / reference - 1), lw=1.8, color="tab:red", label="$J$ (DFT)")
    axes[0].plot(known, 100 * (values / reference - 1), "o", ms=6, color="tab:red")
    axes[0].plot(angles, 100 * (d_values / d_reference - 1), lw=1.8, color="tab:blue",
                 label="$D$ (cluster model)")
    axes[0].axhline(0.0, color="0.75", lw=0.8)
    axes[0].set_ylabel("change relative to BFSO (%)")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    change_j = 100 * (values[-1] / values[0] - 1)
    change_d = 100 * (d_values[-1] / d_values[0] - 1)
    axes[0].set_title(f"$J$ varies {change_j:+.0f}% where $D$ varies {change_d:+.0f}%", fontsize=9)

    ratio = d_values / j_curve
    axes[1].fill_between(angles, ratio, CRITICAL, where=ratio >= CRITICAL,
                         color="tab:purple", alpha=0.13)
    axes[1].plot(angles, ratio, lw=2.0, color="black")
    axes[1].axhline(CRITICAL, color="tab:red", ls="--", lw=1.2)
    axes[1].text(angles[0] + 0.08, CRITICAL * 1.02, r"$(D/J)_c$", fontsize=9, color="tab:red")
    axes[1].text(angles[0] + 0.08, CRITICAL * 1.11, "quantum paramagnet", fontsize=8, color="0.3")
    axes[1].text(angles[-1] - 1.9, CRITICAL * 0.76, "antiferromagnet", fontsize=8, color="0.3")
    axes[1].set_ylabel("$D/J$")
    axes[1].set_title("boundary lies toward $tensile$ strain", fontsize=9)

    crossing = float(np.interp(CRITICAL, ratio[::-1], angles[::-1]))
    axes[1].plot([crossing], [CRITICAL], "o", ms=6, color="tab:red")
    axes[1].annotate(f"{crossing:.2f}$^\\circ$", (crossing, CRITICAL),
                     textcoords="offset points", xytext=(6, -14), fontsize=8)

    for axis in axes:
        axis.axvline(BFSO_ANGLE, color="0.6", ls=":", lw=1.0)
        axis.set_xlabel(r"distortion angle $\Delta\theta$ (deg)")
        axis.tick_params(direction="in", top=True, right=True)
    axes[1].annotate("BFSO", (BFSO_ANGLE, ratio.max() * 0.98),
                     textcoords="offset points", xytext=(4, -4), fontsize=8, color="0.35")
    figure.tight_layout()
    figure.savefig(OUTPUT / "fig_exchange.png", dpi=300)
    plt.close(figure)
    print(f"  critical crossing at {crossing:.2f} deg")


def figure_family(plt) -> None:
    """D across the melilite family, split by whether a QPM is possible at all."""
    ions = [
        {"label": "Mn$^{2+}$\nd$^5$, S=5/2", "d": 5, "D": 0.0018, "kramers": True, "note": "$^6A_1$"},
        {"label": "Fe$^{2+}$\nd$^6$, S=2", "d": 6, "D": 1.539, "kramers": False, "note": "$^5E$"},
        {"label": "Co$^{2+}$\nd$^7$, S=3/2", "d": 7, "D": 3.496, "kramers": True, "note": "$^4A_2$"},
        {"label": "Ni$^{2+}$\nd$^8$, S=1", "d": 8, "D": None, "kramers": False, "note": "$^3T_1$"},
    ]
    figure, axis = plt.subplots(figsize=(6.2, 4.0))
    for index, ion in enumerate(ions):
        if ion["D"] is None:
            axis.bar(index, 3.9, color="0.85", hatch="//", edgecolor="0.6")
            axis.text(index, 2.0, "orbital triplet\nno spin\nHamiltonian",
                      ha="center", va="center", fontsize=7.5, color="0.35")
            continue
        colour = "tab:blue" if not ion["kramers"] else "0.7"
        axis.bar(index, ion["D"], color=colour, edgecolor="0.3", lw=0.6)
        text = f"{ion['D']:.3f}" if ion["D"] < 0.01 else f"{ion['D']:.2f}"
        axis.text(index, ion["D"] + 0.12, text, ha="center", fontsize=8)
        if ion["kramers"]:
            axis.text(index, ion["D"] / 2, "Kramers\nno QPM", ha="center", va="center",
                      fontsize=7.5, color="0.25")
    axis.set_xticks(range(len(ions)))
    axis.set_xticklabels([ion["label"] for ion in ions], fontsize=8)
    for index, ion in enumerate(ions):
        axis.text(index, -0.30, ion["note"], ha="center", fontsize=8, color="0.4")
    axis.set_ylabel("$D$ (meV), at the BFSO geometry")
    axis.set_ylim(-0.45, 4.2)
    axis.axhline(0.0, color="0.5", lw=0.8)
    axis.set_title("only integer spin admits a non-magnetic singlet ground state", fontsize=9)
    axis.tick_params(direction="in", right=True)
    figure.tight_layout()
    figure.savefig(OUTPUT / "fig_family.png", dpi=300)
    plt.close(figure)


def figure_convention(plt) -> None:
    """D(distortion) under the two crystal-field normalization conventions."""
    import run_anisotropy as anisotropy
    import run_xas as runner
    from xtls_py.geometry import _point_charge_field, t2_e_splitting

    runner._load_input_file(runner.ROOT / "inputs" / "Fe_Ba2FeSi2O7.py")
    runner.__dict__["n_analyzed_states"] = 8

    def sign_based_splitting(matrix) -> float:
        """The superseded convention: group levels by the sign of their energy."""
        levels = np.sort(np.real(np.linalg.eigvals(matrix)))
        positive = int(np.sum(levels > 0))
        if positive in (0, 5):
            return float("nan")
        upper = np.sum(levels[5 - positive:]) / positive
        lower = np.sum(levels[: 5 - positive]) / (5 - positive)
        return abs(upper - lower)

    angles = np.arange(4.0, 12.01, 0.25)
    d_sign, d_symmetry = [], []
    for angle in angles:
        bond = anisotropy.bond_length_for(angle)
        runner.__dict__["ligand_angle_offset_deg"] = float(angle)
        runner.__dict__["ligand_radius"] = float(bond)
        positions = runner._ligand_positions()
        raw, _qkm = _point_charge_field(positions, runner.r2, runner.r4)

        scale = runner.ten_dq / sign_based_splitting(raw)
        levels, _sz = anisotropy.multiplet_levels(angle, bond_length_A=bond, scale=scale)
        d_sign.append(anisotropy.fit_anisotropy(levels)["D_meV"])

        scale = runner.ten_dq / abs(t2_e_splitting(raw))
        levels, _sz = anisotropy.multiplet_levels(angle, bond_length_A=bond, scale=scale)
        d_symmetry.append(anisotropy.fit_anisotropy(levels)["D_meV"])

    figure, axis = plt.subplots(figsize=(6.2, 4.0))
    axis.plot(angles, d_sign, lw=1.6, color="tab:red", label="grouped by sign of level energy")
    axis.plot(angles, d_symmetry, lw=1.8, color="tab:blue", label="grouped by cubic symmetry label")
    axis.set_xlabel(r"distortion angle $\Delta\theta$ (deg)")
    axis.set_ylabel("$D$ (meV)")
    axis.legend(frameon=False, fontsize=8, loc="upper left")
    axis.set_title("crystal-field normalization: level crossings break the sign rule", fontsize=9)
    axis.tick_params(direction="in", top=True, right=True)
    figure.tight_layout()
    figure.savefig(OUTPUT / "fig_convention.png", dpi=300)
    plt.close(figure)
    jumps = np.abs(np.diff(d_sign))
    print(f"  sign-based convention: largest step between adjacent points "
          f"{np.nanmax(jumps):.4f} meV, symmetry-based {np.max(np.abs(np.diff(d_symmetry))):.4f}")


def main() -> None:
    import matplotlib.pyplot as plt

    OUTPUT.mkdir(parents=True, exist_ok=True)
    figure_levels(plt)
    figure_geometry(plt)
    figure_exchange(plt)
    figure_family(plt)
    figure_convention(plt)
    print(f"wrote figures to {OUTPUT}")



def figure_levels(plt) -> None:
    """FeO4 geometry and the level splitting that produces D."""
    from matplotlib.patches import Arc

    figure = plt.figure(figsize=(9.8, 4.0))
    left = figure.add_subplot(1, 2, 1)
    right = figure.add_subplot(1, 2, 2)

    # -- (a) the tetrahedron, projected on the xz plane -----------------------
    bond, theta = 1.9875, np.radians(54.7356 + 8.225)
    x, z = bond * np.sin(theta), bond * np.cos(theta)
    ideal = np.radians(54.7356)
    xi, zi = bond * np.sin(ideal), bond * np.cos(ideal)
    for sign_x, sign_z in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
        left.plot([0, sign_x * xi], [0, sign_z * zi], color="0.78", lw=1.1,
                  ls="--", zorder=0)
    left.plot(xi, zi, "o", ms=8, mfc="none", mec="0.7", mew=1.1, zorder=0)
    left.text(xi + 0.10, zi + 0.16, "ideal", fontsize=7.5, color="0.55")
    for sign_x, sign_z in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
        left.plot([0, sign_x * x], [0, sign_z * z], color="0.35", lw=1.4, zorder=1)
        left.plot(sign_x * x, sign_z * z, "o", ms=11, color="tab:red", zorder=2)
    left.plot(0, 0, "o", ms=13, color="tab:blue", zorder=3)
    left.text(0.13, 0.10, "Fe", fontsize=9, color="tab:blue")
    left.text(x + 0.12, z, "O", fontsize=9, color="tab:red")
    left.axvline(0, color="0.75", ls=":", lw=1.0)
    left.annotate("", xy=(0, 1.55), xytext=(0, -1.55),
                  arrowprops=dict(arrowstyle="<->", color="0.6", lw=0.9))
    left.text(0.08, 1.35, "$c$", fontsize=9, color="0.45")
    left.add_patch(Arc((0, 0), 1.5, 1.5, angle=0, theta1=90 - np.degrees(theta),
                       theta2=90, color="tab:green", lw=1.4))
    left.text(0.30, 0.86, r"$\theta$", fontsize=10, color="tab:green")
    left.text(0, -1.85, r"$\Delta\theta = \theta - 54.74^\circ = 8.23^\circ$" "\n"
              r"Fe$-$O $= 1.988\ \mathrm{\AA}$", ha="center", fontsize=8.5, color="0.3")
    left.set_xlim(-2.4, 2.4); left.set_ylim(-2.3, 1.9)
    left.set_aspect("equal"); left.axis("off")
    left.set_title("compressed FeO$_4$ tetrahedron", fontsize=9)

    # -- (b) level splitting, with the computed numbers -----------------------
    stages = [
        (0.00, 0.55, [(0.0, r"$^5D$ (25)")], "free ion"),
        (0.85, 1.40, [(0.55, r"$^5T_2$ (15)"), (-0.55, r"$^5E$ (10)")], "$T_d$ field"),
        (1.70, 2.25, [(0.85, r"$^5B$"), (-0.85, r"$^5A$")], "compression"),
    ]
    for start, stop, levels, caption in stages:
        for height, name in levels:
            right.plot([start, stop], [height, height], lw=2.0, color="0.25")
            right.text(stop + 0.04, height, name, fontsize=8.5, va="center")
        right.text((start + stop) / 2, -1.95, caption, ha="center", fontsize=8.5, color="0.4")
    for a, b in ((0.55, 0.85), (1.40, 1.70)):
        right.plot([a, b], [0.0, 0.55], color="0.7", lw=0.8, ls="--")
        right.plot([a, b], [0.0, -0.55], color="0.7", lw=0.8, ls="--")
    right.plot([1.40, 1.70], [-0.55, -0.85], color="0.7", lw=0.8, ls="--")
    right.plot([1.40, 1.70], [-0.55, 0.85], color="0.7", lw=0.8, ls="--")

    # spin-orbit split manifold, drawn to scale from the computed energies
    base, span = -0.85, 1.30
    for energy, label in ((0.0, "$S_z = 0$"), (1.4517, r"$S_z = \pm 1$"), (5.8269, r"$S_z = \pm 2$")):
        height = base + span * energy / 5.8269
        right.plot([2.60, 3.15], [height, height], lw=2.2, color="tab:blue")
        right.text(3.30, height, f"{label}   {energy:.2f} meV", fontsize=8.5, va="center")
        right.plot([2.25, 2.60], [-0.85, height], color="0.7", lw=0.8, ls="--")
    right.text(2.875, -1.95, "spin$-$orbit", ha="center", fontsize=8.5, color="0.4")
    right.annotate("", xy=(3.15, base), xytext=(3.15, base + span * 1.4517 / 5.8269),
                   arrowprops=dict(arrowstyle="<->", color="tab:blue", lw=1.1))
    right.text(3.06, base + span * 0.72 / 5.8269, "$D$", fontsize=9.5,
               color="tab:blue", ha="right", va="center")
    right.set_xlim(-0.15, 4.6); right.set_ylim(-2.2, 1.75)
    right.axis("off")
    right.set_title("$S = 2$ manifold, energies from this work", fontsize=9)

    figure.tight_layout()
    figure.savefig(OUTPUT / "fig_levels.png", dpi=300)
    plt.close(figure)

if __name__ == "__main__":
    main()

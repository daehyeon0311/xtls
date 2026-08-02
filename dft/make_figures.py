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


def main() -> None:
    import matplotlib.pyplot as plt

    OUTPUT.mkdir(parents=True, exist_ok=True)
    figure_geometry(plt)
    figure_exchange(plt)
    print(f"wrote figures to {OUTPUT}")


if __name__ == "__main__":
    main()

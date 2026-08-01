"""D/J phase diagram against tetrahedral distortion.

    python dft/phase_diagram.py

Combines the two halves of the follow-up: D from the charge-transfer cluster
model (run_anisotropy.py) and J from broken-symmetry DFT (extract_exchange.py).
The XLD paper computed D and took J from neutron scattering, holding it fixed
while the distortion varied. With J(distortion) available the boundary can be
placed properly, since D/J is what decides the ground state, not D alone.

Without DFT J(distortion) yet, this falls back to a constant J and reports
what that assumption implies, which is exactly the paper's position.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

ALPHA_C = 6.35          # critical D/J, from the XLD paper (1.65 / 0.26)
J_NEUTRON = 0.26        # meV, Do et al.
J_DFT_U3 = 0.2329       # meV, this work at U = 3 eV


def load_anisotropy() -> list[dict[str, float]]:
    """D(distortion) from the cluster-model scan."""
    path = ROOT / "outputs" / "anisotropy" / "anisotropy_scan.txt"
    if not path.exists():
        raise SystemExit(f"run run_anisotropy.py first: {path} is missing")
    lines = path.read_text(encoding="utf-8").splitlines()
    keys = lines[0].split()
    rows = []
    for line in lines[1:]:
        values = [float(v) for v in line.split()]
        rows.append(dict(zip(keys, values)))
    return rows


def load_exchange() -> dict[float, float] | None:
    """J(distortion) from DFT, if the strained relaxations have been analysed."""
    path = Path(__file__).resolve().parent / "exchange_vs_strain.json"
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {float(k): float(v) for k, v in payload.items()}


def crossing(angles, values, target: float) -> float | None:
    angles, values = np.asarray(angles), np.asarray(values)
    above = np.where(values >= target)[0]
    if not above.size or above[0] == 0:
        return None
    i = above[0]
    span = values[i] - values[i - 1]
    if abs(span) < 1e-12:
        return float(angles[i])
    weight = (target - values[i - 1]) / span
    return float(angles[i - 1] + weight * (angles[i] - angles[i - 1]))


def main() -> None:
    rows = load_anisotropy()
    angles = np.array([row["delta_theta_deg"] for row in rows])
    d_values = np.array([row["D_meV"] for row in rows])
    exchange = load_exchange()

    print("D from the cluster model, J from DFT (U = 3 eV)\n")
    print("  scenario                     J (meV)   D/J at 7.84 deg   crosses alpha_c at")
    scenarios = {
        "neutron J, held fixed": J_NEUTRON,
        "DFT J (U=3), held fixed": J_DFT_U3,
    }
    reference = float(np.interp(7.8353, angles, d_values))
    for label, j_value in scenarios.items():
        ratio = d_values / j_value
        where = crossing(angles, ratio, ALPHA_C)
        print(
            f"  {label:<28} {j_value:.4f}    {reference / j_value:6.2f}"
            f"          {f'{where:.2f} deg' if where else 'not reached'}"
        )

    if exchange is None:
        print("\n  J(distortion) from DFT not available yet "
              "(dft/exchange_vs_strain.json missing).")
        print("  Holding J fixed is the paper's own assumption; the strained-cell")
        print("  relaxations are what replace it.")
    else:
        print("\n  J(distortion) from DFT:")
        for strain, j_value in sorted(exchange.items()):
            print(f"    strain {strain:+.3f}  J = {j_value:.4f} meV")

    plot(angles, d_values, scenarios, reference)


def plot(angles, d_values, scenarios, reference) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    axes[0].plot(angles, d_values, color="crimson", lw=1.8)
    axes[0].set_ylabel("$D$ (meV)")
    axes[0].axvline(7.8353, color="0.6", ls=":", lw=1.0)
    axes[0].text(7.95, d_values.min(), "BFSO", fontsize=8, color="0.35")

    for (label, j_value), style in zip(scenarios.items(), ["-", "--"]):
        axes[1].plot(angles, d_values / j_value, style, lw=1.6, label=label)
    axes[1].axhline(ALPHA_C, color="0.4", ls=":", lw=1.2)
    axes[1].text(angles[0] + 0.1, ALPHA_C * 1.02, r"$\alpha_c$", fontsize=9, color="0.3")
    axes[1].axvline(7.8353, color="0.6", ls=":", lw=1.0)
    axes[1].set_ylabel("$D/J$")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].text(angles[0] + 0.1, ALPHA_C * 1.25, "quantum paramagnet", fontsize=8, color="0.35")
    axes[1].text(angles[0] + 0.1, ALPHA_C * 0.72, "antiferromagnet", fontsize=8, color="0.35")

    for axis in axes:
        axis.set_xlabel(r"distortion angle $\Delta\theta$ (deg)")
        axis.tick_params(direction="in", top=True, right=True)
    figure.tight_layout()
    output = ROOT / "outputs" / "anisotropy" / "phase_diagram.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=200)
    plt.close(figure)
    print(f"\nsaved: {output}")


if __name__ == "__main__":
    main()

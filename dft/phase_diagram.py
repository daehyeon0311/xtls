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

# Conventions. Do et al. (Nat. Commun. 12, 5331) quote an effective S = 1
# Hamiltonian obtained by projecting out the Sz = +-2 doublet, and state the
# mapping explicitly: J~ = 3J, D~ = D. Everything here is converted to the
# S = 2 Hamiltonian so that DFT energies can be compared directly.
ALPHA_C_EFFECTIVE = 0.158            # critical J~/D~ from Do et al.
CRITICAL_D_OVER_J = 3.0 / ALPHA_C_EFFECTIVE   # = 19.0 in the S = 2 convention

J_NEUTRON = 0.266 / 3.0              # meV, Do et al. converted to S = 2
D_NEUTRON = 1.42                     # meV, D~ = D so this needs no conversion

# Broken-symmetry DFT, E_FM - E_AFM = 8 J S(S+1) with the quantum convention,
# which is what reproduces the measured exchange.
J_DFT_U5 = 4.367 / 48.0              # meV, this work at U = 5 eV


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
    """J(distortion) from DFT, keyed by distortion angle in degrees."""
    path = Path(__file__).resolve().parent / "exchange_vs_distortion.json"
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {float(k): float(v) for k, v in payload.items()}


def interpolate_exchange(angles: np.ndarray, exchange: dict[float, float]) -> np.ndarray:
    """J across the scan, linearly through the computed distortion points."""
    known_angles = np.array(sorted(exchange))
    known_values = np.array([exchange[a] for a in known_angles])
    return np.interp(angles, known_angles, known_values)


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

    print("D from the cluster model, J from DFT, both in the S = 2 convention")
    print(f"critical (D/J)_c = 3/alpha_c = {CRITICAL_D_OVER_J:.1f}\n")
    print("  scenario                     J (meV)   D/J at 7.84 deg   crosses (D/J)_c at")
    scenarios = {
        "neutron J (Do et al.)": J_NEUTRON,
        "DFT J (U = 5 eV)": J_DFT_U5,
    }
    reference = float(np.interp(7.8353, angles, d_values))
    for label, j_value in scenarios.items():
        ratio = d_values / j_value
        where = crossing(angles, ratio, CRITICAL_D_OVER_J)
        print(
            f"  {label:<28} {j_value:.4f}    {reference / j_value:6.2f}"
            f"          {f'{where:.2f} deg' if where else 'beyond 12 deg'}"
        )
    print(f"\n  measured D/J (Do et al.) = {D_NEUTRON / J_NEUTRON:.1f}, "
          f"i.e. {100 * (D_NEUTRON / J_NEUTRON) / CRITICAL_D_OVER_J:.0f}% of critical")

    varying = None
    if exchange is None:
        print("\n  J(distortion) from DFT not available yet "
              "(dft/exchange_vs_distortion.json missing).")
        print("  Holding J fixed is the paper's own assumption.")
    else:
        print("\n  J(distortion) from DFT (U = 5 eV):")
        for angle, j_value in sorted(exchange.items()):
            print(f"    {angle:6.3f} deg   J = {j_value:.4f} meV")
        varying = interpolate_exchange(angles, exchange)
        ratio = d_values / varying
        where = crossing(angles, ratio, CRITICAL_D_OVER_J)
        span = sorted(exchange)
        change = 100.0 * (exchange[span[-1]] / exchange[span[0]] - 1.0)
        print(f"\n  J varies by {change:+.1f}% over the scan, which the paper assumed away.")
        print(f"  With J(distortion) included, D/J crosses the critical ratio at "
              f"{f'{where:.2f} deg' if where else 'beyond the scan'}.")
        fixed = crossing(angles, d_values / J_DFT_U5, CRITICAL_D_OVER_J)
        if where and fixed:
            print(f"  Holding J fixed at its BFSO value would give {fixed:.2f} deg, "
                  f"a shift of {where - fixed:+.2f} deg.")

    plot(angles, d_values, scenarios, reference, varying)


def plot(angles, d_values, scenarios, reference, varying=None) -> None:
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
        axes[1].plot(angles, d_values / j_value, style, lw=1.2, alpha=0.7, label=f"{label}, J fixed")
    if varying is not None:
        axes[1].plot(angles, d_values / varying, "-", lw=2.0, color="black",
                     label="DFT $J(\\Delta\\theta)$")
    axes[1].axhline(CRITICAL_D_OVER_J, color="0.4", ls=":", lw=1.2)
    axes[1].text(angles[0] + 0.1, CRITICAL_D_OVER_J * 1.01, r"$(D/J)_c$", fontsize=9, color="0.3")
    axes[1].axvline(7.8353, color="0.6", ls=":", lw=1.0)
    axes[1].set_ylabel("$D/J$")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].text(angles[0] + 0.1, CRITICAL_D_OVER_J * 1.04, "quantum paramagnet", fontsize=8, color="0.35")
    axes[1].text(angles[0] + 0.1, CRITICAL_D_OVER_J * 0.80, "antiferromagnet", fontsize=8, color="0.35")

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

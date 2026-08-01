"""Exchange coupling from the Quantum ESPRESSO runs.

    python dft/extract_exchange.py

Reads the ferromagnetic and Neel total energies for each Hubbard U, converts
them to J, and compares against the two literature values, which disagree with
each other by a factor of two.

    H = J sum_<ij> S_i . S_j      4 nearest-neighbour bonds per cell
    E_FM - E_AFM = 8 J S^2        S = 2

The classical S^2 convention is the usual one for broken-symmetry DFT; the
S(S+1) form is reported alongside since the choice shifts J by 3/2.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REMOTE = "~/bfso"

RY_MEV = 13605.693122994
SPIN = 2.0
BONDS_FACTOR = 8.0  # E_FM - E_AFM = 8 J S^2

LITERATURE = {
    "neutron (Do et al.)": 0.26,
    "Monte Carlo (Jang et al.)": 0.14,
}
ALPHA_C = 6.35  # critical D/J quoted in the XLD paper (1.65 / 0.26)
D_CLUSTER = 1.45  # our calibrated D for BFSO, in meV


def wsl_read(name: str) -> str | None:
    """Fetch one output file from the WSL working directory."""
    result = subprocess.run(
        ["wsl", "-e", "bash", "-lc", f"cat {REMOTE}/{name} 2>/dev/null"],
        capture_output=True,
        text=True,
        errors="replace",
    )
    return result.stdout or None


def total_energy(text: str | None) -> float | None:
    """Final converged total energy in Ry, or None if the run did not finish."""
    if not text or "convergence has been achieved" not in text:
        return None
    matches = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", text, re.MULTILINE)
    return float(matches[-1]) if matches else None


def magnetisation(text: str | None) -> tuple[float, float] | None:
    if not text:
        return None
    total = re.findall(r"total magnetization\s+=\s+(-?\d+\.\d+)", text)
    absolute = re.findall(r"absolute magnetization\s+=\s+(-?\d+\.\d+)", text)
    if not total or not absolute:
        return None
    return float(total[-1]), float(absolute[-1])


def main() -> None:
    rows = []
    for hubbard_u in (3.0, 4.0, 5.0, 6.0):
        tag = f"u{hubbard_u:g}".replace(".", "p")
        fm_text = wsl_read(f"scf_fm_{tag}.out")
        afm_text = wsl_read(f"scf_afm_{tag}.out")
        e_fm, e_afm = total_energy(fm_text), total_energy(afm_text)
        if e_fm is None or e_afm is None:
            rows.append({"U": hubbard_u, "status": "pending"})
            continue
        difference = (e_fm - e_afm) * RY_MEV
        rows.append(
            {
                "U": hubbard_u,
                "status": "done",
                "E_FM": e_fm,
                "E_AFM": e_afm,
                "dE_meV": difference,
                "J_meV": difference / (BONDS_FACTOR * SPIN**2),
                "J_quantum_meV": difference / (BONDS_FACTOR * SPIN * (SPIN + 1)),
                "mag_fm": magnetisation(fm_text),
                "mag_afm": magnetisation(afm_text),
            }
        )

    print("  U (eV)   E_FM (Ry)        E_AFM (Ry)       dE (meV)   J (meV)   D/J    ground state")
    for row in rows:
        if row["status"] != "done":
            print(f"  {row['U']:5.1f}   {'-- pending --':^47}")
            continue
        ratio = D_CLUSTER / row["J_meV"] if row["J_meV"] else float("nan")
        state = "AFM" if row["dE_meV"] > 0 else "FM"
        verdict = "QPM predicted" if ratio > ALPHA_C else "AFM, below critical"
        print(
            f"  {row['U']:5.1f}   {row['E_FM']:.8f}  {row['E_AFM']:.8f}  "
            f"{row['dE_meV']:+8.3f}   {row['J_meV']:6.4f}   {ratio:5.2f}  {state}, {verdict}"
        )

    done = [row for row in rows if row["status"] == "done"]
    if len(done) >= 2:
        print()
        print("literature:")
        for label, value in LITERATURE.items():
            print(f"  {label:<28} J = {value:.3f} meV")
        print()
        slope = (done[-1]["J_meV"] - done[0]["J_meV"]) / (done[-1]["U"] - done[0]["U"])
        print(f"  dJ/dU = {slope:+.4f} meV/eV")
        for label, target in LITERATURE.items():
            if abs(slope) > 1e-9:
                implied = done[0]["U"] + (target - done[0]["J_meV"]) / slope
                print(f"  U reproducing {label:<28} = {implied:5.2f} eV")
        print()
        print(f"  critical ratio alpha_c = {ALPHA_C:.2f}, our D = {D_CLUSTER:.2f} meV")
        print(f"  -> J must exceed {D_CLUSTER / ALPHA_C:.3f} meV for an ordered ground state")

    if any(row["status"] == "done" and row["mag_fm"] for row in done):
        print()
        print("magnetisation check (should be 8.00 / 0.00 total, ~8 absolute):")
        for row in done:
            print(
                f"  U={row['U']:.0f}  FM total {row['mag_fm'][0]:+.2f} abs {row['mag_fm'][1]:.2f}"
                f"   |  AFM total {row['mag_afm'][0]:+.2f} abs {row['mag_afm'][1]:.2f}"
            )


if __name__ == "__main__":
    main()

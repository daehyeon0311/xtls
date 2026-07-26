from __future__ import annotations

import re
from dataclasses import dataclass

from .basis import FockBasis


_P_SPIN_ORBITALS = 6
_D_SPIN_ORBITALS = 10


@dataclass(frozen=True)
class Configuration:
    """Electron-count description for one configuration sector.

    The current engine targets transition-metal L-edge XAS/XLD, so the
    table-driven shell layout is `2p + 3d + ligand`. A ligand hole count `h`
    represents the compact charge-transfer notation `L^h`.
    """

    d_electrons: int
    p_electrons: int | None = None
    ligand_holes: int = 0
    tag: str = ""

    def __post_init__(self) -> None:
        if not 0 <= self.d_electrons <= _D_SPIN_ORBITALS:
            raise ValueError("d_electrons must satisfy 0 <= n <= 10")
        if self.p_electrons is not None and not 0 <= self.p_electrons <= _P_SPIN_ORBITALS:
            raise ValueError("p_electrons must satisfy 0 <= n <= 6")
        if self.ligand_holes < 0:
            raise ValueError("ligand_holes must be non-negative")
        if self.ligand_holes > _D_SPIN_ORBITALS:
            raise ValueError("ligand_holes must fit the ligand shell")

    @property
    def ligand_electrons(self) -> int:
        return _D_SPIN_ORBITALS - self.ligand_holes

    @property
    def core_holes(self) -> int:
        if self.p_electrons is None:
            return 0
        return _P_SPIN_ORBITALS - self.p_electrons

    @property
    def has_core_shell(self) -> bool:
        return self.p_electrons is not None

    @property
    def n_electrons_without_core(self) -> int:
        return self.d_electrons + self.ligand_electrons

    @property
    def n_electrons(self) -> int:
        return self.n_electrons_without_core + (self.p_electrons or 0)

    def label(self) -> str:
        pieces: list[str] = []
        if self.p_electrons is not None:
            pieces.append(f"2p{self.p_electrons}")
        pieces.append(f"3d{self.d_electrons}")
        if self.ligand_holes == 1:
            pieces.append("L")
        elif self.ligand_holes > 1:
            pieces.append(f"L{self.ligand_holes}")
        if self.tag:
            pieces.append(self.tag)
        return " ".join(pieces)


def parse_configuration(text: str) -> Configuration:
    """Parse compact XTLS-like labels such as `2p6 3d6` or `2p5 3d8L2`."""

    normalized = text.strip().replace("^", "")
    if not normalized:
        raise ValueError("configuration string is empty")
    p_match = re.search(r"(?<![A-Za-z0-9])2p\s*(\d+)(?!\d)", normalized, flags=re.IGNORECASE)
    d_match = re.search(r"(?<![A-Za-z0-9])3d\s*(\d+)(?!\d)", normalized, flags=re.IGNORECASE)
    ligand_match = re.search(r"(?<![A-Za-z0-9])L\s*(\d*)(?!\d)", normalized, flags=re.IGNORECASE)
    if ligand_match is None:
        ligand_match = re.search(r"(?<=\d)L\s*(\d*)(?!\d)", normalized, flags=re.IGNORECASE)
    if d_match is None:
        raise ValueError(f"configuration must contain a 3d electron count: {text!r}")
    ligand_holes = 0
    if ligand_match is not None:
        ligand_holes = int(ligand_match.group(1) or "1")
    return Configuration(
        p_electrons=None if p_match is None else int(p_match.group(1)),
        d_electrons=int(d_match.group(1)),
        ligand_holes=ligand_holes,
    )


def xas_initial_configurations(
    n_d_electrons: int,
    max_ligand_holes: int = 1,
) -> tuple[Configuration, ...]:
    """Return `2p6 3d^n`, `2p6 3d^(n+1)L`, ... sectors."""

    _validate_hole_request(n_d_electrons, max_ligand_holes)
    max_holes = min(max_ligand_holes, _D_SPIN_ORBITALS - n_d_electrons)
    return tuple(
        Configuration(
            p_electrons=_P_SPIN_ORBITALS,
            d_electrons=n_d_electrons + holes,
            ligand_holes=holes,
        )
        for holes in range(max_holes + 1)
    )


def xas_final_configurations(
    n_d_electrons: int,
    max_ligand_holes: int = 1,
) -> tuple[Configuration, ...]:
    """Return `2p5 3d^(n+1)`, `2p5 3d^(n+2)L`, ... sectors."""

    if not 0 <= n_d_electrons <= _D_SPIN_ORBITALS - 1:
        raise ValueError("n_d_electrons must satisfy 0 <= n <= 9 for L-edge XAS")
    _validate_hole_request(n_d_electrons, max_ligand_holes)
    max_holes = min(max_ligand_holes, _D_SPIN_ORBITALS - (n_d_electrons + 1))
    return tuple(
        Configuration(
            p_electrons=_P_SPIN_ORBITALS - 1,
            d_electrons=n_d_electrons + 1 + holes,
            ligand_holes=holes,
        )
        for holes in range(max_holes + 1)
    )


def charge_transfer_configurations(
    n_d_electrons: int,
    max_ligand_holes: int = 2,
) -> tuple[Configuration, ...]:
    """Return valence-only `3d^n + 3d^(n+1)L + ...` sectors."""

    _validate_hole_request(n_d_electrons, max_ligand_holes)
    max_holes = min(max_ligand_holes, _D_SPIN_ORBITALS - n_d_electrons)
    return tuple(
        Configuration(d_electrons=n_d_electrons + holes, ligand_holes=holes)
        for holes in range(max_holes + 1)
    )


def basis_from_configurations(configurations: tuple[Configuration, ...] | list[Configuration]) -> FockBasis:
    """Build a Fock basis from compatible `2p + 3d + ligand` configurations."""

    configs = tuple(configurations)
    if not configs:
        raise ValueError("at least one configuration is required")
    has_core_shell = configs[0].has_core_shell
    total_electrons = configs[0].n_electrons
    for config in configs:
        if config.has_core_shell != has_core_shell:
            raise ValueError("cannot mix configurations with and without a 2p shell")
        if config.n_electrons != total_electrons:
            raise ValueError("all configurations must have the same total electron count")

    states: list[int] = []
    for config in configs:
        states.extend(configuration_states(config))
    n_orbitals = _P_SPIN_ORBITALS + 2 * _D_SPIN_ORBITALS if has_core_shell else 2 * _D_SPIN_ORBITALS
    return FockBasis.from_states(n_orbitals, total_electrons, states)


def configuration_states(config: Configuration) -> tuple[int, ...]:
    """Enumerate bit states for one configuration sector."""

    d_basis = FockBasis.fixed_n(_D_SPIN_ORBITALS, config.d_electrons)
    ligand_basis = FockBasis.fixed_n(_D_SPIN_ORBITALS, config.ligand_electrons)
    states: list[int] = []
    if config.has_core_shell:
        assert config.p_electrons is not None
        p_basis = FockBasis.fixed_n(_P_SPIN_ORBITALS, config.p_electrons)
        for p_state in p_basis.states:
            for d_state in d_basis.states:
                for ligand_state in ligand_basis.states:
                    states.append(
                        p_state
                        | (d_state << _P_SPIN_ORBITALS)
                        | (ligand_state << (_P_SPIN_ORBITALS + _D_SPIN_ORBITALS))
                    )
    else:
        for d_state in d_basis.states:
            for ligand_state in ligand_basis.states:
                states.append(d_state | (ligand_state << _D_SPIN_ORBITALS))
    return tuple(states)


def configuration_energy(
    config: Configuration,
    delta: float,
    u_charge_transfer: float = 0.0,
    core_hole_potential: float = 0.0,
) -> float:
    """Simple diagonal sector energy used by the current XAS/XLD solver.

    Without a core hole this gives `0`, `Delta`, `2Delta + Udd`, ...
    for `L0`, `L1`, `L2`, ... sectors. With one 2p core hole this gives
    `0`, `Delta + Udd - Upd`, `2Delta + 3Udd - 2Upd`, ... .
    """

    holes = config.ligand_holes
    if config.core_holes:
        energy = holes * delta + holes * (holes + 1) * u_charge_transfer / 2.0
        energy -= config.core_holes * holes * core_hole_potential
    else:
        energy = holes * delta + holes * (holes - 1) * u_charge_transfer / 2.0
    return float(energy)


def _validate_hole_request(n_d_electrons: int, max_ligand_holes: int) -> None:
    if not 0 <= n_d_electrons <= _D_SPIN_ORBITALS:
        raise ValueError("n_d_electrons must satisfy 0 <= n <= 10")
    if max_ligand_holes < 0:
        raise ValueError("max_ligand_holes must be non-negative")

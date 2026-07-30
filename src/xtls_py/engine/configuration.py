from __future__ import annotations

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


def core_xps_final_configurations(
    n_d_electrons: int,
    max_ligand_holes: int = 1,
) -> tuple[Configuration, ...]:
    """Return `2p5 3d^n`, `2p5 3d^(n+1)L`, ... sectors for core-level XPS.

    The photoelectron leaves the cluster, so unlike XAS the d count does not
    grow: only the 2p shell loses one electron. These are the `#f` sectors of
    the XTLS `Mode=XPS` X-card.
    """

    _validate_hole_request(n_d_electrons, max_ligand_holes)
    max_holes = min(max_ligand_holes, _D_SPIN_ORBITALS - n_d_electrons)
    return tuple(
        Configuration(
            p_electrons=_P_SPIN_ORBITALS - 1,
            d_electrons=n_d_electrons + holes,
            ligand_holes=holes,
        )
        for holes in range(max_holes + 1)
    )


def valence_xps_final_configurations(
    n_d_electrons: int,
    max_ligand_holes: int = 1,
) -> tuple[Configuration, ...]:
    """Return `2p6 3d^(n-1)`, `2p6 3d^n L`, ... sectors for valence-band XPS.

    Here the photoelectron is removed from the valence shell, so the core shell
    stays filled and the d count drops by one.
    """

    if not 1 <= n_d_electrons <= _D_SPIN_ORBITALS:
        raise ValueError("n_d_electrons must satisfy 1 <= n <= 10 for valence XPS")
    _validate_hole_request(n_d_electrons, max_ligand_holes)
    max_holes = min(max_ligand_holes, _D_SPIN_ORBITALS - n_d_electrons)
    return tuple(
        Configuration(
            p_electrons=_P_SPIN_ORBITALS,
            d_electrons=n_d_electrons - 1 + holes,
            ligand_holes=holes,
        )
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


def sector_energy(
    ligand_holes: int,
    core_holes: int,
    delta: float,
    u_charge_transfer: float = 0.0,
    core_hole_potential: float = 0.0,
    d_electron_offset: int = 0,
) -> float:
    """Diagonal energy of one charge-transfer sector.

    This is the single definition of the sector energy used everywhere: by the
    diagonal of the many-body Hamiltonian, and by the configuration-energy
    table written to the output files.

    Following the configuration-centroid formula of the XTLS manual, writing
    `h` for the ligand-hole count, `c` for the core-hole count and `dd` for the
    d-electron offset of this final state relative to the initial `3d^n`
    reference,

        E = h*Delta + Udd * h * (2*dd + h - 1) / 2 - c * h * Udc.

    `d_electron_offset` is what distinguishes the spectroscopies, because
    `Delta` is defined on the initial configuration and each sector adds a
    different number of d-d pairs:

    =====================  ====  ====  ===========================
    sector                   dd     c  L1 energy
    =====================  ====  ====  ===========================
    initial `3d^n`            0     0  `Delta`
    XAS `2p5 3d^(n+1)`        1     1  `Delta + Udd - Udc`
    core XPS `2p5 3d^n`       0     1  `Delta - Udc`
    valence XPS `3d^(n-1)`   -1     0  `Delta - Udd`
    =====================  ====  ====  ===========================
    """

    holes = ligand_holes
    energy = holes * delta
    energy += u_charge_transfer * holes * (2 * d_electron_offset + holes - 1) / 2.0
    energy -= core_holes * holes * core_hole_potential
    return float(energy)


def _validate_hole_request(n_d_electrons: int, max_ligand_holes: int) -> None:
    if not 0 <= n_d_electrons <= _D_SPIN_ORBITALS:
        raise ValueError("n_d_electrons must satisfy 0 <= n <= 10")
    if max_ligand_holes < 0:
        raise ValueError("max_ligand_holes must be non-negative")

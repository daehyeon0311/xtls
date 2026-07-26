from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations


@dataclass(frozen=True)
class FockBasis:
    """Fixed-particle-number fermionic Fock basis.

    States are stored as integer bit patterns. Orbital index 0 corresponds to
    the least significant bit.
    """

    n_orbitals: int
    n_electrons: int
    states: tuple[int, ...]
    index: dict[int, int]

    @classmethod
    def fixed_n(cls, n_orbitals: int, n_electrons: int) -> "FockBasis":
        if n_orbitals < 0:
            raise ValueError("n_orbitals must be non-negative")
        if not 0 <= n_electrons <= n_orbitals:
            raise ValueError("n_electrons must satisfy 0 <= n_electrons <= n_orbitals")
        states = []
        for occupied in combinations(range(n_orbitals), n_electrons):
            state = 0
            for orbital in occupied:
                state |= 1 << orbital
            states.append(state)
        ordered = tuple(states)
        return cls(
            n_orbitals=n_orbitals,
            n_electrons=n_electrons,
            states=ordered,
            index={state: idx for idx, state in enumerate(ordered)},
        )

    @classmethod
    def from_states(
        cls,
        n_orbitals: int,
        n_electrons: int,
        states: list[int] | tuple[int, ...],
    ) -> "FockBasis":
        if n_orbitals < 0:
            raise ValueError("n_orbitals must be non-negative")
        if not 0 <= n_electrons <= n_orbitals:
            raise ValueError("n_electrons must satisfy 0 <= n_electrons <= n_orbitals")
        ordered = tuple(states)
        for state in ordered:
            if state < 0 or state >= (1 << n_orbitals):
                raise ValueError("state is outside the orbital space")
            if state.bit_count() != n_electrons:
                raise ValueError("all states must have n_electrons occupied orbitals")
        if len(set(ordered)) != len(ordered):
            raise ValueError("states must be unique")
        return cls(
            n_orbitals=n_orbitals,
            n_electrons=n_electrons,
            states=ordered,
            index={state: idx for idx, state in enumerate(ordered)},
        )

    def __len__(self) -> int:
        return len(self.states)

    def occupations(self, state: int) -> tuple[int, ...]:
        return tuple((state >> orbital) & 1 for orbital in range(self.n_orbitals))

    def label(self, state: int) -> str:
        return "".join(str(bit) for bit in reversed(self.occupations(state)))

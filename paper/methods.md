# 2. Methods

> Draft. Korean translation in `methods_ko.md`.

---

## 2.1 Charge-transfer cluster model

The Fe L₂,₃ spectra and the single-ion anisotropy are obtained from a
full-multiplet charge-transfer cluster model of the FeO₄ unit, in the same
framework as the XTLS code [Tanaka]. The basis spans 2p⁶3dⁿ, 2p⁶3dⁿ⁺¹L̲ and
2p⁶3dⁿ⁺²L̲² configurations, where L̲ denotes a ligand hole, and the many-body
Hamiltonian is diagonalized exactly within it.

The Hamiltonian contains the intra-atomic Coulomb and exchange integrals, the
3d and 2p spin–orbit couplings, the crystal field from the ligand positions,
and the d–ligand hybridization. Slater integrals and spin–orbit constants are
taken from the Appendix-A tables of Hartree–Fock values and reduced to 82.5%,
following common practice. Charge-transfer parameters are Δ = 5.5 eV,
U_dd = 4.25 eV and U_dc = 5.75 eV, as in Ref. [XLD].

The sector energies follow the configuration-centroid formula

    E(h) = h·Δ + U_dd · h · (2δ_d + h − 1)/2 − c · h · U_dc,

with h the ligand-hole count, c the core-hole count and δ_d the d-electron
offset of the sector relative to the initial 3dⁿ reference. Writing it this way
keeps the initial (δ_d = 0), XAS-final (δ_d = +1) and photoemission-final
(δ_d = 0 with c = 1) sectors on one expression.

### Crystal-field normalization

The point-charge field is scaled so that its t₂–e separation equals the
specified 10Dq. **The two groups are identified by their cubic symmetry labels**
— (xy, yz, zx) against (x²−y², 3z²−r²) — rather than by the sign of their
energies.

This matters under distortion. Sorting levels by sign is equivalent while the
field is close to cubic, but as the tetrahedron is compressed the yz, zx and
x²−y² levels cross zero, the grouping changes composition, and the
normalization jumps discontinuously. Across a scan from 4° to 12° the
sign-based rule produces steps of up to 0.140 meV between adjacent points,
against 0.009 meV for the symmetry-based rule (Fig. S1), and it also runs
systematically steeper — an artifact that inflates any distortion dependence
extracted from it.

Because 10Dq is a fitted quantity rather than an observable, its numerical
value depends on the convention it was fitted under. The two conventions differ
by a factor that itself varies with the distortion (2.04 at the Ba₂FeSi₂O₇
angle), so a quoted 10Dq is meaningful only together with its convention. The
value used here, 0.1508 eV, is calibrated to reproduce D = 1.45 meV for
Ba₂FeSi₂O₇ under the symmetry-based rule; it corresponds to 0.30 eV under the
sign-based one.

For a distortion scan the field is calibrated **once** at the reference
structure and that scale is then held fixed, allowing the point-charge model to
supply the absolute variation with the ligand positions. Renormalizing at every
geometry would discard precisely that variation and retain only the change in
shape.

### Extracting D

The five lowest eigenstates form the S = 2 manifold, separated from the next
multiplet by 338 meV at the experimental structure. They are fitted to

    E(S_z) = D S_z² + A S_z⁴,

solved on the doublet centroids as D = (16E₁ − E₂)/12 and A = (E₂ − 4E₁)/12,
where E₁ and E₂ are the S_z = ±1 and ±2 levels. For the family screening
(§3.4) the fit is generalized to any spin by regressing the 2S+1 lowest levels
onto the S_z² ladder.

---

## 2.2 Broken-symmetry density-functional theory

The exchange coupling is obtained with Quantum ESPRESSO 7.5 [QE] using the PBE
functional with a Hubbard correction on the Fe 3d shell, PAW datasets, plane-wave
cutoffs of 60 Ry for wavefunctions and 480 Ry for the density, a 4×4×6
Monkhorst–Pack mesh and Gaussian smearing of 0.01 Ry. Convergence is 10⁻⁸ Ry.

The cell is the 24-atom unit cell of Ba₂FeSi₂O₇ built from the Rietveld
refinement of Ref. [Jang] at 1.7 K, space group P-42₁m with a = 8.3194 Å and
c = 5.3336 Å. The resulting FeO₄ tetrahedron has Fe–O = 1.9875 Å with the four
oxygens exactly equivalent, and a deviation of 8.225° from the ideal
tetrahedral angle.

J follows from the energy difference between ferromagnetic and Néel alignments
of the two Fe sites, over the four nearest-neighbour bonds per cell:

    E_FM − E_AFM = 8 J S(S+1),   S = 2.

The magnetisation is a useful check: every run gives 8.00 μ_B per cell
ferromagnetic and 0.00 total with 7.99 absolute for Néel, as required for two
high-spin d⁶ sites.

The Hubbard parameter is fixed by experiment rather than convention. Scanning
U = 3, 4 and 5 eV gives J = 0.1552, 0.1188 and 0.0910 meV; linear interpolation
reproduces the measured 0.0887 meV at U = 5.07 eV, and U = 5 eV is used
throughout.

*[pending: interlayer J′ from a doubled cell along c, 48 atoms, in-plane order
held antiferromagnetic while the relative sign between layers is switched]*

---

## 2.3 Alignment of conventions

Comparison with Ref. [Do] requires care. That work quotes an effective S = 1
Hamiltonian obtained by projecting the S = 2 problem onto the S_z = 0, ±1
subspace, and states the mapping beneath its Eq. (2) as J̃ = 3J, Δ̃ = Δ/3 and
D̃ = D. Its fitted values J̃ = 0.266 meV and D̃ = 1.42 meV therefore correspond
to J = 0.0887 meV and D = 1.42 meV in the S = 2 Hamiltonian, and the critical
ratio α_c = J̃/D̃ = 0.158 becomes

    (D/J)_c = 3/α_c = 19.0.

All values in this work are quoted in the S = 2 convention. The distinction is
not cosmetic: comparing a classical S² exchange from broken-symmetry DFT
against an effective S = 1 ratio differs by a factor of order three and inverts
conclusions about which Hubbard parameter is admissible.

We further note that α_c is not a universal constant. Quantum Monte Carlo gives
0.18 in two dimensions and 0.10 in three [Ref], with 0.158 corresponding to the
interlayer ratio J′/J = 0.1 appropriate to Ba₂FeSi₂O₇. Since c-axis strain
changes the interlayer separation, α_c moves with the deformation; this is
addressed in §3.3.3.

---

## 2.4 Deformation geometry

Two deformations are distinguished throughout.

**Chemical (A-site substitution).** A straight line through the crystallographic
structures of Ba₂FeSi₂O₇ (Δθ = 8.23°, Fe–O = 1.988 Å) and Sr₂FeSi₂O₇
(Δθ = 4.86°, Fe–O = 1.969 Å), pinned at the former. Here the Fe–O distance
grows with the distortion, since a and c contract together under the smaller
A-site cation. This is the constraint Ref. [XLD] expresses as a fixed Poisson
ratio.

**Uniaxial (c-axis strain).** The c parameter is scaled with fractional
coordinates held fixed, at constant a. The in-plane projection of the Fe–O
vector is then fixed and the bond shortens as the distortion increases —
opposite to the chemical route.

Uniaxial strain is used for all deformed calculations, since it is what strain
control means experimentally. Scaling the cell deforms every polyhedron
together and keeps Si–O within its normal range (1.60–1.65 Å across
6°–10°). Displacing only the oxygens coordinating Fe, while leaving Si fixed,
stretches Si–O to 1.76 Å at 4° and compresses it to 1.49 Å at 12°, since O3
bridges the two cations; the self-consistent field does not converge on such
structures.

---

## 2.5 Code availability

The cluster-model implementation used here is an independent Python
reimplementation of the XTLS physics, and reproduces the original Fe L₂,₃
spectrum to 5.8% RMS (§3.1). It is available at [repository], together with the
input files for every calculation reported. Eigensolver starting vectors are
pinned so that runs are bit-reproducible; note that within a degenerate
subspace this fixes the basis without making any particular choice correct, so
quantities that are not subspace invariants — per-state moments, spin-resolved
channels — should be read with that in mind.

---

## Notes for revision

- [ ] Fill in §2.2 once the interlayer runs finish
- [ ] Reference tags: [Tanaka] Tanaka & Jo JPSJ 63, 2788 (1994); [QE] Giannozzi
      et al.; [Jang] PRB 104, 214434 (2021); [Do] Nat. Commun. 12, 5331 (2021);
      [XLD] Choi et al. JKPS 88, 933 (2026)
- [ ] Decide whether §2.1's normalization discussion belongs in the main text
      or Supplementary — it is a methodological contribution but long

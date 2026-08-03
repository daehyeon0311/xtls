# 3. Results

> Draft. Numbers are final except where marked *pending*.
> Korean translation in `results_ko.md`.
>
> On tone: statements about Ref. [XLD] are kept factual and specific — what was
> assumed, what was computed, what follows. The paper is not wrong in its
> measurements or in D; what changes is a prediction that rested on an
> uncomputed quantity, and saying so plainly is both accurate and the strongest
> position. No hedging, no editorialising.

---

## 3.1 Validation of the cluster model

We first establish that the present implementation reproduces the published
spectra and anisotropies before applying it to deformed structures.

Our charge-transfer cluster calculation reproduces the Fe L₂,₃ absorption
spectrum obtained with XTLS to 5.8% RMS over the L₃ and L₂ regions, once the
two codes' energy zeros are aligned (the offset is 4.93 eV and carries no
physical content). The Appendix-A Slater integrals used here agree with those
printed in the XTLS manual to three decimal places.

For Ba₂FeSi₂O₇ the five lowest states form a well-isolated S = 2 manifold, the
next multiplet lying 338 meV above it. Fitting those five levels to
E(S_z) = D S_z² + A S_z⁴ gives

    D = 1.450 meV,   A = +0.0017 meV,

in agreement with D = 1.42 meV from inelastic neutron scattering [Do]. The
quartic term is 0.1% of D and was not resolved previously.

We also find that the S_z = ±2 doublet is split by 0.0138 meV. Ref. [Mai] noted
that these are two distinct states of B symmetry, degenerate only to second
order in the spin–orbit coupling; the present calculation quantifies that
splitting.

**Sr₂FeSi₂O₇ as an independent test.** Calibrating only on Ba₂FeSi₂O₇ and
leaving Sr₂FeSi₂O₇ as a prediction, we obtain D = 1.318 meV. Terahertz
absorption on Sr₂FeSi₂O₇ [Mai] resolves modes at 1.0 and 1.4 THz which, from
their temperature dependence and relative intensities, correspond to the
E(±1) → B(±2) and A(0) → B(±2) transitions at 3D and 4D respectively. Those
place D between 1.379 and 1.447 meV. Our value deviates by −4.4% and −9.0%
from the two assignments, against −11.5% and −15.7% for the 1.22 meV reported
in Ref. [XLD] (Table S1).

The lowest mode at 0.6 THz is reproduced by neither value, consistent with the
statement in Ref. [Mai] that its strong temperature dependence falls outside a
single-ion description.

**Prediction for Ba₂FeSi₂O₇.** No terahertz measurement has been reported for
Ba₂FeSi₂O₇. With D = 1.45 meV the corresponding modes lie at 1.052 and
1.402 THz, about 5% above those of Sr₂FeSi₂O₇.

---

## 3.2 Exchange coupling from first principles

Ref. [XLD] computed D and took J from neutron scattering. We compute both.

Broken-symmetry calculations on the 24-atom cell give the Néel arrangement
lower in energy than the ferromagnetic one at every Hubbard U examined,
consistent with the ordering observed below 5.2 K. The magnetisation is exact
in each run: 8.00 μ_B per cell ferromagnetic, and 0.00 total with 7.99 absolute
for Néel, i.e. two S = 2 sites.

Comparison with experiment requires care over conventions. Ref. [Do] quotes an
effective S = 1 Hamiltonian obtained by projecting out the S_z = ±2 doublet and
states the mapping as J̃ = 3J, D̃ = D. Converted to the S = 2 Hamiltonian, their
values are J = 0.0887 meV and D = 1.42 meV, and the critical ratio
α_c = J̃/D̃ = 0.158 becomes (D/J)_c = 3/α_c = 19.0. All values below are in the
S = 2 convention, with J extracted from E_FM − E_AFM = 8 J S(S+1).

| U (eV) | J (meV) | deviation from [Do] |
|---|---|---|
| 3 | 0.1552 | +75% |
| 4 | 0.1188 | +34% |
| 5 | 0.0910 | +2.6% |

U = 5 eV reproduces the measured exchange to 2.6%, and linear interpolation
places exact agreement at U = 5.07 eV. This lies in the range ordinarily
adopted for Fe²⁺, and we use U = 5 eV throughout.

Combining with the cluster-model anisotropy,

    D/J = 1.4445 / 0.0910 = 15.9,

against 16.0 from experiment. Two independent calculations — a charge-transfer
cluster model for D, broken-symmetry DFT for J — therefore reproduce the
measured ratio and place Ba₂FeSi₂O₇ at 84% of the critical value.

---

## 3.3 Distortion dependence of D and J

This section contains the central result.

### 3.3.1 Two deformations, two D(Δθ)

"Distortion" of the FeO₄ tetrahedron can be varied in two physically distinct
ways, and they are not equivalent (Fig. 2).

Interpolating between the crystallographic structures of Ba₂FeSi₂O₇
(Δθ = 8.23°, Fe–O = 1.988 Å) and Sr₂FeSi₂O₇ (Δθ = 4.86°, Fe–O = 1.969 Å)
describes A-site substitution: the Fe–O distance *grows* with the distortion,
because the a axis contracts along with c. Uniaxial strain along c at fixed a
does the opposite — the in-plane projection of the Fe–O vector is fixed, so the
bond *shortens* as the distortion increases.

The resulting anisotropies differ accordingly. Along the substitution trend D
rises monotonically by 33% between 4° and 12°. Under uniaxial strain it traces
a shallow minimum near 8° and varies by only 2.5% over the same range.

This distinction matters because Ref. [XLD] evaluated D(Δθ) along the
substitution trend while proposing strain as the experimental route. The
mechanism it invokes — increasing compression raises D — holds for the geometry
in which D was computed, but not for the one proposed.

We adopt uniaxial c-axis strain throughout, since that is what the proposal
means. Fractional coordinates are held and c is scaled, so every polyhedron
deforms together; displacing only the oxygens breaks the Si–O bonds they also
belong to (1.62 → 1.76 Å at 4°) and the calculation does not converge.

### 3.3.2 J dominates the distortion dependence

Table 3 gives both quantities on identical geometries.

| Δθ | c scale | Fe–O (Å) | J (meV) | D (meV) | D/J |
|---|---|---|---|---|---|
| 6.00° | 1.0979 | 2.0292 | 0.0662 | 1.4730 | 22.25 |
| 8.225° | 1.0000 | 1.9875 | 0.0910 | 1.4445 | 15.88 |
| 10.00° | 0.9247 | 1.9575 | 0.1111 | 1.4554 | 13.10 |

Over this range **J changes by +67.8% while D changes by −2.5%** (Fig. 3a). The
quantity Ref. [XLD] assumed to be negligible is the dominant one.

The mechanism is direct: compressing c shortens Fe–O from 2.029 to 1.958 Å,
strengthening the Fe–O–O–Fe superexchange. D, by contrast, barely responds to
this particular deformation. Ref. [XLD] argued qualitatively that the two
effects on J — reduced hybridization from the longer bond, increased hopping
from the changed angle — would cancel. Our calculation shows they do not.

### 3.3.3 The phase boundary and its direction

Because J rises faster than D, **D/J decreases with compression** (Fig. 3b).
The boundary at (D/J)_c = 19.0 is crossed at Δθ = 6.95°, with the quantum
paramagnet on the *low*-distortion side. Ba₂FeSi₂O₇ at 8.23° lies on the
ordered side, as observed.

Reaching the transition therefore requires **tensile** strain along c, of about
+4.7% (c: 5.334 → 5.587 Å), and not the compression proposed previously.

Two caveats sharpen rather than weaken this. First, α_c itself depends on
dimensionality: quantum Monte Carlo gives 0.18 in two dimensions and 0.10 in
three [Ref], and the 0.158 used here corresponds to J′/J = 0.1. Tensile strain
separates the layers and pushes the system toward the two-dimensional limit,
which *lowers* (D/J)_c from 19.0 toward 16.7 while D/J is rising. Both move
toward the transition, so +4.7% is an upper bound; the two-dimensional limit
would require only +1.2%. *[pending: J′(strain) will replace this range with a
value]*

Second, the same reasoning excludes the bulk route. A-site substitution applies
chemical pressure but moves Fe–O the other way, and indeed Sr₂FeSi₂O₇ is *less*
distorted than Ba₂FeSi₂O₇, not more. Uniaxial strain is not merely convenient
here; it is the only available route.

Among cubic substrates, MgAl₂O₄(001) matches the 8.3194 Å in-plane lattice at
−2.84% misfit, implying +2.44% tensile strain along c for a Poisson ratio of
0.30 — mid-range for the target. SrTiO₃ reaches the conservative target but at
−6.1% misfit, which tends to relax. We are not aware of epitaxial growth of any
melilite, so the synthetic difficulty is untested.

---

## 3.4 Which ion can host the transition

A quantum-paramagnetic ground state must be a non-magnetic singlet. Kramers
degeneracy keeps every level of a half-integer-spin system at least doubly
degenerate, so no value of D produces one. This restricts the transition to
integer-spin members of the family regardless of how large the anisotropy
becomes (Fig. 4).

Substituting each divalent ion into the Ba₂FeSi₂O₇ geometry, with charge-transfer
parameters held common so that only the d count varies:

| ion | d | S | term | D (meV) | QPM possible |
|---|---|---|---|---|---|
| Mn²⁺ | 5 | 5/2 | ⁶A₁ | 0.0018 | no (Kramers) |
| **Fe²⁺** | 6 | 2 | ⁵E | 1.539 | **yes** |
| Co²⁺ | 7 | 3/2 | ⁴A₂ | 3.496 | no (Kramers) |
| Ni²⁺ | 8 | 1 | ³T₁ | not defined | see below |

Mn²⁺ serves as a check on the method: the half-filled shell is an orbital
singlet, the second-order spin–orbit route is closed, and D vanishes to
0.002 meV with the manifold isolated by 2.5 eV.

Co²⁺ has the largest anisotropy, 2.3 times that of Fe²⁺, and much the strongest
distortion dependence (+594% between 4° and 12°). The enhancement factor
separates into 1.61 from λ² (ζ₃d = 0.066 against 0.052 eV) and 1.41 from the
smaller ⁴A₂–⁴T₂ separation, reproducing the computed 2.27. This is consistent
with the large anisotropy of Ba₂CoGe₂O₇, but being a Kramers ion it orders
rather than becoming paramagnetic — as it does experimentally, despite
Λ/J ≈ 8 exceeding the ratio at which Ba₂FeSi₂O₇ sits.

Ni²⁺ is the only other integer-spin candidate, but ³T₁ retains its orbital
degeneracy: the S = 1 manifold is not isolated, and fitting the three lowest
levels to D S_z² leaves a residual of 70% of D. A single-ion anisotropy is not
defined for this site.

**Fe²⁺ is therefore effectively the only site in the melilite family where this
transition can occur**, which justifies the choice of material in Ref. [XLD]
after the fact.

---

## Notes for revision

- [ ] §3.1: decide whether the S_z = ±2 splitting belongs here or in Discussion
- [ ] §3.3.3: replace the α_c range with the computed J′ once available
- [ ] Ref. tags to be resolved: [XLD] = Choi et al. JKPS 88, 933 (2026);
      [Do] = Do et al. Nat. Commun. 12, 5331 (2021);
      [Mai] = Mai et al. PRB 94, 224416 (2016)
- [ ] Confirm the substrate misfit convention with the growth group before
      committing to MgAl₂O₄

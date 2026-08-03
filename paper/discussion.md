# 4. Discussion

> Draft. Korean translation in `discussion_ko.md`.

---

## 4.1 Why the exchange dominates

The result that decides everything else is the disparity between how J and D
respond to the same deformation: +68% against −2.5% over four degrees of
distortion. The asymmetry has a straightforward origin.

Compressing c at fixed a shortens the Fe–O bond from 2.029 to 1.958 Å. The
intralayer exchange proceeds through an Fe–O–O–Fe path, so it inherits the
strong distance dependence of the two Fe–O hybridizations at either end, and J
rises accordingly. The single-ion anisotropy, by contrast, arises from
second-order spin–orbit mixing between the ⁵A ground state and the ⁵E excited
state, D ~ λ²/ΔE. Under this particular deformation the two competing
contributions to ΔE — the tetragonal splitting grows while the overall crystal
field weakens with the longer in-plane projection — very nearly cancel, leaving
D almost unchanged and tracing a shallow minimum near the experimental
structure.

Ref. [XLD] anticipated a cancellation, but placed it in J rather than in D:
"the Heisenberg exchange J decreases by the former effect but increases by the
latter effect, hence the variation in J with the distortion can be minimal."
The present calculation finds the cancellation is in D, and that J varies
strongly. The two effects invoked for J — reduced hybridization from a longer
bond, enhanced hopping from the changed angle — do not in fact offset, because
under uniaxial strain the bond does not lengthen at all; it shortens, and both
effects then push J the same way.

That the observed T_N values of Ba₂FeSi₂O₇ (5.2 K) and Sr₂FeSi₂O₇ (5.0 K) are
so similar was cited as evidence for a distortion-insensitive J. It is not
evidence for the uniaxial case: those two compounds differ by A-site
substitution, along which Fe–O moves the other way, and their similar T_N is
consistent with the compensation that the chemical route happens to provide.

## 4.2 Geometry is part of the model

The broader lesson concerns what "distortion" denotes. A single angle Δθ does
not specify a structure. Two deformations reaching the same Δθ — one by
chemical substitution, one by uniaxial strain — differ in the Fe–O distance and
therefore in every quantity derived from it. Here they differ in the *sign* of
dD/dΔθ, not merely its size.

This is why the prediction reversed. Ref. [XLD] computed D(Δθ) along the
substitution trend, which was the natural choice given that its two compounds
define exactly that line, and then proposed strain as the route to larger Δθ.
Each step is reasonable in isolation; taken together they mix two geometries.
The result is that the published D(Δθ) belongs to a deformation that cannot be
applied experimentally, since one cannot continuously substitute the A-site
beyond the Ba end-member, while the deformation that can be applied gives a
different D(Δθ) and, once J is included, the opposite conclusion.

We would suggest that reports of anisotropy-versus-distortion state the
deformation path explicitly, as one would state a strain tensor rather than a
scalar.

## 4.3 What an experiment would need

The transition requires tensile c-axis strain of at most +4.7%, and less once
the change in α_c is accounted for. Among cubic substrates MgAl₂O₄(001) matches
the 8.3194 Å in-plane lattice at −2.84% misfit, implying +2.44% along c for a
Poisson ratio of 0.30 — comfortably inside the window and small enough to hold
below the critical thickness. SrTiO₃ would reach the conservative end of the
target but at −6.1% misfit, where relaxation through misfit dislocations is
likely.

We are not aware of epitaxial growth of any melilite, so the synthesis is
untested. The structure is layered with Ba slabs separating FeSi₂O₇ sheets,
which is not obviously unfavourable for layer-by-layer growth, but this remains
a genuine open question rather than a formality.

The signature to look for is the collapse of the ordered moment and the closing
of the longitudinal-mode gap at the zone centre, both already characterized in
the bulk material [Do]. Terahertz absorption would be the more accessible
probe for a thin film, and §3.1 provides the reference frequencies: 1.052 and
1.402 THz for unstrained Ba₂FeSi₂O₇, moving down as the transition is
approached.

A bulk route is not available. Chemical pressure through A-site substitution
moves Fe–O in the wrong direction — Sr₂FeSi₂O₇ is *less* distorted than
Ba₂FeSi₂O₇, not more — so no member of the substitution series approaches the
boundary. Uniaxial strain is not merely convenient here but necessary, which
raises the value of settling the epitaxy question.

## 4.4 The family, and why Fe²⁺

Kramers degeneracy is a hard constraint: a half-integer-spin site cannot have a
non-magnetic singlet ground state, so Mn²⁺ (S = 5/2) and Co²⁺ (S = 3/2) are
excluded however large their anisotropy. Ba₂CoGe₂O₇ illustrates this directly —
its Λ/J ≈ 8 exceeds the corresponding ratio in Ba₂FeSi₂O₇, yet it orders
magnetically, and its large anisotropy manifests instead as magnetoelectric
coupling and spin-stretching modes.

That leaves the integer-spin members. Ni²⁺ (S = 1) is nominally eligible but is
a ³T₁ orbital triplet in tetrahedral coordination: the S = 1 manifold is not
isolated and no single-ion anisotropy is defined, as the 70% fit residual
shows. Whether a Jahn–Teller distortion quenches the orbital degeneracy
sufficiently to restore a spin description is a separate question we have not
addressed, and one that would need the distorted structure rather than the
Ba₂FeSi₂O₇ geometry used here.

Fe²⁺ therefore stands alone in this family, which retrospectively justifies the
choice of material — a point worth making since the reasoning was not available
when the choice was made.

## 4.5 Limitations

**Two methods.** D comes from a cluster model and J from periodic DFT. The
conventions are aligned (§2.3) and both reproduce their measured counterparts,
but the two are not derived from one calculation, and a residual inconsistency
cannot be excluded. The agreement of D/J with experiment to 1% is reassuring
rather than conclusive.

**The Hubbard parameter.** J varies by a factor of 1.7 across U = 3 to 5 eV. We
fix U by requiring the measured exchange, which lands at 5.07 eV, well inside
the usual range for Fe²⁺; but this makes J a calibrated rather than a predicted
quantity. What is predicted is its *variation* with strain, which is what the
argument uses, and which is far less sensitive to U than the absolute value.

**Uniform c-axis scaling.** Internal coordinates were held rather than relaxed.
Under real epitaxial strain the ions would relax within the strained cell,
which would reduce the effective distortion at a given lattice strain and so
increase the strain required. The direction of the conclusion is unaffected,
but +4.7% should be read as a lower bound on the lattice strain for a given
target distortion.

**The critical ratio.** α_c is taken from quantum Monte Carlo for a model
Hamiltonian and depends on the interlayer coupling, which strain also changes
(§3.3.3). We have bounded the effect and shown it works in favour of the
conclusion; the calculation that would remove the bound is in progress.
*[pending: interlayer J′]*

---

## Notes for revision

- [ ] §4.3: soften or strengthen the epitaxy paragraph after talking to a
      growth group
- [ ] §4.5: fold in J′ when available and remove the pending marker
- [ ] Consider moving §4.2 to the Conclusion if the paper runs long — it is the
      most transferable point but not the headline result

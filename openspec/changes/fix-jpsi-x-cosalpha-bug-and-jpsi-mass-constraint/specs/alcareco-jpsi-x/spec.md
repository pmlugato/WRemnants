## MODIFIED Requirements

### Requirement: B-Level Kalman Vertex Fit Under Preset C

The B-level vertex fit applied to non-V0 candidates under preset C SHALL operate on physically meaningful quantities. Specifically, the cosα cut SHALL use the matched primary vertex (PV) — not the beamspot — as the origin of the displacement direction. The Lxy/σ cut SHALL be enforced unconditionally when active, with no silent bypass for `lxy < 1e-10 cm`. The PV used for both Lxy and cosα SHALL be the closest-in-z PV from the offline PV collection, with a fallback to `pvs->front()` only when the closest PV is more than 1 cm away in z.

The cosα variant SHALL be the **2D transverse** cosα: the angle between the B candidate's transverse momentum (px, py) and the transverse displacement vector from the matched PV to the B vertex (dvx, dvy). The 3D form is explicitly NOT used, because residual PV z-resolution propagates into the 3D angle in a way that degrades the cut's discrimination power.

Under preset C, the producer SHALL use `KinematicConstrainedVertexFitter` (not the previous `KalmanVertexFitter`) to perform both the B-vertex fit AND the J/ψ mass constraint in one call. The cut THRESHOLDS (`minBVtxProb`, `maxMotherAlphaBS`, `minBLxyOverSigma`) are unchanged from `add-jpsi-x-vertex-fit-and-low-pt`; only the underlying computations change.

#### Scenario: cosα origin is the matched PV, not the beamspot
- **WHEN** the producer evaluates the cosα cut under preset C
- **THEN** the displacement vector `(dvx, dvy)` is computed as `(vx_B − vx_PV, vy_B − vy_PV)` using the matched PV's position, AND the cosα is the 2D transverse cosα `(px_B · dvx + py_B · dvy) / (pT_B · sqrt(dvx² + dvy²))`, AND no z-component (`dvz`, `pz`) enters the calculation

#### Scenario: PV matching uses closest-in-z
- **WHEN** the producer selects a PV for Lxy and cosα computations
- **THEN** the PV chosen SHALL be `argmin_pv |vz_B − vz_pv|` over the offline PV collection. If the minimum exceeds 1 cm, the producer SHALL fall back to `pvs->front()` and emit a warning-level log message

#### Scenario: Lxy/σ cut applies unconditionally
- **WHEN** the producer evaluates the Lxy/σ cut under preset C, and the fitted B vertex coincides with the matched PV to within machine precision (`lxy < 1e-10 cm`)
- **THEN** the candidate SHALL be rejected (the cut MUST fail) — a vertex coincident with the PV is by definition not displaced and cannot satisfy `Lxy/σ > minBLxyOverSigma`

#### Scenario: Preset C uses the constrained vertex fitter
- **WHEN** the producer runs under preset C and `applyBVtxFit_ = true`
- **THEN** the B-level fit is performed by `KinematicConstrainedVertexFitter` with `MassKinematicConstraint(3.0969 GeV)` applied to the muon pair; the fit returns the constrained B 4-momentum and vertex in one call. The previous `KalmanVertexFitter` invocation is removed

### Requirement: B-Candidate Mass Computation

The B candidate 4-momentum SHALL be computed with the J/ψ mass constraint `m(μμ) = 3.0969 GeV` applied to the muon pair. The constraint is a property of the producer's mass calculation, not a property of any specific preset, and SHALL be applied to candidates emitted under preset B and preset C alike. Preset A's behavior is not affected by this requirement because preset A is no longer actively tested.

Two code paths are used depending on whether the B-level vertex fit runs:

- Under preset B (`applyBVtxFit_ = false`), the constraint is applied via a dimuon-only `KinematicParticleFitter` + `MassKinematicConstraint`; the constrained J/ψ 4-momentum is summed with the bachelor 4-momentum to give `lvM`; the B vertex remains the existing midpoint of input vertices.
- Under preset C (`applyBVtxFit_ = true`), the constraint is applied as part of the multi-track `KinematicConstrainedVertexFitter`, which returns both the constrained 4-momentum AND the B vertex.

On fit failure under either path (expected ~1 % of valid candidates), the producer SHALL fall back to the unconstrained `jpsi.p4() + bachelor.p4()` (track mode) or `jpsi.p4() + xCand.p4()` (VCC mode) sum AND the midpoint-of-input-vertices vertex, log the fallback at debug level, and NOT drop the candidate. Fallback frequency SHALL be exposed via the producer's summary log.

#### Scenario: B+ mass under preset B uses dimuon-only J/ψ constraint
- **WHEN** `produceTrackMode()` computes the B+ candidate 4-momentum under preset B for a passing (J/ψ, bachelor) pair
- **THEN** the producer invokes a dimuon-only `KinematicParticleFitter` with `MassKinematicConstraint(3.0969 GeV)` on the two muon TransientTracks, retrieves the constrained J/ψ 4-momentum, and computes `lvM = jpsi_constrained.p4() + lvBach`. The B vertex stays as the existing midpoint of input vertices. No B-level vertex fit is performed

#### Scenario: B+ mass under preset C uses the constrained vertex fitter
- **WHEN** `produceTrackMode()` computes the B+ candidate 4-momentum under preset C for a passing (J/ψ, bachelor) pair
- **THEN** the producer invokes `KinematicConstrainedVertexFitter` with the two muon TransientTracks + the bachelor TransientTrack and `MassKinematicConstraint(3.0969 GeV)` applied to the muon pair; the resulting fit's mother-particle 4-momentum and decay vertex become `lvM` and `vtxM`. The existing three cuts (vtxProb, cosα via the new PV-based 2D helper, Lxy/σ) are applied on the fit output

#### Scenario: VCC mode J/ψ constraint
- **WHEN** `produceVccMode()` computes a VCC-mode B candidate 4-momentum (e.g., B0 → J/ψ K*0)
- **THEN** under preset B the dimuon-only constraint applies to the muon pair before summing with the X-side `xCand.p4()`; under preset C the multi-track `KinematicConstrainedVertexFitter` includes all leaf tracks (μ⁺, μ⁻, plus the 2 X-side daughters) with the J/ψ mass constraint on the (μ⁺, μ⁻) pair

#### Scenario: Fit-failure fallback preserves candidates
- **WHEN** either the dimuon-only `KinematicParticleFitter` (preset B) or the multi-track `KinematicConstrainedVertexFitter` (preset C) fails to converge or returns an invalid result
- **THEN** the producer SHALL fall back to the unconstrained 4-momentum sum AND the midpoint-of-daughter-vertices vertex, log the fallback at debug level, and NOT drop the candidate. Fallback counters SHALL be exposed via the producer's summary log

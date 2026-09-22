# Participant-sign confidence bounds for population mean rank

`multivar.inference.TaskMeanRank` implements two separately named procedures:
`FrobeniusFeasibility` uses per-transformation low-rank feasibility at every
rank, and `OperatorFirstThenBound` uses the direct operator-norm first-root
test followed by the conservative operator bound. Neither procedure estimates
and removes a nuisance component. Do not select the smaller observed p-value
across procedures.

The target is the rank of stacked positive-weighted population group means
after fixed row maps. Participants are independent units. Each participant's
complete profile must be jointly centrally symmetric about its group center;
featurewise marginal symmetry is insufficient. Finite first moments are needed
to interpret centers as population means. Participant distributions and group
sizes may differ. The API records coordinate and source declarations but does
not establish these assumptions or authorize a data-dependent coordinate map.

`TaskRankGroup` retains participant identities and raw finite matrices, one
fixed `MatrixEnclosure` row map, and a positive `RealInterval` weight. The
intervals must contain the intended fixed scientific map and weight. Gale
arithmetic encloses every map product, mean, weighted sign sum, centered sum,
Frobenius feasibility matrix, comparison radius and norm. Constant group sign
patterns are recognized by integer counts. The active-block value
`h=4(nPlus/n)(nMinus/n)` encloses the exact rational without cancellation.

A candidate SVD is validated through Gale's factor residual and orthogonality
certificate. Tail energies are sums of enclosed tail singular values squared;
leading energies are never subtracted from a total. A draw is excluded only
when interval bounds prove the required strict inequality. An ambiguous
comparison counts. Arithmetic/certificate failures refuse the entire result.
A per-draw monotone envelope preserves numerical rank monotonicity, recording
any additional count as `ConservativeCarry`.

`TaskRankPlan.exhaustive` constructs every full-group sign pattern within an
explicit materialization budget. Its p-values are counts divided by group
size. `fixedUniform` generates a fixed budget through Resample4s integer
SplitMix64 draws; seed and budget must be set independently of response data.
Its numerator and denominator both include the extra identity. `replayFixed`
retains supplied signs and unique replicate identities; its validity is
conditional on the caller's fixed-budget independent uniform full-group
sampling declaration. Repeated random sign patterns are permitted. An
arbitrary supplied array cannot establish its own sampling law.

The significance level is an exact positive rational below one. Decisions
compare integer count products, not rounded display p-values. Testing starts
at rank zero and stops at the first non-rejection. The returned rank is a lower
confidence bound, not an exact dimension estimate or a significance claim for
individual axes. No division of alpha across ordered roots is needed.
The finite Monte Carlo value is not a guaranteed numerical upper bound on the
exhaustive p-value. Adaptive budgets or seed searching require another proof.

All transformations are evaluated before a receipt is returned. Cancellation
before preparation, between draws, or before delivery returns no partial
inferential evidence. A receipt retains both the tested prefix and all
computational comparisons, counts, identities, declared coordinates, maps,
weights, sampling convention, method version and exact significance level.
An application must additionally bind the actual provider closure and its
input-model qualification. This candidate establishes neither realistic power
nor a product default; those are distinct downstream qualification gates.

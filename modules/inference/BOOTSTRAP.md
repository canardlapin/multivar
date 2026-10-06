# Bootstrap means of complete matrix profiles

`MeanProfileBootstrap` refits the singular value decomposition of the mean of
independent, equally weighted matrix profiles. Each sampling unit contributes
the same `d × P` matrix. Store those matrices as consecutive row blocks in one
`(N d) × P` input; a bootstrap count applies to the complete block.

This is a general inference primitive for fixed matrix profiles. Applications
define the scientific contrasts, sampling-unit identities, reference weights and
feature geometry before calling it. It does not resample individual profile rows,
propagate uncertainty from upstream estimation, or certify confidence coverage.

## Fit in a fixed coordinate system

```scala
import gale.linalg.Matrix
import multivar.inference.*

val profiles = Matrix(2, 3)(1, 1, 0, 1, -1, 0)
val core = MeanProfileBootstrap.prepare(
  profiles,
  coordinatesPerUnit = 1,
  options = ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced)
).toOption.get

val draw = ProfileBootstrapCounts.fromIndices(2, Vector(0, 0)).toOption.get
val fit = core.refit(draw).toOption.get
val featureBlock = fit.loadingBlock(0, 3).toOption.get
```

The ordinary resampling interface composes with the core:
`BootstrapAction.rows(RowCount(core.unitCount).toOption.get)` draws the sampling
units, and `ProfileBootstrapCounts.fromDraw` turns its `Draw` into multiplicities.
Use explicit indices for independent comparisons. The order of identical copies
does not affect this mean estimator; indices are never interpreted as independent
within-profile observations.

`Reduced` prepares an untruncated Householder QR of the transpose of the stacked
profiles: `Z' = Q T`, then contracts the counts into `L_b T'`. The per-draw SVD
has only `d × min(N d, P)` entries. A one-row mean is normalized directly without
an eigensolver or SVD. Loadings are lifted through Q only when requested.

Gale owns QR and the full bidiagonal SVD. The thin Q is formed by applying the
compact QR to the first coordinate columns, avoiding `QR.q` and its `P × P`
allocation. QR rank diagnostics do not truncate preparation. The observed group
rank cannot restrict bootstrap directions. Every replicate applies the declared
absolute-zero and relative singular-value cutoff after aggregation.

`Direct` refits the same estimator on the original feature axis. `Auto` uses a
conservative work estimate, expected replicate count and working-byte limit to
choose a route; it does not benchmark the machine during analysis. Its estimate
includes input and factorization buffers, not arbitrary simultaneous application
studies or unrequested downstream outputs. Callers must budget those separately.

Multiplicity contraction sums complete copies with compensated addition and
divides once afterward. Fractional per-unit weights would introduce product
rounding that compensation cannot recover when large values cancel. Unrepresentable
accumulation or scaling fails explicitly; no underflowed nonzero source is
silently admitted as an exact zero. Profile entries are scaled before QR to avoid
overflow. Reconstruction and
orthogonality are checked. When cancellation exposes a mean below the prepared
error allowance, the core aggregates the original entries and verifies their
representation at the mean's own scale; the fit records `QrDirectRefinement`.
It refuses an unsupported residual or a singular value within numerical
uncertainty of the rank cutoff. The current implementation does not use the
unguarded Gram shortcut. Exact means the same estimator in real arithmetic,
with qualified floating-point error, not bitwise equivalence across backends.

## Map moments and linear readouts

`fit.loadingCoordinates` carries directions in the core's coordinate system.
`fit.coordinates(scaled = true)` multiplies each direction by its singular value.
The fit's `leftDirections`, `singularValues` and loadings reconstruct its mean.
Rank zero is a valid empty fit; it is not filled with fabricated directions.

Use the existing `ComponentAlignment` on reduced coordinates when the scientific
contract admits isolated axes. Apply its permutation/sign to all associated
curves and singular values where relevant. For tied components, use
`PrincipalAngles` on the reduced bases and `SubspaceStabilityReducer`, with
explicit rank-loss/failure accounting. A rotation does not make tied axes
identifiable, and bootstrap stability does not test their significance.

`ReducedLoadingMoments` accumulates an isolated axis's full coordinate covariance
with Welford updates. `addBatch` uses one local accumulator allocation per batch;
`combine` merges immutable states with the between-batch mean correction.
Keys and coordinate dimensions must match. Fewer than two valid replicates
return unavailable covariance. Summaries expose:

```text
mean(v) = Q mean(a)
Var(v_j) = Q[j,:] Cov(a) Q[j,:]'
```

Call `summary.projectRows(core.basisBlock(from, until))` with the unwrapped
`Either` result to recover means, sample variances and standard deviations in
feature blocks. Retaining only the diagonal of `Cov(a)` is incorrect. There is
no feature-by-feature covariance allocation. Full reduced covariance still costs
`O(r²)`; the constructor checks its working budget. A direct fit with a very
large feature axis can instead use ordinary streamed map moments, or the caller
can select the reduced route for that output workload.

Moment summaries do not produce percentile intervals. Exact quantiles require
retaining/replaying draws and lifting each requested feature block before
taking quantiles. Keep scaled and unscaled salience conventions consistent
when computing any bootstrap ratio; a zero SE is not voxelwise significance.

For repeated readouts, `core.project(readout)` caches `readout * Q` once, and its
`scores(fit)` multiplies only small matrices. It rejects a fit from a different
core. `featureFrom` supports a contiguous feature subset, such as a time slice;
the application retains the subset's physical identity. Original-observation
scores project the original observations through each fitted direction. Means
of the resampled observations are a different quantity and can be constructed
by applying the same counts to cached readout rows.

Prepared cores and moment states are immutable and can be shared across worker
chunks. The caller controls scheduling and derives draws by replicate identity
through the existing resampling API. No thread pool or platform runtime lives in
this module. Cancellation callbacks are checked between stages and at bounded
profile/block/batch boundaries; an individual Gale factorization is synchronous.

## Evidence

`MeanProfileBootstrapSuite` contains 416 explicit draws from independently
computed NumPy raw-cell marginalizations and full SVDs. It compares complete
fits, map moments, original-observation readouts and time-slice contractions.
Additional cases cover cancellation, omitted observed directions, zero/tied
fits, rank loss, rank boundaries, invalid shapes/budgets and covariance merging.
Shared tests run on JVM and Scala.js. `MeanProfileBootstrapSmoke` uses published
artifacts, and the opt-in JVM `MeanProfileBootstrapBenchmark` compares direct and
reduced execution with identical profiles, counts and map-moment outputs.

The scientific precedent for a fixed bootstrap span is
[Fisher et al. (2016)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5014451/).
This matrix-profile mean contract and its downstream interpretations have their
own numerical and statistical admission requirements.

## Complete-profile readout means

`ProfileBootstrapCounts.meanProfiles(values, coordinatesPerUnit)` exposes the
same compensated complete-copy aggregation used by the refitter. Rows are
unit-major, and all units remain present even when their multiplicity is zero.
The method validates shape, finite inputs, output-memory admission and
cancellation. It produces a mean without an additional decomposition; downstream
adapters can reuse the exact draw for score and response-context readouts.

## Delete-one refits and studentized scalar intervals

`counts.deleteOne(unit)` creates a distinct `ProfileDeletionCounts`: one copy
of a present unit is removed from an N-of-N bootstrap multiset. Its domain still
contains all N original units, but its sample size and mean denominator are
N−1. The ordinary `ProfileBootstrapCounts` constructor still requires N draws.
`core.refitDeletion(deletion)` reuses the prepared basis and returns a
`MeanProfileDeletionFit`; loading blocks and cached readout projections accept
that fit with the same owner checks as bootstrap fits. No retained direction is
assumed to remain identifiable after deletion.

`JackknifeStandardError.compute(values, expectedDeletions)` computes
`sqrt((N−1)/N * sum((value_i − mean(value))²))` from exactly N finite deletion
values. When a bootstrap multiset repeats a unit, deleting each copy gives the
same fit, but that deletion value must occur with its original multiplicity in
the N-position vector. A constant vector returns zero; the interval constructor
refuses zero standard errors. Stable scaling avoids unnecessary overflow and
preserves small differences around a large common location.

`StudentizedValue` contains an estimate and a strictly positive **standard
error**, not a variance. `StudentizedInterval.compute` consumes the complete,
ordered stream of `StudentizedDraw`s for the supplied replicate identities. It
forms `(bootstrapEstimate − observedEstimate) / bootstrapSE`, computes the two
R type-7 pivot quantiles, and inverts them using the observed standard error.
Missing, extra, reordered, repeated or failed draws cannot silently become an
available-case interval. The caller supplies any stricter SE floor, confidence
level, identity set, cancellation callback and working-memory budget. The
reducer retains one pivot per draw; it does not retain maps or allocate a
feature-by-feature covariance. Sorting is synchronous, with cancellation
checks before and after it. Nonfinite arithmetic is rejected explicitly.

These are numerical capabilities, not a statistical admission rule. The caller
owns the estimand, independent sampling units, alignment, axis identifiability,
deletion failures, interval coverage qualification and multiplicity claims.
In particular, neither a successful fit nor an observed spectral gap proves
that intervals for population axes are valid. There is no percentile fallback,
BCa correction, automatic admission or implicit exclusion of failed draws.

Shared JVM/Scala.js tests compare every three-unit bootstrap multiset and its
valid deletions against literal direct refits. Fifty-six scalar fixtures retain
explicit estimates and draw standard errors with independently computed R
type-7 bounds; additional tests check analytic mean jackknives, equivariance,
extreme scales and strict failure accounting. Published-artifact smoke exercises
the new public deletion, SE and interval entry points.

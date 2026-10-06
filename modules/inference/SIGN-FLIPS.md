# Row and block sign-flip actions

`RowSignFlip` is an immutable vector of row multipliers restricted to -1 and +1.
Its smart constructor rejects an empty vector or any other multiplier. `applyTo`
checks the input row count and multiplies every column of a row by the same sign.
`compose` combines two compatible actions; applying an action twice is identity.

`SignFlipAction.forPartition` draws one independent sign for each `RowPartition`
group and expands it to all rows in that group. Groups may be noncontiguous.
`SignFlipAction.independent` uses one group per row. Neither action requires
balanced signs: all-positive and all-negative draws are legitimate.

```scala
import multivar.inference.*
import resample4s.kernel.Seed

val draw = for
  rows <- RowCount(4)
  partition <- RowPartition.from(rows, Vector(Vector(0, 2), Vector(1, 3)))
  replicate <- ReplicateId(3)
  signs <- SignFlipAction.forPartition(partition).draw(Seed.fromLong(17L), replicate)
yield signs
```

Draws use existing resample4s two-element permutations, independently derived
for each group. Stream domain 1003 is reserved for this action; replicate
identity and the root seed determine each draw independently of evaluation
order. Existing permutation streams are unchanged. Retain the partition, seed,
replicate, explicit resulting signs and provider version for replay. The public
surface is recorded separately from ordinary estimator facade exports in
`tools/public-surface/sign-flip-api.txt`.

These are computational actions, not a null-model compiler, significance test
or validity certificate. A caller must establish that the joint distribution of
each flipped block vector is symmetric under its complete sign reversal and
that the joint distribution is invariant under independent block flips.
Independence between blocks plus joint central symmetry within each block is a
sufficient model assumption. Scalar marginal symmetry alone is insufficient.
The caller also owns centering, nuisance treatment, statistic, resampling budget,
sequential inference and multiplicity. No normalization, centering, row
permutation, feature-wise randomization or distributional pretest is performed
by these actions.

The distinction between exchangeability and independent/symmetric errors is
explained by [Winkler et al. (2014)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4010955/).
That reference does not by itself certify a downstream fitted-response pipeline
or an estimated-deflation procedure.

## Verification, 6 September 2026

`compileAll` and `testAll` pass in
`/private/tmp/multivar-sign-final-tests.log`: 576 core tests on each platform,
48 IR JVM/46 IR JS, and 82 inference tests on each platform (1,410 total).
The five shared sign-action tests cover exact all-column transformation,
composition/involution over the complete three-row orbit, invalid values/shapes,
Scala constructor access, noncontiguous group identity, deterministic replay
independent of evaluation order, and unconstrained four-row sign patterns.
Constructor invariants additionally protect JVM callers from bypassing the
Scala-private smart-constructor boundary.

Canonical `smokeCheck` and `smoke/runMain multivar.smoke.SignFlipSmoke` pass in
`/private/tmp/multivar-sign-final-smoke.log` using the documented JDK22 Scaladoc
workaround. Documentation is not disabled. Existing repeated-classpath
Scaladoc warnings remain; the JDK25 documentation issue has a separate ticket.
The consumer uses publishedLocal artifacts without source-project dependencies.

Both ordinary and sign-flip public-surface checks pass in
`/private/tmp/multivar-sign-final-surface.log`. The ordinary estimator surface
has no diff; the explicit expert action additions are recorded in the separate
snapshot. No permutation action was modified. These checks establish the
primitive's behavior and publication graph, not downstream statistical validity.

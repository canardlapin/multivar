# Typed Multivariate Inference

## Decision

`multivar-inference` owns perturbation-based inference for fitted multivariate
structures. Its dependency boundary is:

```text
Gale -> multivar-core -> multivar-inference
                       -> Resample4s
```

Gale owns matrices, operators, general solvers, and numerical certificates.
Multivar core owns fitted geometry and statistical interpretation. Resample4s
owns lawful resampling plans, canonical label structures, and deterministic
seed-path derivation. Inference composes those capabilities into typed targets,
null actions, Monte Carlo programs, latent units, evidence, and provenance.

The module must not acquire ScalaFIM, neuroimaging, dataset, scheduler, storage,
Spark, filesystem, or platform-runtime dependencies. It must not recreate a
private numerical-solver or random-number layer.

The central abstraction is:

> An inferential program combines an ordered fitted structure, an invariant
> target statistic, a design-preserving null action, and a continuation rule.

Bootstrap stability is a second interpretation of the same fitted structure.
It shares units, refitting, deterministic replicate planning, and provenance,
but it is not significance testing.

## Scientific contracts retained from the ScalaFIM prototype

The port preserves the prototype's scientific design while changing ownership
and infrastructure:

1. Ordered latent structures share a ladder, not a universal statistic.
   Removal of selected structure remains family-specific.
2. Evidence attaches to latent units. A separated axis is orientable; tied or
   near-tied axes form an unoriented subspace.
3. Row permutation, paired-block independence, within-group permutation,
   whole-cluster permutation, bootstrap sampling, and nuisance residual
   randomization are different lawful actions.
4. Exact reduced updates are optional proved capabilities. They never silently
   replace full refits.
5. Exact, conditional, asymptotic, and heuristic validity remain explicit.
6. Bootstrap ratios, selection frequencies, and principal angles are stability
   evidence, not automatically p-values.

The R fixture bundle under `modules/inference/fixtures/v1` remains the external
scientific oracle. Historical wire and fixture identifiers beginning with
`scalafim-` are compatibility data and are intentionally unchanged.

## Typed plan and execution boundary

`InferenceSpec[F, K, N, D]` is immutable plan data. It contains no functions,
fitted objects, matrices, open resources, or scheduler handles. It is prepared
for a future explicit wire codec; it is not described as serializable until
such a codec and schema exist.

`InferenceCompiler` requires typed evidence for the fit/target, target/null,
and null/design relationships. Static callers receive missing-given errors for
unsupported combinations. Dynamic callers receive
`InferenceError.UnsupportedProblem`.

`InferenceProgram` is interpreted separately. `InferenceExecutor.runProgram`
makes the compiled design, null hypothesis, seed, unit policy, and sequential
Monte Carlo policy authoritative. The operational configuration adds the
global budget and maximum ladder length. A mismatch is a typed refusal.
`runProgram` is specifically the significance-ladder interpreter: it accepts
`RequestedEvidence.SignificanceOnly` and returns
`InferenceError.UnsupportedEvidence` for stability-only or combined requests.
Those requests remain valid plan data for bootstrap/stability interpreters;
compilation does not claim that one executor implements every evidence mode.
The initial significance protocols expose axis statistics only. If a
`UnitPolicy` actually forms a multi-axis subspace, execution returns
`InferenceError.UnsupportedUnitPolicy` until that protocol provides a
rotation-invariant subspace statistic; it never relabels axis-wise evidence as
subspace evidence.

Conditioned plans remain data-free by naming their conditioning resource.
Execution must not pretend that a name is a projector: callers prepare
`ResidualPermutationAction` with an actual `RowProjector` and optional
`RowWhitening`. The ordinary ladder executor refuses conditioned plans until
those resources are supplied through an appropriate prepared interpreter.

## Deterministic resampling

All resampling uses Resample4s:

- `Seed` is the root deterministic seed value;
- `StreamPath` derives replicate and independent-lane streams;
- `Seed.derivationAlgorithm` records the versioned `seed-path/v1` contract;
- `Permutation`, `Draw`, `Labels`, `Blocks`, `Groups`, `Strata`, and
  `Split[Selection]` are the canonical action and design values;
- replicate identity, not scheduling order, determines each action.

Multivar inference owns semantic lowering from `ResamplingDesign` to a lawful
Resample4s plan. Resample4s remains general and has no Multivar-specific
targets, protocols, or provenance schema.

## Numerical interpretation

Inference uses Multivar core solver capabilities for rank-revealing fitting and
deflation. Two operations deliberately require a different spectral contract:

- ordered-root inference retains the complete spectrum, including numerical
  zero roots, because zero-root subspaces remain inferential units;
- principal angles retain zero singular values of the cross-basis matrix,
  because a zero singular value represents an angle of `pi / 2`.

Those full-spectrum operations call Gale's exact `SingularSelection.All`
boundary directly. They must not reinterpret a rank-revealing result as if it
contained discarded zero singular values.

## Monte Carlo and result contracts

Fixed and sequential Monte Carlo use the Phipson-Smyth correction:

```text
p = (1 + exceedances) / (draws + 1)
```

Sequential early non-rejection is a policy over the same accumulator. Receipts
record allocation, consumption, exceedance count, batch schedule, stopping
boundary and reason, p-value, Monte Carlo standard error, seed, random
algorithm, and seed-derivation algorithm.

The ladder stops at the first non-selection. Global budget is a separate pure
state transition, and unused draws roll forward. Budget exhaustion is an
explicit result state.

Evidence uses:

```scala
enum Evidence[+A]:
  case NotRequested
  case Unavailable(reason: UnavailableReason)
  case Computed(value: A)
```

Absence, refusal, and a computed result are never collapsed into one sentinel.

## Built-in protocol scope

The initial supported families are:

- PCA variance roots under independent within-column row actions;
- PLSC covariance roots under an explicitly chosen break-X or break-Y null;
- CCA canonical correlations under paired independence;
- non-negative generalized-eigen roots;
- held-out RRR predictive gain;
- one named CPCA block under a declared constraint-relative null;
- multiblock consensus roots under independent non-reference-block actions;
- method-native feature evidence with explicit multiplicity correction.

Merely exposing a spectrum is insufficient to claim inferential support. A new
family needs a target, a compatible null and design, a removal or predictive
continuation, a validity statement, and independent evidence.

## Verification

Every numerical or inferential change should provide the applicable layers:

1. smart-constructor and invalid-state tests;
2. algebraic laws for action shape, determinism, invariance, and reduction;
3. independent analytic or R fixtures;
4. end-to-end scenario tests through compilation and execution.

The focused cross-platform gate is:

```sh
sbt inferenceJVM/test inferenceJS/test
```

Module completion additionally requires:

```sh
./tools/check-boundary.sh
sbt compileAll testAll
./tools/inference-surface.sh --check
sbt mimaCheck docsCheck smokeCheck
```

The published-consumer smoke must compile an `InferenceSpec` from the
publishedLocal `multivar-inference` artifact without a source-project
`dependsOn`.

## Non-goals

- Reimplementing fitted PCA, PLSC, CCA, RRR, CPCA, or multiblock models.
- Owning general resampling primitives already supplied by Resample4s.
- Owning general numerical solvers already supplied by Gale.
- A universal p-value abstraction that erases target/null provenance.
- Approximate inference presented as exact.
- Feature significance inferred from ordinary bootstrap stability.
- Distributed runtime dependencies in shared code.

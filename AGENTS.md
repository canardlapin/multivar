# AGENTS.md

Guidance for coding agents working in **multivar**, a general Scala 3 library
for typed multivariate analysis. Read [README.md](README.md) for the public
surface and [docs/plans/multivar-package-hierarchy.md](docs/plans/multivar-package-hierarchy.md)
before changing ownership boundaries.

## Layout and dependency boundary

- `modules/core` is a JVM/Scala.js `crossProject` containing the `multivar.*`
  semantic, numerical, solver, model-family, and workflow packages.
- `modules/ir` is a JVM/Scala.js `crossProject` containing `multivar.ir.*`,
  JSON codecs, schemas, and conformance fixtures.
- Core depends on the pinned Gale source revision in `build.sbt`. IR depends on
  core. Keep this graph acyclic.
- Gale owns portable matrices, linear operators, general numerical solvers, and
  first-order stopping certificates. Multivar owns semantic lowering and
  statistical interpretation; do not recreate a `multivar.numerics` layer.
- Do not introduce `scalafim.*`, neuroimaging, dataset, scheduler, storage, or
  platform-runtime dependencies. Those belong in downstream adapter projects.
- Historical `scalafim-*` wire identifiers are compatibility contracts. Do not
  rename an existing schema or media type without a new version and migration
  path.

## Build and tests

- Toolchain: Scala 3.7.4, sbt 1.10.5, MUnit.
- Run `sbt compileAll testAll` before declaring a feature complete.
- Shared tests must pass on both JVM and Scala.js. Platform-only tests belong in
  the corresponding `jvm` or `js` source tree.
- Keep `-deprecation -feature -unchecked` warning-clean.
- Numerical changes need law, adversarial, or independent-oracle evidence.
  Floating-point assertions require explicit tolerances.

## Scala style and numerical work

- Use Scala 3 significant indentation with two spaces.
- Prefer closed enums, precise records, opaque domain values, smart
  constructors, exhaustive matches, and typed errors over strings or sentinels.
- Keep public APIs immutable and inspectable. Primitive arrays and `while`
  loops are appropriate inside allocation-sensitive kernels.
- Keep model declarations, solver requests, achieved guarantees, certificates,
  provenance, and fitted capabilities distinct.
- Do not weaken invariants or add escape hatches merely to preserve a call
  site. Make compatibility explicit and tested.

## GitHub identity

- The canonical repository is `canardlapin/multivar`.
- Use repo-local identity `canardlapin
  <307091466+canardlapin@users.noreply.github.com>` and
  `github.account=canardlapin`.
- Git pushes use the `github-canardlapin` SSH alias. GitHub API operations use
  an isolated `gh-canardlapin` profile or a repository-local wrapper; never
  switch the machine-wide GitHub account.

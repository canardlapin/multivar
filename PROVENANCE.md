# Provenance

`multivar` was extracted from
[`canardlapin/scalafim`](https://github.com/canardlapin/scalafim) at committed
source revision `7c590372dd3cea733f5835b9ca1bb6f62e92292f` on 2026-07-26.

The extraction boundary contains:

- ScalaFIM's former `modules/multivar` sources, tests, and fixtures, now in
  `modules/core`;
- ScalaFIM's former `modules/multivar-ir` schemas, codecs, tests, and
  conformance corpus, now in `modules/ir`;
- the multivar design records and independent R/LowRankModels.jl oracle tools.

The extraction deliberately:

- renamed Scala packages from `scalafim.multivar.*` to `multivar.*`;
- removed `LocusSelectionAdapter`, which depends on ScalaFIM's finite-locus
  model and remains a downstream ScalaFIM responsibility;
- moved the reusable first-order optimization contracts, algorithms, typed
  failures, and tests from `scalafim.linalg` into `gale.optim` at immutable
  Gale revision `d55fe2f97196a76ab7879e1a12f1e92403aeba06`;
- replaced the extracted ScalaFIM matrix/map bridge with Gale-native `DMat` and
  `DoubleLinearOperator` usage throughout;
- retained existing `scalafim-*` schema names, media types, and oracle fixture
  versions as stable wire-compatibility identifiers.

No source under `modules/core` or `modules/ir` depends on a `scalafim.*`
package or a copied `multivar.numerics` implementation. The standalone
repository is the authority for future multivar development; ScalaFIM consumes
a pinned standalone revision through its own adapter boundary.

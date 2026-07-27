# Generalized low-rank models

A generalized low-rank model (GLRM) can combine feature-specific losses,
missing cells, and a shared latent representation. The feature domains and
losses are part of the model, so fitting begins with a program rather than a
single matrix call.

## Declare the observed data

Build an `ObservationPattern` whose cells are observed, missing, structurally
inapplicable, or censored. Then build a `GlrmFeatureLayout`. Each
`GlrmFeatureSpec` gives one feature:

- a stable `GlrmFeatureId`;
- a domain such as real, non-negative, binary, ordinal, or categorical; and
- an `EntryLoss` compatible with that domain.

`GeneralizedLowRankProgram.from` checks the feature layout against the
observation pattern, checks each observed value against its domain, and records
the missingness statement and prediction target.

## Fit through the lifecycle

The fitting sequence is:

1. `GeneralizedLowRankLifecycle.declare` binds the objective to its
   mathematical contract.
2. `GeneralizedLowRankLifecycle.compile` receives the latent space, initial
   factors, penalties, and first-order configuration.
3. `GeneralizedLowRankLifecycle.fit` runs the admitted PALM problem and returns
   the factors, objective value, convergence receipt, and fitted encoder.

This separation keeps three facts distinct: the model that was requested, the
solver plan that was admitted, and the guarantee achieved by the fitted
result.

## Encode a new row

`FittedLatentEncoder` estimates a latent code from a new row while holding the
fitted decoder fixed. It accepts partial observations and reports whether the
available features identify the requested code. Decoding and encoding are
separate operations; a decoder does not imply that a new row can be identified
from every subset of features.

Use PCA for complete real-valued data with squared loss and Euclidean geometry.
Use a GLRM when feature domains, losses, or missingness change what is being
estimated.

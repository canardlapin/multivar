# Choose a method

Start with the question the analysis must answer. The table separates methods
with direct dense entry points from methods that require an explicit model
definition.

## Direct dense analyses

| Question | Start with | Result to inspect |
| --- | --- | --- |
| Which directions summarize one matrix as supplied? | `Svd.fit` | `SvdFit` |
| Which directions explain centered variation? | `Pca.fit` | `PcaFit` |
| Which paired directions have high covariance? | `Plsc.fit` | `PlscFit` |
| Which paired directions have high correlation? | `Cca.fit` | `CcaFit` |
| How can `y` be predicted from `x` with limited rank? | `ReducedRankRegression.fit` | `ReducedRankRegressionFit` |
| How can `y` be predicted from `x` through PLS components? | `PlsRegression.fit` | `PlsRegressionFit` |
| Which linear directions separate labelled classes? | `FisherDiscriminant.fit` | `FisherDiscriminantFit` |
| How does PCA change under row and feature metrics? | `Gpca.fit` | `GpcaFit` |
| What remains inside supplied row and feature design spaces? | `Cpca.fit` | `CpcaFit` |
| How can a kernel eigensystem be approximated from landmarks? | `Nystrom.fit` or `fitRbf` | `NystromFit` |

## Explicit model definitions

| Question | Build | Fit or evaluate |
| --- | --- | --- |
| How do sparse or smooth penalties change a factor? | `RankOneStructuredFactorization` | `solve` |
| Should several structured factors be greedy or joint? | `GreedyStructuredFactorization` or `JointRankKStructuredFactorization` | `solve` |
| How do feature-specific losses and missing cells share a latent model? | `GeneralizedLowRankProgram` | `GeneralizedLowRankLifecycle` |
| How do declared effect and residual operators differ canonically? | `CanonicalEffectProblem` | `fit` or `fitSpectrum` |
| How do frame constraints change a canonical effect? | `ConstrainedCanonicalProblem` | `fit` |
| How do independent views maximize a declared association? | `DirectSumStudy` and `MaximizeAssociation` | `MultisetAssociation.fit` |
| How do aligned blocks share row scores under different losses? | `AlignedSharedScoreGlrm` | `evaluate` and `fittedEncoder` |

## Reuse a fitted model

- Use `transform` for SVD, PCA, and Nyström, and `inverseTransform` or
  `reconstruct` to return from component coordinates to feature coordinates.
- Use `project` for GPCA and LDA.
- Use `transformX` or `transformY` for PLSC and CCA.
- Use `predict` for reduced-rank regression and PLS regression; use `transform`
  on a PLS fit for X-scores.
- Use fitted encoders for GLRM and aligned multiblock models.
- Use synthesis capabilities only when a reconstruction contract has been
  fitted explicitly.

Read the task page before moving to an explicit model definition. The required
inputs define the model; they are not wrapper code around the solver.

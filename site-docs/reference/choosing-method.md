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
| Which linear directions separate labelled classes? | `Lda.fit` | `LdaFit` |
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

- Use `project` for SVD, PCA, GPCA, and LDA.
- Use `projectX` or `projectY` for PLSC and CCA.
- Use `predict` for reduced-rank regression.
- Use `transform` for Nyström.
- Use fitted encoders for GLRM and aligned multiblock models.
- Use synthesis capabilities only when a reconstruction contract has been
  fitted explicitly.

Read the task page before moving to an explicit model definition. The required
inputs define the model; they are not wrapper code around the solver.

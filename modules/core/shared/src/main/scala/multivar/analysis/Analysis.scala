package multivar
package analysis

/** Curated ordinary API for dense multivariate analysis.
  *
  * This package is an editorial statement, not a mirror of the ownership tree.
  * Import it when you want the supported analyst surface:
  *
  * {{{
  * import multivar.analysis.*
  * }}}
  *
  * Semantic owners under `multivar.core`, `multivar.family.*`, and the rest of
  * the hierarchy remain the authoritative homes of each type. Expert surfaces
  * live in [[multivar.advanced]] and [[multivar.syntax.unsafe]].
  */
export multivar.core.{MultivarError, PreprocessSpec, VarianceConvention}
export multivar.family.spectral.{Gpca, GpcaCentering, GpcaFit, Pca, PcaFit, Svd, SvdFit}
export multivar.family.paired.{
  Cca,
  CcaFit,
  Plsc,
  PlscFit,
  PlsAlgorithm,
  PlsRegression,
  PlsRegressionFit,
  ReducedRankRegression,
  ReducedRankRegressionFit
}
export multivar.family.canonical.{FisherDiscriminant, FisherDiscriminantFit, WithinScatterPolicy}
export multivar.family.cpca.{Cpca, CpcaFit}
export multivar.family.kernel.{KernelSymmetryPolicy, Nystrom, NystromFit}

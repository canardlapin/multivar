package multivar
package family.spectral


import multivar.core.*
import multivar.capability.*

import gale.linalg.DMat
import gale.linalg.DVec

/** Shared machinery behind [[PcaFit]] and [[SvdFit]].
  *
  * The two results deliberately do not share a public supertype. They agree on how
  * the numbers are computed and disagree on what the numbers mean: an SVD of
  * uncentered data decomposes a sum of squares about the origin, which is not a
  * variance. Keeping the arithmetic here and the vocabulary in each result lets the
  * public surface of each be read off its own class.
  *
  * This carrier is public only because Scaladoc crashes when a public result's
  * constructor mentions a package-private type. It is not part of the ordinary
  * surface: `multivar.analysis` does not export it, and `tools/public-surface`
  * does not snapshot it.
  */
final class SpectralCore private[multivar] (
    val frameTransform: FittedFrameTransform,
    val result: SvdResult,
    val totalSumSquares: Double,
    val preprocessing: Option[ColumnAffineSummary]
):
  def rows: Int =
    frameTransform.trainingRowSpace.descriptor.size

  def nFeatures: Int =
    frameTransform.featureSpace.descriptor.size

  def requestedComponents: Int =
    frameTransform.diagnostics.requestedComponents.value

  def effectiveComponents: Int =
    frameTransform.diagnostics.effectiveComponents

  def scores: DMat =
    frameTransform.scores

  def loadings: DMat =
    frameTransform.weights

  def components: DMat =
    frameTransform.weights.transpose

  def singularValues: DVec =
    result.singularValues

  def squaredSingularValues: DVec =
    mapSingularValues(value => value * value)

  def inertiaRatio: DVec =
    mapSingularValues(value => value * value / totalSumSquares)

  /** Only meaningful for a fit that centres and has at least two observations, which
    * is why `Pca.fit` requires both and `SvdFit` does not offer this.
    */
  def explainedVariance: DVec =
    val denominator = (rows - 1).toDouble
    mapSingularValues(value => value * value / denominator)

  def transform(input: MatrixView): Either[MultivarError, DMat] =
    frameTransform.project(input)

  def inverseTransform(scores: DMat): Either[MultivarError, DMat] =
    synthesis.flatMap(_.synthesizeOriginal(scores).map(_.values))

  def reconstruct(input: MatrixView): Either[MultivarError, DMat] =
    synthesis.flatMap(_.reconstruct(input).map(_.values))

  /** Synthesis is the transpose of the fitted loadings, which are orthonormal by
    * construction of the SVD. Building it on demand keeps a tolerance check that
    * should never fire out of the fitting path, and caching it keeps the check from
    * being repeated on every reconstruction.
    */
  private lazy val synthesis: Either[MultivarError, FittedBidirectionalTransform] =
    frameTransform.withOrthonormalSynthesis()

  private def mapSingularValues(f: Double => Double): DVec =
    val values = result.singularValues
    val out = new Array[Double](values.length)
    var index = 0
    while index < values.length do
      out(index) = f(values(index))
      index += 1
    GaleNumerics.vectorFromArray(out)

/** A fitted singular value decomposition of the preprocessed input.
  *
  * SVD does not centre, so it decomposes the sum of squares about the origin rather
  * than the variance about the mean. The retained shares are reported as inertia to
  * keep that distinction visible.
  */
final class SvdFit private[multivar] (private val core: SpectralCore):
  /** Training observations in component coordinates, one row per observation. */
  def scores: DMat =
    core.scores

  /** Right singular vectors, one component per column (features by components). */
  def loadings: DMat =
    core.loadings

  /** The loadings transposed, one component per row. */
  def components: DMat =
    core.components

  def singularValues: DVec =
    core.singularValues

  /** Sum of squares about the origin carried by each retained component. */
  def inertia: DVec =
    core.squaredSingularValues

  /** Each component's share of the total sum of squares of the preprocessed input. */
  def inertiaRatio: DVec =
    core.inertiaRatio

  /** Column offsets removed by the fitted preprocessing, in `(x - center) / scale`
    * form. `None` when the preprocessing collapsed a column and so has no such
    * description.
    */
  def center: Option[DVec] =
    core.preprocessing.map(_.center)

  /** Column divisors applied by the fitted preprocessing. See [[center]]. */
  def scale: Option[DVec] =
    core.preprocessing.map(_.scale)

  def nFeatures: Int =
    core.nFeatures

  def requestedComponents: Int =
    core.requestedComponents

  def effectiveComponents: Int =
    core.effectiveComponents

  /** Express new observations in component coordinates, applying the fitted
    * preprocessing first.
    */
  def transform(input: DMat): Either[MultivarError, DMat] =
    core.transform(MatrixView.dense(input))

  def transform(input: MatrixView): Either[MultivarError, DMat] =
    core.transform(input)

  /** Map component coordinates back to original feature coordinates. */
  def inverseTransform(scores: DMat): Either[MultivarError, DMat] =
    core.inverseTransform(scores)

  /** The rank-`effectiveComponents` approximation of the given observations, in
    * original feature coordinates.
    */
  def reconstruct(input: DMat): Either[MultivarError, DMat] =
    core.reconstruct(MatrixView.dense(input))

  def reconstruct(input: MatrixView): Either[MultivarError, DMat] =
    core.reconstruct(input)

object SvdFit:
  private[multivar] def coreOf(fit: SvdFit): SpectralCore =
    fit.core

object Svd:
  def fit(
      input: DMat,
      components: Int
  ): Either[MultivarError, SvdFit] =
    fit(input, components, PreprocessSpec.Pass)

  def fit(
      input: DMat,
      components: Int,
      preproc: PreprocessSpec
  ): Either[MultivarError, SvdFit] =
    ComponentCount(components).flatMap: checked =>
      fit(MatrixView.dense(input), checked, preproc)

  def fit(
      input: MatrixView,
      components: ComponentCount,
      preproc: PreprocessSpec = PreprocessSpec.Pass,
      solver: SvdSolver = DenseSolvers.svd
  ): Either[MultivarError, SvdFit] =
    fitFrame(input, components, preproc, solver, "svd").map(SvdFit(_))

/** A fitted principal component analysis of the preprocessed input.
  *
  * Every reported quantity describes the preprocessed data, not the raw input. Under
  * the default column centring the two agree about variation; under `Standardize()`
  * the variances are those of the standardized columns, which is the correlation-PCA
  * convention.
  */
final class PcaFit private[multivar] (private val core: SpectralCore):
  /** Training observations in principal component coordinates. */
  def scores: DMat =
    core.scores

  /** Principal axes, one component per column (features by components). */
  def loadings: DMat =
    core.loadings

  /** The loadings transposed, one component per row. */
  def components: DMat =
    core.components

  def singularValues: DVec =
    core.singularValues

  /** Variance of the preprocessed data along each retained component, under the
    * sample convention that divides by `n - 1`.
    */
  def explainedVariance: DVec =
    core.explainedVariance

  /** Each component's share of the total variance of the preprocessed data. The
    * shares sum to one only when every component is retained.
    */
  def explainedVarianceRatio: DVec =
    core.inertiaRatio

  /** Column means removed by the fitted preprocessing, in `(x - center) / scale`
    * form. `None` when the preprocessing collapsed a column and so has no such
    * description.
    */
  def center: Option[DVec] =
    core.preprocessing.map(_.center)

  /** Column divisors applied by the fitted preprocessing. See [[center]]. */
  def scale: Option[DVec] =
    core.preprocessing.map(_.scale)

  def nFeatures: Int =
    core.nFeatures

  def requestedComponents: Int =
    core.requestedComponents

  def effectiveComponents: Int =
    core.effectiveComponents

  /** Express new observations in component coordinates, applying the fitted
    * preprocessing first.
    */
  def transform(input: DMat): Either[MultivarError, DMat] =
    core.transform(MatrixView.dense(input))

  def transform(input: MatrixView): Either[MultivarError, DMat] =
    core.transform(input)

  /** Map component coordinates back to original feature coordinates. */
  def inverseTransform(scores: DMat): Either[MultivarError, DMat] =
    core.inverseTransform(scores)

  /** The rank-`effectiveComponents` approximation of the given observations, in
    * original feature coordinates.
    */
  def reconstruct(input: DMat): Either[MultivarError, DMat] =
    core.reconstruct(MatrixView.dense(input))

  def reconstruct(input: MatrixView): Either[MultivarError, DMat] =
    core.reconstruct(input)

object PcaFit:
  private[multivar] def coreOf(fit: PcaFit): SpectralCore =
    fit.core

object Pca:
  /** Fit ordinary dense PCA while validating the requested component count.
    *
    * This is the direct entry point for dense data. It expands exactly to the
    * `MatrixView` and `ComponentCount` API below.
    */
  def fit(
      input: DMat,
      components: Int
  ): Either[MultivarError, PcaFit] =
    fit(input, components, PreprocessSpec.Center, DenseSolvers.svd)

  def fit(
      input: DMat,
      components: Int,
      preproc: PreprocessSpec
  ): Either[MultivarError, PcaFit] =
    fit(input, components, preproc, DenseSolvers.svd)

  def fit(
      input: DMat,
      components: Int,
      preproc: PreprocessSpec,
      solver: SvdSolver
  ): Either[MultivarError, PcaFit] =
    ComponentCount(components).flatMap: checked =>
      fit(MatrixView.dense(input), checked, preproc, solver)

  /** Fit PCA.
    *
    * At least two observations are required. PCA reports variances, and a variance
    * estimated from one observation is not a number this API is willing to invent.
    * Use `Svd.fit` for a single-observation decomposition.
    */
  def fit(
      input: MatrixView,
      components: ComponentCount,
      preproc: PreprocessSpec = PreprocessSpec.Center,
      solver: SvdSolver = DenseSolvers.svd
  ): Either[MultivarError, PcaFit] =
    if input.rows < 2 then Left(MultivarError.InsufficientRows("pca observations", 2, input.rows))
    else fitFrame(input, components, preproc, solver, "pca").map(PcaFit(_))

private def fitFrame(
    input: MatrixView,
    components: ComponentCount,
    preproc: PreprocessSpec,
    solver: SvdSolver,
    method: String
): Either[MultivarError, SpectralCore] =
  for
    fitted <- preproc.fit(input)
    transformed <- fitted.transform(input)
    stats <- transformed.columnStats
    total <- requireInertia(stats.totalSumSquares, method)
    svd <- solver.decompose(transformed, components)
    _ <- requireComponents(svd)
    transform <- FittedFrameTransform.fromTraining(
      input,
      svd.v,
      fitted,
      method,
      components,
      Some(svd.singularValues)
    )
  yield SpectralCore(transform, svd, total, ColumnAffineSummary.of(fitted))

/** The denominator of every reported share. A successful decomposition of data with
  * no spread would leave it at zero, so it is checked rather than assumed.
  */
private def requireInertia(total: Double, method: String): Either[MultivarError, Double] =
  if total > 0.0 && total.isFinite then Right(total)
  else Left(MultivarError.SolverFailed(s"$method found no variation in the preprocessed input"))

private def requireComponents(svd: SvdResult): Either[MultivarError, Unit] =
  if svd.singularValues.length == 0 then Left(MultivarError.SolverFailed("no singular values above tolerance"))
  else Right(())

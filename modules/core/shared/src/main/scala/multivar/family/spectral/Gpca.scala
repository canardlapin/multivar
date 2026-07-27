package multivar
package family.spectral

import multivar.capability.FittedFrameTransform
import multivar.core.*

import gale.linalg.DMat
import gale.linalg.DVec

/** How dense GPCA removes a constant direction before the generalized eigenproblem.
  *
  * Generic column preprocessing and GPCA row centering are distinct. This policy
  * lowers onto the same arithmetic as the semantic diagram centering plans without forcing every
  * dense caller through a typed duality diagram.
  */
enum GpcaCentering:
  /** Leave the table as supplied. */
  case None

  /** Subtract ordinary (uniform) column means. */
  case Ordinary

  /** Subtract means weighted by a non-negative row measure induced by the row metric. */
  case ByRowMeasure

  /** Remove the constant direction in the row-metric inner product. */
  case OrthogonalToConstant

  /** Assert the table is already centered; apply no further shift. */
  case AlreadyCentered

  /** Ordinary centering under an identity row metric; typed error otherwise. */
  case Auto

/** Opaque dense GPCA result. Typed operator internals open through [[multivar.advanced]]. */
final class GpcaFit private[multivar] (
    private val operator: GpcaOperatorFit[? <: SemanticSpace, ? <: SemanticSpace, ? <: SemanticSpace],
    private val frame: FittedFrameTransform,
    private val featureCentre: Option[DVec]
):
  def scores: DMat = frame.scores

  /** Generalized component weights, one component per column. */
  def weights: DMat = frame.weights

  def eigenvalues: DVec = operator.generalizedEigenvalues

  def singularValues: DVec = operator.singularValues

  /** Training feature centre applied before projection, when centering was fitted. */
  def center: Option[DVec] = featureCentre

  def transform(input: DMat): Either[MultivarError, DMat] =
    transform(MatrixView.dense(input))

  def transform(input: MatrixView): Either[MultivarError, DMat] =
    frame.project(input)

object GpcaFit:
  private[multivar] def operatorOf(
      fit: GpcaFit
  ): GpcaOperatorFit[? <: SemanticSpace, ? <: SemanticSpace, ? <: SemanticSpace] =
    fit.operator

  private[multivar] def frameOf(fit: GpcaFit): FittedFrameTransform =
    fit.frame

object Gpca:
  def fit(
      input: DMat,
      components: Int,
      rowMetric: MetricSpec,
      featureMetric: MetricSpec,
      centering: GpcaCentering = GpcaCentering.Auto
  ): Either[MultivarError, GpcaFit] =
    for
      checked <- ComponentCount(components)
      _ <-
        if rowMetric.dim == input.rows then Right(())
        else Left(MultivarError.MetricShapeMismatch(IndexAxis.Row, input.rows, rowMetric.dim))
      _ <-
        if featureMetric.dim == input.cols then Right(())
        else Left(MultivarError.MetricShapeMismatch(IndexAxis.Feature, input.cols, featureMetric.dim))
      fitted <- fitCentering(centering, rowMetric, MatrixView.dense(input))
      prepared <- fitted.transform(MatrixView.dense(input))
      rows <- MvSpace.of("gpca.rows", SpaceRole.Samples, input.rows)
      features <- MvSpace.of("gpca.features", SpaceRole.Observed, input.cols)
      problem <- DynamicGpcaProblem.from(
        prepared,
        rows,
        features,
        rowMetric,
        featureMetric,
        ValueIdentity.source(ValueId.unsafe("gpca.input")),
        SemanticProvenance.source("dense-gpca")
      )
      operator <- problem.fit(checked)
      weights <- gpcaWeights(operator)
      transform <- FittedFrameTransform.fromTraining(
        MatrixView.dense(input),
        weights,
        fitted,
        "gpca",
        checked,
        Some(operator.singularValues)
      )
      centre =
        centering match
          case GpcaCentering.None | GpcaCentering.AlreadyCentered => None
          case _ => ColumnAffineSummary.of(fitted).map(_.center)
    yield new GpcaFit(operator, transform, centre)

  def fit(
      input: DMat,
      components: Int,
      rowWeights: DVec,
      featureWeights: DVec
  ): Either[MultivarError, GpcaFit] =
    fit(input, components, rowWeights, featureWeights, GpcaCentering.Auto)

  def fit(
      input: DMat,
      components: Int,
      rowWeights: DVec,
      featureWeights: DVec,
      centering: GpcaCentering
  ): Either[MultivarError, GpcaFit] =
    for
      rowMetric <- MetricSpec.diagonal(rowWeights)
      featureMetric <- MetricSpec.diagonal(featureWeights)
      fit <- fit(input, components, rowMetric, featureMetric, centering)
    yield fit

  /** Fit the column affine that realizes `centering` against `rowMetric`. */
  private[multivar] def fitCentering(
      centering: GpcaCentering,
      rowMetric: MetricSpec,
      input: MatrixView
  ): Either[MultivarError, FittedPreprocessor] =
    centering match
      case GpcaCentering.None | GpcaCentering.AlreadyCentered =>
        PreprocessSpec.Pass.fit(input)
      case GpcaCentering.Ordinary =>
        PreprocessSpec.Center.fit(input)
      case GpcaCentering.Auto =>
        rowMetric match
          case _: MetricSpec.Identity => PreprocessSpec.Center.fit(input)
          case _ =>
            Left(
              MultivarError.InvalidRowGeometry(
                "GpcaCentering.Auto requires an identity row metric; choose Ordinary, ByRowMeasure, OrthogonalToConstant, None, or AlreadyCentered explicitly"
              )
            )
      case GpcaCentering.ByRowMeasure =>
        measureWeights(rowMetric).flatMap(weightedCenter(input, _))
      case GpcaCentering.OrthogonalToConstant =>
        orthogonalWeights(rowMetric).flatMap(weightedCenter(input, _))

  private def weightedCenter(input: MatrixView, weights: DVec): Either[MultivarError, FittedPreprocessor] =
    if input.cols <= 0 then Left(MultivarError.InvalidDimension("preprocessing input columns", input.cols))
    else if input.rows != weights.length then
      Left(MultivarError.MatrixShapeMismatch(s"centering weights have length ${weights.length}, expected ${input.rows}"))
    else
      input.toDense(StoragePolicy.AllowDense).map { dense =>
        val means = weightedColumnMeans(dense, weights)
        FittedColumnAffine(input.cols, MatrixView.ones(input.cols), MatrixView.negate(means))
      }

  private def measureWeights(rowMetric: MetricSpec): Either[MultivarError, DVec] =
    rowMetric match
      case _: MetricSpec.Identity =>
        Right(uniformWeights(rowMetric.dim))
      case diagonal: MetricSpec.Diagonal =>
        normalizeNonNegative(diagonal.weights, "GPCA row measure")
      case _: MetricSpec.DenseSymmetric | _: MetricSpec.SparseSymmetric =>
        Left(
          MultivarError.InvalidRowGeometry(
            "GpcaCentering.ByRowMeasure requires an identity or diagonal row metric; use OrthogonalToConstant for a general SPD metric"
          )
        )

  private def orthogonalWeights(rowMetric: MetricSpec): Either[MultivarError, DVec] =
    rowMetric match
      case _: MetricSpec.Identity =>
        Right(uniformWeights(rowMetric.dim))
      case diagonal: MetricSpec.Diagonal =>
        normalizeNonNegative(diagonal.weights, "GPCA metric-orthogonal measure")
      case dense: MetricSpec.DenseSymmetric =>
        val ones = GaleNumerics.matrixFromRowMajor(dense.dim, 1, Array.fill(dense.dim)(1.0))
        dense.matvec(ones).flatMap { applied =>
          val raw = applied.copyData
          var denominator = 0.0
          var index = 0
          while index < raw.length do
            denominator += raw(index)
            index += 1
          if !denominator.isFinite || Math.abs(denominator) <= 1e-14 then
            Left(MultivarError.InvalidRowGeometry(s"row metric gives invalid centering denominator $denominator"))
          else
            index = 0
            while index < raw.length do
              raw(index) /= denominator
              index += 1
            Right(GaleNumerics.vectorFromArray(raw))
        }
      case _: MetricSpec.SparseSymmetric =>
        Left(
          MultivarError.InvalidRowGeometry(
            "GpcaCentering.OrthogonalToConstant on a sparse row metric requires SemanticGpca"
          )
        )

  private def uniformWeights(dim: Int): DVec =
    val weight = 1.0 / dim.toDouble
    GaleNumerics.vectorFromArray(Array.fill(dim)(weight))

  private def normalizeNonNegative(weights: DVec, role: String): Either[MultivarError, DVec] =
    val raw = weights.copyData
    var mass = 0.0
    var index = 0
    var error = Option.empty[MultivarError]
    while index < raw.length && error.isEmpty do
      val weight = raw(index)
      if !weight.isFinite then error = Some(MultivarError.NonFiniteValue(role, index, weight))
      else if weight < 0.0 then error = Some(MultivarError.InvalidRowGeometry(s"$role weight $index is negative: $weight"))
      else mass += weight
      index += 1
    error match
      case Some(value) => Left(value)
      case None if !mass.isFinite || mass <= 0.0 =>
        Left(MultivarError.InvalidRowGeometry(s"$role must have finite positive mass, got $mass"))
      case None =>
        index = 0
        while index < raw.length do
          raw(index) /= mass
          index += 1
        Right(GaleNumerics.vectorFromArray(raw))

  private def weightedColumnMeans(matrix: DMat, weights: DVec): DVec =
    val out = new Array[Double](matrix.cols)
    var col = 0
    while col < matrix.cols do
      var acc = 0.0
      var row = 0
      while row < matrix.rows do
        acc += weights(row) * matrix(row, col)
        row += 1
      out(col) = acc
      col += 1
    GaleNumerics.vectorFromArray(out)

private def gpcaWeights(
    fit: GpcaOperatorFit[? <: SemanticSpace, ? <: SemanticSpace, ? <: SemanticSpace]
): Either[MultivarError, DMat] =
  fit.functionalFrame.weights.toDense.left.map:
    case SemanticError.MultivarFailure(error)  => error
    case SemanticError.LinearMapFailure(error) => LinalgErrorAdapter.toMultivarError(error)
    case error                                 => MultivarError.SolverFailed(error.message)

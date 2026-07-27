package multivar
package family.paired

import multivar.core.*

import gale.linalg.DMat

/** Shared preparation for a dense paired (X, Y) analysis.
  *
  * Owns preprocessing, working views, and sample geometry that both the
  * cross-spectral estimators (PLSC, CCA, RRR) and predictive latent estimators
  * (PLS regression) need. Solver-specific semantics stay outside this value.
  * Dense materialization is deferred so storage policy remains the caller's
  * choice at the operator or SIMPLS boundary.
  */
final class PreparedPair private[paired] (
    val xOriginal: MatrixView,
    val yOriginal: MatrixView,
    val xPreprocessor: FittedPreprocessor,
    val yPreprocessor: FittedPreprocessor,
    val xWorking: MatrixView,
    val yWorking: MatrixView,
    val rowMetric: Option[MetricSpec],
    val moments: PairedMoments
):
  /** Response preprocessor required by predictive maps that restore original Y. */
  def invertibleResponse: Either[MultivarError, FittedInvertiblePreprocessor] =
    yPreprocessor match
      case invertible: FittedInvertiblePreprocessor => Right(invertible)
      case other =>
        Left(
          MultivarError.InvalidMap(
            s"predictive paired fit requires an invertible response preprocessor, got ${other.getClass.getName}"
          )
        )

  /** Materialize finite dense working matrices under the caller's storage policy. */
  def workingDense(
      policy: StoragePolicy = StoragePolicy.AllowDense
  ): Either[MultivarError, (DMat, DMat)] =
    for
      xDense <- xWorking.toDense(policy)
      yDense <- yWorking.toDense(policy)
      _ <- MatrixOps.checkFinite("paired predictors", xDense)
      _ <- MatrixOps.checkFinite("paired responses", yDense)
    yield (xDense, yDense)

final case class PairedMoments(sampleCount: Int, xFeatures: Int, yFeatures: Int):
  require(sampleCount > 0, "paired sample count must be positive")
  require(xFeatures > 0, "paired predictor features must be positive")
  require(yFeatures > 0, "paired response features must be positive")

  /** Scale for covariance-style operators that divide by `n - 1`. */
  def covarianceScale: Double =
    1.0 / Math.max(1, sampleCount - 1).toDouble

object PreparedPair:
  /** Prepare a paired analysis with ordinary (not necessarily invertible) response fit. */
  def from(
      x: MatrixView,
      y: MatrixView,
      xPreproc: PreprocessSpec,
      yPreproc: PreprocessSpec,
      rowMetric: Option[MetricSpec] = None,
      policy: StoragePolicy = StoragePolicy.AllowDense
  ): Either[MultivarError, PreparedPair] =
    prepare(x, y, xPreproc, yPreproc, invertibleResponse = false, rowMetric, policy)

  /** Prepare a predictive paired analysis; response preprocessing must be invertible. */
  def fromPredictive(
      x: MatrixView,
      y: MatrixView,
      xPreproc: PreprocessSpec,
      yPreproc: PreprocessSpec,
      rowMetric: Option[MetricSpec] = None,
      policy: StoragePolicy = StoragePolicy.AllowDense
  ): Either[MultivarError, PreparedPair] =
    prepare(x, y, xPreproc, yPreproc, invertibleResponse = true, rowMetric, policy)

  private def prepare(
      x: MatrixView,
      y: MatrixView,
      xPreproc: PreprocessSpec,
      yPreproc: PreprocessSpec,
      invertibleResponse: Boolean,
      rowMetric: Option[MetricSpec],
      policy: StoragePolicy
  ): Either[MultivarError, PreparedPair] =
    if x.rows != y.rows then
      Left(MultivarError.MatrixShapeMismatch(s"paired input expected equal rows, got ${x.rows} and ${y.rows}"))
    else if x.rows <= 0 then Left(MultivarError.InvalidDimension("paired sample count", x.rows))
    else if x.cols <= 0 then Left(MultivarError.InvalidDimension("paired predictor columns", x.cols))
    else if y.cols <= 0 then Left(MultivarError.InvalidDimension("paired response columns", y.cols))
    else
      for
        _ <- rowMetric match
          case Some(metric) if metric.dim != x.rows =>
            Left(MultivarError.MetricShapeMismatch(IndexAxis.Row, x.rows, metric.dim))
          case _ =>
            Right(())
        fittedX <- xPreproc.fit(x)
        fittedY: FittedPreprocessor <-
          if invertibleResponse then yPreproc.fitInvertible(y).map(p => p: FittedPreprocessor)
          else yPreproc.fit(y)
        xp <- fittedX.transform(x, policy = policy)
        yp <- fittedY.transform(y, policy = policy)
      yield
        new PreparedPair(
          x,
          y,
          fittedX,
          fittedY,
          xp,
          yp,
          rowMetric,
          PairedMoments(x.rows, x.cols, y.cols)
        )

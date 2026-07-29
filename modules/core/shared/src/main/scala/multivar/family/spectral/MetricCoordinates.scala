package multivar
package family.spectral

import multivar.core.*

import gale.linalg.DMat

/** Whitened matrix together with the metric roots that formed it.
  *
  * Two conventions are supported:
  *
  * - [[MetricCoordinates.Form.Table]]: `matrix = left.half · X · right.half`
  *   (GMD / CPCA table whitening).
  * - [[MetricCoordinates.Form.Cross]]: `matrix = left.pinvHalf · C · right.half`
  *   where `right` is the target *metric* \(R\) (so `right.half = R^{1/2}`).
  *
  * This is shared spectral substrate, not a public estimator.
  */
private[multivar] final class MetricCoordinates private (
    val matrix: DMat,
    val left: MetricRoots,
    val right: MetricRoots,
    val form: MetricCoordinates.Form
):
  /** Map whitened left singular vectors back to original source coordinates. */
  def sourceWeights(u: DMat): DMat =
    left.pinvHalf.applyLeft(u)

  /** Map whitened right singular vectors back to original target coordinates. */
  def targetWeights(v: DMat): DMat =
    form match
      case MetricCoordinates.Form.Table => right.pinvHalf.applyLeft(v)
      case MetricCoordinates.Form.Cross => right.half.applyLeft(v)

  /** Reconstruct an original-space block from a whitened approximation.
    *
    * Table: \(X \approx L^{-1/2} \widehat K R^{-1/2}\).
    * Cross: \(C \approx L^{1/2} \widehat K R^{-1/2}\).
    */
  def decode(block: DMat): DMat =
    form match
      case MetricCoordinates.Form.Table =>
        left.pinvHalf.applyLeft(right.pinvHalf.applyRight(block))
      case MetricCoordinates.Form.Cross =>
        left.half.applyLeft(right.pinvHalf.applyRight(block))

private[multivar] object MetricCoordinates:
  enum Form:
    case Table
    case Cross

  /** Whitened table \(Q^{1/2} X R^{1/2}\). */
  def table(
      input: DMat,
      rowMetric: MetricSpec,
      featureMetric: MetricSpec,
      eigenSolver: SymmetricEigenSolver,
      tolerance: Double,
      policy: StoragePolicy = StoragePolicy.AllowDense,
      rowRole: String = "row metric",
      featureRole: String = "feature metric"
  ): Either[MultivarError, MetricCoordinates] =
    for
      left <- MetricSqrt.factor(rowMetric, eigenSolver, tolerance, policy, rowRole)
      right <- MetricSqrt.factor(featureMetric, eigenSolver, tolerance, policy, featureRole)
      _ <- requireCompatible(input.rows, left.half.dim, "table rows", rowRole)
      _ <- requireCompatible(input.cols, right.half.dim, "table columns", featureRole)
    yield
      val matrix = left.half.applyLeft(right.half.applyRight(input))
      new MetricCoordinates(matrix, left, right, Form.Table)

  /** Whitened table from already-factored roots. */
  def tableFromRoots(
      input: DMat,
      rowRoots: MetricRoots,
      featureRoots: MetricRoots
  ): Either[MultivarError, MetricCoordinates] =
    for
      _ <- requireCompatible(input.rows, rowRoots.half.dim, "table rows", "row metric roots")
      _ <- requireCompatible(input.cols, featureRoots.half.dim, "table columns", "feature metric roots")
    yield
      val matrix = rowRoots.half.applyLeft(featureRoots.half.applyRight(input))
      new MetricCoordinates(matrix, rowRoots, featureRoots, Form.Table)

  /** Canonical cross \(G^{-1/2} C R^{1/2}\) from metric roots of \(G\) and \(R\). */
  def cross(
      cross: DMat,
      sourceRoots: MetricRoots,
      targetMetricRoots: MetricRoots
  ): Either[MultivarError, MetricCoordinates] =
    for
      _ <- requireCompatible(cross.rows, sourceRoots.half.dim, "cross rows", "source geometry")
      _ <- requireCompatible(cross.cols, targetMetricRoots.half.dim, "cross columns", "target metric")
    yield
      val matrix = sourceRoots.pinvHalf.applyLeft(targetMetricRoots.half.applyRight(cross))
      new MetricCoordinates(matrix, sourceRoots, targetMetricRoots, Form.Cross)

  private def requireCompatible(
      observed: Int,
      expected: Int,
      observedRole: String,
      expectedRole: String
  ): Either[MultivarError, Unit] =
    if observed == expected then Right(())
    else
      Left(
        MultivarError.MatrixShapeMismatch(
          s"$observedRole has size $observed but $expectedRole has dimension $expected"
        )
      )

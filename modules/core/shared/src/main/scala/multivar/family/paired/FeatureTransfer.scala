package multivar
package family.paired

import multivar.core.*

import gale.linalg.DMat

/** Right-side feature transfer \(F = Z^{+} Z_b\) for a resolved left block. */
private[paired] object FeatureTransfer:
  /** Moore–Penrose action \(Z^{+} B\) via the supported SVD of \(Z\). */
  def applyPseudoInverse(
      z: DMat,
      block: DMat,
      solver: SvdSolver,
      tolerance: Double
  ): Either[MultivarError, DMat] =
    if z.rows != block.rows || z.cols != block.cols then
      Left(
        MultivarError.MatrixShapeMismatch(
          s"feature transfer expected matching shapes, got Z ${z.rows}x${z.cols} and block ${block.rows}x${block.cols}"
        )
      )
    else if z.rows == 0 || z.cols == 0 || TaskMath.frobeniusNorm2(z) <= 0.0 ||
        TaskMath.frobeniusNorm2(block) <= 0.0
    then Right(DMat.zeros(math.max(z.cols, block.cols), math.max(z.cols, block.cols)))
    else
      val limit = math.min(z.rows, z.cols)
      ComponentCount(limit) match
        case Left(error) =>
          Left(error)
        case Right(request) =>
          solver.decompose(MatrixView.dense(z), request) match
            case Left(error) =>
              Left(error)
            case Right(svd) =>
              val kept = supportedRank(svd.singularValues, tolerance)
              Right:
                if kept == 0 then DMat.zeros(z.cols, z.cols)
                else
                  val u = MatrixOps.takeColumns(svd.u, kept)
                  val v = MatrixOps.takeColumns(svd.v, kept)
                  val inv = invertSingular(svd.singularValues, kept)
                  // Z^+ B = V Σ^+ U^T B
                  val utBlock = GaleNumerics.multiply(u.transpose, block)
                  val scaled = MetricOperator.scaleColumnsDense(utBlock.transpose, inv).transpose
                  GaleNumerics.multiply(v, scaled)

  def encoder(transfer: DMat, rightVectors: DMat): DMat =
    if rightVectors.cols == 0 then DMat.zeros(transfer.rows, 0)
    else GaleNumerics.multiply(transfer, rightVectors)

  def canonicalScores(svd: SvdResult): DMat =
    if svd.singularValues.length == 0 then DMat.zeros(svd.u.rows, 0)
    else MetricOperator.scaleColumnsDense(svd.u, svd.singularValues)

  private def supportedRank(values: gale.linalg.DVec, tolerance: Double): Int =
    var largest = 0.0
    var i = 0
    while i < values.length do
      val absolute = math.abs(values(i))
      if absolute > largest then largest = absolute
      i += 1
    val cutoff = tolerance * math.max(1.0, largest)
    var kept = 0
    i = 0
    while i < values.length do
      if math.abs(values(i)) > cutoff then kept += 1
      i += 1
    kept

  private def invertSingular(values: gale.linalg.DVec, kept: Int): gale.linalg.DVec =
    val out = new Array[Double](kept)
    var i = 0
    while i < kept do
      out(i) = 1.0 / values(i)
      i += 1
    GaleNumerics.vectorFromArray(out)

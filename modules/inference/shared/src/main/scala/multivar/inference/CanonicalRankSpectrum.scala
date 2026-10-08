package multivar.inference

import gale.backend.Backend.given
import gale.linalg.{DMat, QROptions, QRPivoting}
import gale.spectral.SpectralBackend.given
import gale.spectral.SpectralDiagnostics
import multivar.core.MatrixOps

/** All classical canonical roots in already prepared, zero-mean coordinates.
  * No centering, root truncation, or choice of a deflated reference law. */
final class CanonicalRankSpectrum private[inference] (
    val rows: Int, val leftColumns: Int, val rightColumns: Int,
    val correlations: Vector[Double], val observedWilks: Vector[Double],
    val observedCanonicalDiagnostics: SpectralDiagnostics, val plannedElements: Long
):
  def candidateRank: Int = correlations.size
  val qrRankTolerance: Double = CanonicalRankSpectrum.qrRankTolerance

private[inference] final case class PreparedCanonicalRank(
    spectrum: CanonicalRankSpectrum, leftBasis: DMat, rightBasis: DMat,
    leftDirections: DMat, rightDirections: DMat
)

object CanonicalRankSpectrum:
  val qrRankTolerance: Double = 1e-12

  def from(left: DMat, right: DMat, maximumElements: Long = 4_000_000L): Either[InferenceError, CanonicalRankSpectrum] =
    prepare(left, right, maximumElements).map(_.spectrum)

  private[inference] def prepare(left: DMat, right: DMat, maximumElements: Long): Either[InferenceError, PreparedCanonicalRank] =
    val n = BigInt(left.rows); val p = BigInt(left.cols); val q = BigInt(right.cols)
    val required = 12 * n * (p + q) + 40 * (p * p + q * q + p * q)
    if left.rows != right.rows then Left(InferenceError.RowCountMismatch("canonical paired coordinates", left.rows, right.rows))
    else if left.cols < 1 || right.cols < 1 || n <= p + q then Left(InferenceError.InvalidPartition("classical rank confirmation requires rows > left columns + right columns"))
    else if maximumElements < 0L || required > maximumElements || Vector(n * p, n * q, p * p, q * q, p * q).exists(_ > Int.MaxValue) then
      Left(InferenceError.UnsupportedProblem(s"canonical spectrum requires $required owned elements, allowed $maximumElements"))
    else
      for
        _ <- MatrixOps.checkFinite("left canonical input", left).left.map(e => InferenceError.NumericalFailure("canonical input", e.message))
        _ <- MatrixOps.checkFinite("right canonical input", right).left.map(e => InferenceError.NumericalFailure("canonical input", e.message))
        qx <- orthonormal(left)
        qy <- orthonormal(right)
        fit <- (qx.t * qy).svd.flatMap(_.requireConverged).left.map(e => InferenceError.NumericalFailure("complete canonical SVD", e.toString))
        roots = fit.singularValues.toVector
        _ <- if roots.size == math.min(left.cols, right.cols) then Right(()) else Left(InferenceError.InvalidSpectrum("complete CCA omitted candidate roots"))
        _ <- wilks(roots)
        tails = Vector.tabulate(roots.size)(k => -roots.drop(k).map(r => math.log1p(-r * r)).sum)
      yield PreparedCanonicalRank(new CanonicalRankSpectrum(left.rows, left.cols, right.cols, roots, tails, fit.diagnostics, required.toLong), qx, qy, fit.u, fit.vt.t)

  private[inference] def orthonormal(input: DMat): Either[InferenceError, DMat] =
    val qr = input.qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(qrRankTolerance)))
    val rank = qr.diagnostics.rank.getOrElse(0)
    if rank != input.cols then Left(InferenceError.RankLoss(input.cols, rank))
    else
      val selectors = DMat.tabulate(input.rows, input.cols)((i, j) => if i == j then 1.0 else 0.0)
      qr.applyQ(selectors).left.map(e => InferenceError.NumericalFailure("canonical QR basis", e.toString))

  private[inference] def wilks(values: Vector[Double]): Either[InferenceError, Double] =
    if values.isEmpty || values.exists(x => !x.isFinite || x < 0.0 || x >= 1.0 || 1.0 - x * x <= 1e-12) then
      Left(InferenceError.InvalidSpectrum("classical canonical correlations require finite [0,1) values with 1-rho^2 > 1e-12; unit or numerically unit correlations are unavailable"))
    else
      val result = -values.map(x => math.log1p(-x * x)).sum
      if result.isFinite then Right(result) else Left(InferenceError.NonFiniteStatistic("canonical Wilks", result))

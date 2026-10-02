package multivar.inference

import gale.backend.Backend.given
import gale.linalg.{DMat, QROptions, QRPivoting}
import gale.spectral.SpectralBackend.given
import multivar.core.{LinalgErrorAdapter, MatrixOps}

enum CanonicalResidualMethod:
  case HuhJhun
  /** Source row selection must be frozen outside this numerical adapter. */
  case Theil(selectedSourceRows: Vector[Int])

/** Geometry of nuisance residual coordinates, not an exchangeability proof.
  * The caller supplies the error law and legal actions separately. */
final class CanonicalResidualBasis private (
    val nuisanceRank: Int, val matrix: DMat, val method: CanonicalResidualMethod,
    val plannedElements: Long
):
  val qrRankTolerance: Double = 1e-12
  val lawTolerance: Double = 1e-8
  def project(input: DMat, maximumOutputElements: Long = 4_000_000L): Either[InferenceError, DMat] =
    if input.rows != matrix.rows then Left(InferenceError.RowCountMismatch("canonical residual projection", matrix.rows, input.rows))
    else if maximumOutputElements < 0L || BigInt(matrix.cols) * input.cols > maximumOutputElements || BigInt(matrix.cols) * input.cols > Int.MaxValue then Left(InferenceError.UnsupportedProblem("canonical residual output cardinality exceeds its explicit ceiling"))
    else
      MatrixOps.checkFinite("canonical residual input", input).left.map(error => InferenceError.NumericalFailure("residual projection", error.message)).flatMap: _ =>
        val output = matrix.t * input
        MatrixOps.checkFinite("canonical residual output", output).left.map(error => InferenceError.NumericalFailure("residual projection", error.message)).map(_ => output)

object CanonicalResidualBasis:
  /** The design must include the actual centering column in its working
    * coordinates. Rank is measured, including a redundant supplied intercept.
    * Owned-array bounds exclude private Gale kernel workspace. */
  def from(nuisance: DMat, method: CanonicalResidualMethod = CanonicalResidualMethod.HuhJhun,
      maximumElements: Long = 4_000_000L): Either[InferenceError, CanonicalResidualBasis] =
    val n = BigInt(nuisance.rows); val z = BigInt(nuisance.cols)
    val required = 12 * n * n + 12 * n * z + 16 * z * z
    if nuisance.rows < 3 || nuisance.cols < 1 || nuisance.cols >= nuisance.rows then Left(InferenceError.InvalidPartition("nuisance design requires positive columns below its row count"))
    else if maximumElements < 0L || required > maximumElements || n * n > Int.MaxValue || n * z > Int.MaxValue then Left(InferenceError.UnsupportedProblem(s"residual basis requires $required owned elements, allowed $maximumElements"))
    else
      for
        _ <- adapt("nuisance", MatrixOps.checkFinite("canonical nuisance", nuisance))
        qr = nuisance.qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(1e-12)))
        rank = qr.diagnostics.rank.getOrElse(0)
        _ <- if rank > 0 && rank < nuisance.rows then Right(()) else Left(InferenceError.RankLoss(nuisance.cols, rank))
        count = nuisance.rows - rank
        basis <- method match
          case CanonicalResidualMethod.HuhJhun =>
            val selectors = DMat.tabulate(nuisance.rows, count)((i, j) => if i == j + rank then 1.0 else 0.0)
            linear("Huh-Jhun orthogonal complement", qr.applyQ(selectors))
          case CanonicalResidualMethod.Theil(rows) =>
            if rows.size != count || rows.distinct.size != rows.size || rows.exists(i => i < 0 || i >= nuisance.rows) then
              Left(InferenceError.InvalidPartition(s"Theil selection must contain $count distinct actual source rows"))
            else
              val selectors = DMat.tabulate(nuisance.rows, count)((i, j) => if i == rows(j) then 1.0 else 0.0)
              for
                residual <- linear("Theil residual selection", qr.residualize(selectors))
                gram = selectors.t * residual
                _ <- linear("Theil selected residual SPD", gram.cholesky).map(_ => ())
                spectral <- linear("Theil selected residual SVD", gram.svd.flatMap(_.requireConverged))
                _ <- if spectral.rank == count && (0 until count).forall(i => spectral.singularValues(i).isFinite && spectral.singularValues(i) > 0.0) then Right(()) else Left(InferenceError.RankLoss(count, spectral.rank))
                inverseHalf = spectral.u * DMat.tabulate(count, count)((i, j) => if i == j then 1.0 / math.sqrt(spectral.singularValues(i)) else 0.0) * spectral.u.t
              yield residual * inverseHalf
        _ <- adapt("residual basis", MatrixOps.checkFinite("canonical residual basis", basis))
        orthogonal = basis.t * basis
        _ <- if (0 until count).forall(i => (0 until count).forall(j => math.abs(orthogonal(i, j) - (if i == j then 1.0 else 0.0)) <= 1e-8)) then Right(()) else Left(InferenceError.NumericalFailure("residual basis", "orthonormality residual exceeds 1e-8"))
        orthogonality = basis.t * nuisance
        _ <- adapt("nuisance orthogonality", MatrixOps.checkFinite("canonical nuisance orthogonality", orthogonality))
        _ <- if (0 until nuisance.cols).forall: j =>
          var scale = 0.0
          var i = 0
          while i < nuisance.rows do
            scale = math.max(scale, math.abs(nuisance(i, j)))
            i += 1
          (0 until count).forall(i => math.abs(orthogonality(i, j)) / math.max(1.0, scale) / nuisance.rows <= 1e-8)
        then Right(()) else Left(InferenceError.NumericalFailure("residual basis", "nuisance orthogonality residual exceeds scaled tolerance"))
      yield new CanonicalResidualBasis(rank, basis, method, required.toLong)

  private def adapt[A](role: String, result: Either[multivar.core.MultivarError, A]): Either[InferenceError, A] =
    result.left.map(error => InferenceError.NumericalFailure(role, error.message))

  private def linear[A](role: String, result: Either[gale.linalg.LinAlgError, A]): Either[InferenceError, A] =
    adapt(role, LinalgErrorAdapter.adapt(result))

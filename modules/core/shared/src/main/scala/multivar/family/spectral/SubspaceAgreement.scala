package multivar
package family.spectral

import gale.backend.Backend.given
import gale.linalg.DMat
import gale.spectral.SpectralBackend.given
import multivar.core.*

/** Subspace angles and fixed-coordinate agreement answer different questions.
  * Absolute axis cosines ignore signs but preserve declared column order;
  * they do not search a permutation or choose a favorable rotation. */
final case class SubspaceAgreement(
    leftRank: Int,
    rightRank: Int,
    principalAngles: Vector[Double],
    fixedAxisAbsoluteCosines: Option[Vector[Double]],
    plannedElements: Long
):
  require(leftRank > 0 && rightRank > 0 && principalAngles.size == math.min(leftRank, rightRank))
  require(principalAngles.forall(a => a.isFinite && a >= 0.0 && a <= math.Pi / 2.0))
  require(fixedAxisAbsoluteCosines.forall(xs => leftRank == rightRank && xs.size == leftRank && xs.forall(x => x.isFinite && x >= 0.0 && x <= 1.0)))
  require(plannedElements >= 0L)


object SubspaceAgreement:
  /** Reuses the existing spectral range extractor and Gale's converged compact
    * cross-basis SVD. No observation-square projector is constructed. The
    * conservative owned-array bound excludes private kernel scratch. */
  def compare(left: DMat, right: DMat, rankTolerance: Double = 1e-10,
      maximumElements: Long = 4_000_000L): Either[MultivarError, SubspaceAgreement] =
    val n = BigInt(left.rows); val p = BigInt(left.cols); val q = BigInt(right.cols)
    val required = 10 * n * (p + q) + 20 * (p * p + q * q + p * q)
    if left.rows != right.rows || left.rows <= 0 || left.cols <= 0 || right.cols <= 0 || left.cols > left.rows || right.cols > right.rows then
      Left(MultivarError.MatrixShapeMismatch("subspace generators require the same positive ambient dimension and narrow positive columns"))
    else if !rankTolerance.isFinite || rankTolerance <= 0.0 || rankTolerance > 1e-4 then Left(MultivarError.InvalidTolerance("subspace rank", rankTolerance))
    else if maximumElements < 0L || required > maximumElements || Vector(n * p, n * q, p * q, p * p, q * q).exists(_ > Int.MaxValue) then
      Left(MultivarError.InvalidMap(s"subspace owned-array admission requires $required elements, allowed $maximumElements"))
    else
      for
        _ <- MatrixOps.checkFinite("left subspace", left)
        _ <- MatrixOps.checkFinite("right subspace", right)
        leftNorms = Vector.tabulate(left.cols)(i => left.col(i).norm2)
        rightNorms = Vector.tabulate(right.cols)(i => right.col(i).norm2)
        _ <- if (leftNorms ++ rightNorms).forall(x => x.isFinite && x > 0.0) then Right(()) else Left(MultivarError.InvalidMap("subspace generators require nonzero finite column norms"))
        normalizedLeft = DMat.tabulate(left.rows, left.cols)((i, j) => left(i, j) / leftNorms(j))
        normalizedRight = DMat.tabulate(right.rows, right.cols)((i, j) => right(i, j) / rightNorms(j))
        a <- OrthonormalSubspace.span(normalizedLeft, DenseSolvers.svd, rankTolerance)
        b <- OrthonormalSubspace.span(normalizedRight, DenseSolvers.svd, rankTolerance)
        _ <- if a.rank == left.cols && b.rank == right.cols then Right(()) else Left(MultivarError.SolverFailed(s"subspace generators are numerically rank deficient: ${a.rank}/${left.cols}, ${b.rank}/${right.cols}"))
        decomposition <- LinalgErrorAdapter.adapt((a.basis.t * b.basis).svd.flatMap(_.requireConverged))
        singular = decomposition.singularValues
        _ <- if (0 until singular.length).forall(i => singular(i).isFinite && singular(i) >= 0.0 && singular(i) <= 1.0 + 1e-8) then Right(()) else Left(MultivarError.SolverFailed("cross-subspace singular values outside [0,1]"))
        angles = Vector.tabulate(singular.length)(i => math.acos(math.min(1.0, singular(i))))
        axes <-
          if left.cols != right.cols then Right(None)
          else
            val values = Vector.newBuilder[Double]
            var column = 0
            var invalid = false
            while column < left.cols && !invalid do
              var cosine = 0.0; var row = 0
              while row < left.rows do
                cosine += normalizedLeft(row, column) * normalizedRight(row, column)
                row += 1
              if !cosine.isFinite || math.abs(cosine) > 1.0 + 1e-8 then invalid = true
              else values += math.min(1.0, math.abs(cosine))
              column += 1
            if invalid then Left(MultivarError.SolverFailed("fixed-axis normalized cosine")) else Right(Some(values.result()))
      yield SubspaceAgreement(a.rank, b.rank, angles, axes, required.toLong)

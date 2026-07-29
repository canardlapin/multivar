package multivar.inference

import gale.linalg.DMat
import gale.linalg.DVec
import multivar.core.{ComponentCount, DenseSolvers, MatrixView, SvdResult, SvdSolver}
import resample4s.kernel.Permutation
import resample4s.spi.AlgorithmId

final case class PairedMatrixData private (
    x: DMat,
    y: DMat
)

object PairedMatrixData:
  def from(
      x: DMat,
      y: DMat
  ): Either[InferenceError, PairedMatrixData] =
    if x.rows != y.rows then
      Left(InferenceError.RowCountMismatch("paired matrix data", x.rows, y.rows))
    else if x.rows < 2 then Left(InferenceError.InvalidCount("paired matrix rows", x.rows))
    else if x.cols < 1 then Left(InferenceError.InvalidCount("paired X columns", x.cols))
    else if y.cols < 1 then Left(InferenceError.InvalidCount("paired Y columns", y.cols))
    else
      firstNonFinite(x, "paired X")
        .orElse(firstNonFinite(y, "paired Y"))
        .toLeft(PairedMatrixData(x, y))

  private def firstNonFinite(
      matrix: DMat,
      role: String
  ): Option[InferenceError] =
    val values = matrix.copyData
    var i = 0
    while i < values.length do
      if !values(i).isFinite then
        return Some(InferenceError.NonFiniteStatistic(s"$role entry $i", values(i)))
      i += 1
    None

final case class CrossCovarianceFit private (
    roots: DVec
)

object CrossCovarianceFit:
  private[inference] def unsafe(values: Vector[Double]): CrossCovarianceFit =
    CrossCovarianceFit(InferenceNumerics.vectorFromSeq(values))

final case class PairedCrossCore private[inference] (
    xRows: DMat,
    xValues: DVec,
    yRows: DMat,
    yValues: DVec
)

final case class ReducedCrossFit private[inference] (
    roots: Vector[Double]
)

final case class ExactCrossCovarianceReduction(
    solver: SvdSolver = DenseSolvers.svd
) extends ExactRefitReduction[PairedMatrixData, CrossCovarianceFit, Permutation]:
  override type Core = PairedCrossCore
  override type ReducedFit = ReducedCrossFit

  override val algorithm: AlgorithmId =
    AlgorithmId
      .of("paired-thin-svd-cross/v1")
      .fold(
        error => throw IllegalStateException(error.message),
        (value: AlgorithmId) => value
      )

  override def core(
      data: PairedMatrixData,
      observed: CrossCovarianceFit
  ): Either[InferenceError, PairedCrossCore] =
    for
      x <- decompose(data.x)
      y <- decompose(data.y)
    yield PairedCrossCore(x.u, x.singularValues, y.u, y.singularValues)

  override def update(
      core: PairedCrossCore,
      action: Permutation
  ): Either[InferenceError, ReducedCrossFit] =
    for
      permutedY <- applyPermutation(core.yRows, action)
      overlap = InferenceNumerics.transposeMultiply(core.xRows, permutedY)
      reduced = scaleCross(overlap, core.xValues, core.yValues)
      fit <- roots(reduced)
    yield ReducedCrossFit(fit)

  override def lift(reduced: ReducedCrossFit): Either[InferenceError, CrossCovarianceFit] =
    Right(CrossCovarianceFit.unsafe(reduced.roots))

  private def decompose(matrix: DMat): Either[InferenceError, SvdResult] =
    val rank = Math.min(matrix.rows, matrix.cols)
    ComponentCount(rank)
      .left.map(error =>
        InferenceError.NumericalFailure("exact reduction component count", error.message)
      )
      .flatMap(requested =>
        solver.decompose(MatrixView.dense(matrix), requested)
          .left.map(error => InferenceError.NumericalFailure("exact reduction SVD", error.message))
      )

  private def roots(matrix: DMat): Either[InferenceError, Vector[Double]] =
    decompose(matrix).map { fit =>
      Vector.tabulate(fit.singularValues.length) { index =>
        val value = fit.singularValues(index)
        value * value
      }
    }

  private def scaleCross(
      overlap: DMat,
      xValues: DVec,
      yValues: DVec
  ): DMat =
    val out = new Array[Double](overlap.rows * overlap.cols)
    var row = 0
    while row < overlap.rows do
      var col = 0
      while col < overlap.cols do
        out(row * overlap.cols + col) =
          xValues(row) * overlap(row, col) * yValues(col)
        col += 1
      row += 1
    InferenceNumerics.matrixFromRowMajor(overlap.rows, overlap.cols, out)

object CrossCovarianceRefit:
  def fit(
      data: PairedMatrixData,
      solver: SvdSolver = DenseSolvers.svd
  ): Either[InferenceError, CrossCovarianceFit] =
    val cross = InferenceNumerics.transposeMultiply(data.x, data.y)
    val rank = Math.min(cross.rows, cross.cols)
    ComponentCount(rank)
      .left.map(error =>
        InferenceError.NumericalFailure("cross-covariance component count", error.message)
      )
      .flatMap(requested =>
        solver.decompose(MatrixView.dense(cross), requested)
          .left.map(error => InferenceError.NumericalFailure("cross-covariance SVD", error.message))
      )
      .map { fit =>
        CrossCovarianceFit.unsafe(Vector.tabulate(fit.singularValues.length) { index =>
          val value = fit.singularValues(index)
          value * value
        })
      }

  given refit: Refit[PairedMatrixData, CrossCovarianceFit] with
    override def fit(data: PairedMatrixData): Either[InferenceError, CrossCovarianceFit] =
      CrossCovarianceRefit.fit(data)

  given rowPermutationAction: DataAction[PairedMatrixData, Permutation] with
    override def apply(
        data: PairedMatrixData,
        action: Permutation
    ): Either[InferenceError, PairedMatrixData] =
      applyPermutation(data.y, action)
        .flatMap(PairedMatrixData.from(data.x, _))

  given exactReduction: ExactCrossCovarianceReduction =
    ExactCrossCovarianceReduction()

private def applyPermutation(
    matrix: DMat,
    permutation: Permutation
): Either[InferenceError, DMat] =
  if permutation.domain != matrix.rows then
    Left(
      InferenceError.RowCountMismatch(
        "row permutation input",
        permutation.domain,
        matrix.rows
      )
    )
  else
    Right(
      matrix.selectRows(permutation.toIArray.toIndexedSeq)
    )

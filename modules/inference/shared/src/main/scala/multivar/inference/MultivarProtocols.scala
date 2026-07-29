package multivar.inference

import gale.linalg.DMat
import gale.spectral.{SingularSelection, Svds}
import multivar.core.{ComponentCount, DenseSolvers, MatrixView, SvdResult, SvdSolver}
import resample4s.kernel.{Permutation, Seed}

final case class PcaVarianceState private[inference] (
    residual: DMat,
    removed: Int
)

object PcaVarianceState:
  def from(input: DMat): Either[InferenceError, PcaVarianceState] =
    ProtocolMatrices.validate("PCA input", input).map { _ =>
      PcaVarianceState(ProtocolMatrices.centerColumns(input), removed = 0)
    }

final case class PcaVarianceProtocol(
    solver: SvdSolver = DenseSolvers.svd
) extends ExactLadderProtocol[
      PcaVarianceState,
      PcaFitFamily,
      TargetKind.VarianceRoots,
      NullKind.RowPermutation
    ]:
  type Action = Vector[Permutation]

  override val fit: FitDescriptor[PcaFitFamily] = FitDescriptor.Pca
  override val target: TargetSpec[TargetKind.VarianceRoots] = TargetSpec.VarianceRoots
  override val nullHypothesis: NullSpec[NullKind.RowPermutation] =
    NullSpec.PermuteRows
  override val validity: ValidityClaim = ValidityClaim.Exact

  override def roots(initial: PcaVarianceState): Either[InferenceError, Vector[Double]] =
    ProtocolMatrices.fullSingularValues(initial.residual)
      .map(_.map(value => value * value))

  override def observed(state: PcaVarianceState): Either[InferenceError, Double] =
    val total = ProtocolMatrices.frobeniusSquared(state.residual)
    if total <= ProtocolMatrices.zeroTolerance(total) then Right(0.0)
    else ProtocolMatrices.singularValues(state.residual, 1, solver).map { values =>
      values.head * values.head / total
    }

  override def action(
      state: PcaVarianceState,
      replicate: ReplicateId,
      seed: Seed
  ): Either[InferenceError, Action] =
    ResamplingPlans.pcaColumns(
      seed,
      replicate,
      state.residual.rows,
      state.residual.cols
    )

  override def actionForDesign(
      state: PcaVarianceState,
      design: ResamplingDesign[?],
      replicate: ReplicateId,
      seed: Seed
  ): Either[InferenceError, Action] =
    for
      permutation <- ProtocolDesign.permutation(
        design,
        state.residual.rows,
        "PCA resampling design"
      )
      actions <- ResamplingPlans.independent(
        permutation,
        seed,
        replicate,
        state.residual.cols,
        ResamplingPlans.PcaColumnDomain
      )
    yield actions

  override def validateDesign(
      state: PcaVarianceState,
      design: ResamplingDesign[?]
  ): Either[InferenceError, Unit] =
    ProtocolDesign
      .permutation(design, state.residual.rows, "PCA resampling design")
      .map(_ => ())

  override def nullStatistic(
      state: PcaVarianceState,
      step: ComponentIx,
      replicate: ReplicateId,
      action: Action
  ): Either[InferenceError, Double] =
    for
      permuted <- ProtocolMatrices.permuteColumns(state.residual, action)
      values <-
        if ProtocolMatrices.rankLimit(permuted) <= 0 then Right(Vector.empty[Double])
        else ProtocolMatrices.singularValues(permuted, 1, solver)
    yield
      if values.isEmpty then 0.0
      else
        val total = ProtocolMatrices.frobeniusSquared(state.residual)
        if total <= ProtocolMatrices.zeroTolerance(total) then 0.0
        else values.head * values.head / total

  override def remove(state: PcaVarianceState): Either[InferenceError, PcaVarianceState] =
    ProtocolMatrices.svd(state.residual, 1, solver).map { fit =>
      val residual = ProtocolMatrices.removeRankOne(
        state.residual,
        fit.u,
        fit.singularValues(0),
        fit.v
      )
      PcaVarianceState(ProtocolMatrices.zeroSmall(residual), state.removed + 1)
    }

final case class PlscCovarianceState private[inference] (
    x: DMat,
    y: DMat,
    removed: Int
)

object PlscCovarianceState:
  def from(x: DMat, y: DMat): Either[InferenceError, PlscCovarianceState] =
    for
      _ <- ProtocolMatrices.validate("PLSC X", x)
      _ <- ProtocolMatrices.validate("PLSC Y", y)
      _ <-
        if x.rows == y.rows then Right(())
        else Left(InferenceError.RowCountMismatch("PLSC paired blocks", x.rows, y.rows))
    yield PlscCovarianceState(
      ProtocolMatrices.centerColumns(x),
      ProtocolMatrices.centerColumns(y),
      removed = 0
    )

final case class PlscCovarianceProtocol(
    override val nullHypothesis: NullSpec[NullKind.PairedIndependence] =
      NullSpec.BreakY,
    solver: SvdSolver = DenseSolvers.svd
) extends ExactLadderProtocol[
      PlscCovarianceState,
      PlscFitFamily,
      TargetKind.CovarianceRoots,
      NullKind.PairedIndependence
    ]:
  type Action = Permutation

  override val fit: FitDescriptor[PlscFitFamily] = FitDescriptor.Plsc
  override val target: TargetSpec[TargetKind.CovarianceRoots] = TargetSpec.CovarianceRoots
  override val validity: ValidityClaim = ValidityClaim.Exact

  override def roots(initial: PlscCovarianceState): Either[InferenceError, Vector[Double]] =
    val cross = InferenceNumerics.transposeMultiply(initial.x, initial.y)
    ProtocolMatrices.fullSingularValues(cross)
      .map(_.map(value => value * value))

  override def observed(state: PlscCovarianceState): Either[InferenceError, Double] =
    leadingCrossRoot(state.x, state.y)

  override def action(
      state: PlscCovarianceState,
      replicate: ReplicateId,
      seed: Seed
  ): Either[InferenceError, Action] =
    ResamplingPlans.permutation(seed, replicate, state.y.rows)

  override def actionForDesign(
      state: PlscCovarianceState,
      design: ResamplingDesign[?],
      replicate: ReplicateId,
      seed: Seed
  ): Either[InferenceError, Action] =
    ProtocolDesign
      .permutation(design, state.y.rows, "PLSC resampling design")
      .flatMap(_.draw(seed, replicate))

  override def validateDesign(
      state: PlscCovarianceState,
      design: ResamplingDesign[?]
  ): Either[InferenceError, Unit] =
    ProtocolDesign
      .permutation(design, state.y.rows, "PLSC resampling design")
      .map(_ => ())

  override def nullStatistic(
      state: PlscCovarianceState,
      step: ComponentIx,
      replicate: ReplicateId,
      action: Action
  ): Either[InferenceError, Double] =
    nullHypothesis match
      case NullSpec.BreakX =>
        leadingCrossRoot(
          state.x.selectRows(action.toIArray.toIndexedSeq),
          state.y
        )
      case NullSpec.BreakY =>
        leadingCrossRoot(
          state.x,
          state.y.selectRows(action.toIArray.toIndexedSeq)
        )

  override def remove(state: PlscCovarianceState): Either[InferenceError, PlscCovarianceState] =
    val cross = InferenceNumerics.transposeMultiply(state.x, state.y)
    ProtocolMatrices.svd(cross, 1, solver).map { fit =>
      val xNext = ProtocolMatrices.removeFeatureDirection(state.x, fit.u)
      val yNext = ProtocolMatrices.removeFeatureDirection(state.y, fit.v)
      PlscCovarianceState(
        ProtocolMatrices.zeroSmall(xNext),
        ProtocolMatrices.zeroSmall(yNext),
        state.removed + 1
      )
    }

  private def leadingCrossRoot(x: DMat, y: DMat): Either[InferenceError, Double] =
    ProtocolMatrices.singularValues(InferenceNumerics.transposeMultiply(x, y), 1, solver).map { values =>
      values.headOption.fold(0.0)(value => value * value)
    }

private object ProtocolMatrices:
  def validate(role: String, matrix: DMat): Either[InferenceError, Unit] =
    if matrix.rows < 2 then Left(InferenceError.InvalidCount(s"$role rows", matrix.rows))
    else if matrix.cols < 1 then Left(InferenceError.InvalidCount(s"$role columns", matrix.cols))
    else
      val values = matrix.copyData
      var i = 0
      while i < values.length do
        if !values(i).isFinite then return Left(InferenceError.NonFiniteStatistic(s"$role entry $i", values(i)))
        i += 1
      Right(())

  def centerColumns(matrix: DMat): DMat =
    val out = matrix.copyData
    var col = 0
    while col < matrix.cols do
      var mean = 0.0
      var row = 0
      while row < matrix.rows do
        mean += matrix(row, col)
        row += 1
      mean /= matrix.rows
      row = 0
      while row < matrix.rows do
        out(row * matrix.cols + col) -= mean
        row += 1
      col += 1
    InferenceNumerics.matrixFromRowMajor(matrix.rows, matrix.cols, out)

  def rankLimit(matrix: DMat): Int =
    Math.min(matrix.rows, matrix.cols)

  def singularValues(
      matrix: DMat,
      rank: Int,
      solver: SvdSolver
  ): Either[InferenceError, Vector[Double]] =
    svd(matrix, rank, solver).map(_.singularValues.toVector)

  def fullSingularValues(
      matrix: DMat
  ): Either[InferenceError, Vector[Double]] =
    Svds
      .svd(matrix, SingularSelection.All)
      .left
      .map(error =>
        InferenceError.NumericalFailure(
          "full-spectrum SVD",
          error.getMessage
        )
      )
      .map(_.singularValues.toVector)

  def svd(
      matrix: DMat,
      rank: Int,
      solver: SvdSolver
  ): Either[InferenceError, SvdResult] =
    ComponentCount(rank)
      .left.map(error => InferenceError.NumericalFailure("SVD component count", error.message))
      .flatMap { requested =>
        solver.decompose(MatrixView.dense(matrix), requested)
          .left.map(error => InferenceError.NumericalFailure("dense SVD", error.message))
      }

  def frobeniusSquared(matrix: DMat): Double =
    val values = matrix.copyData
    var total = 0.0
    var i = 0
    while i < values.length do
      total += values(i) * values(i)
      i += 1
    total

  def zeroTolerance(scale: Double): Double =
    Math.max(1.0, scale) * Math.ulp(1.0)

  def zeroSmall(matrix: DMat): DMat =
    val scale = frobeniusSquared(matrix)
    if scale <= zeroTolerance(scale) then DMat.zeros(matrix.rows, matrix.cols)
    else matrix

  def permuteColumns(
      matrix: DMat,
      permutations: Vector[Permutation]
  ): Either[InferenceError, DMat] =
    if permutations.length != matrix.cols then
      Left(
        InferenceError.RowCountMismatch(
          "PCA column permutation lanes",
          matrix.cols,
          permutations.length
        )
      )
    else
      val out = new Array[Double](matrix.rows * matrix.cols)
      var col = 0
      while col < matrix.cols do
        val permutation = permutations(col).toIArray
        var row = 0
        while row < matrix.rows do
          out(row * matrix.cols + col) = matrix(permutation(row), col)
          row += 1
        col += 1
      Right(
        InferenceNumerics.matrixFromRowMajor(
          matrix.rows,
          matrix.cols,
          out
        )
      )

  def removeRankOne(
      matrix: DMat,
      u: DMat,
      singularValue: Double,
      v: DMat
  ): DMat =
    val out = matrix.copyData
    var row = 0
    while row < matrix.rows do
      var col = 0
      while col < matrix.cols do
        out(row * matrix.cols + col) -= singularValue * u(row, 0) * v(col, 0)
        col += 1
      row += 1
    InferenceNumerics.matrixFromRowMajor(matrix.rows, matrix.cols, out)

  def removeFeatureDirection(matrix: DMat, direction: DMat): DMat =
    val scores = InferenceNumerics.multiply(matrix, direction)
    val out = matrix.copyData
    var row = 0
    while row < matrix.rows do
      var col = 0
      while col < matrix.cols do
        out(row * matrix.cols + col) -= scores(row, 0) * direction(col, 0)
        col += 1
      row += 1
    InferenceNumerics.matrixFromRowMajor(matrix.rows, matrix.cols, out)

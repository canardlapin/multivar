package multivar
package family.paired

import multivar.core.*
import multivar.family.spectral.MetricCoordinates
import multivar.family.spectral.OrthonormalSubspace

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.Matrix

/** Canonical least-squares task operator in Euclidean working coordinates.
  *
  * Predictive form: \(K = G_\rho^{-1/2} C R^{1/2}\) ([[MetricCoordinates.Form.Cross]]).
  * Decomposition form: \(K = Q^{1/2} X R^{1/2}\) ([[MetricCoordinates.Form.Table]]).
  */
private[paired] final class CanonicalTask private (
    val coordinates: MetricCoordinates,
    val sourceGram: DMat,
    val cross: DMat,
    val zeroObjective: Double,
    val sourceAxis: TaskAxis,
    val targetAxis: TaskAxis
):
  def matrix: DMat =
    coordinates.matrix

  def availableGain: Double =
    TaskMath.frobeniusNorm2(matrix)

  def fullLoss: Double =
    math.max(0.0, zeroObjective - availableGain)

  /** Working coefficient map for a predictive (cross) task. */
  def workingMap(svd: SvdResult): DMat =
    if svd.singularValues.length == 0 then
      DMat.zeros(coordinates.left.half.dim, coordinates.right.half.dim)
    else
      val uScaled =
        MetricOperator.scaleColumnsDense(coordinates.sourceWeights(svd.u), svd.singularValues)
      val vDecoded = coordinates.right.pinvHalf.applyLeft(svd.v)
      GaleNumerics.multiply(uScaled, vDecoded.transpose)

  /** Whitened block reconstruction \(\widehat K = U\Sigma V^\top\). */
  def whitenedReconstruction(svd: SvdResult): DMat =
    TaskMath.reconstruction(svd)

  /** Original-axis reconstruction of a whitened block (table or cross). */
  def decodeBlock(block: DMat): DMat =
    coordinates.decode(block)

  def sourceScores(workingX: DMat, svd: SvdResult): DMat =
    if svd.singularValues.length == 0 then DMat.zeros(workingX.rows, 0)
    else
      coordinates.form match
        case MetricCoordinates.Form.Cross =>
          val weights =
            MetricOperator.scaleColumnsDense(coordinates.sourceWeights(svd.u), svd.singularValues)
          GaleNumerics.multiply(workingX, weights)
        case MetricCoordinates.Form.Table =>
          MetricOperator.scaleColumnsDense(coordinates.sourceWeights(svd.u), svd.singularValues)

private[paired] final case class TaskAxis(label: String, dim: Int)

private[paired] object CanonicalTask:
  def fromPrepared(
      prepared: PreparedPair,
      regularization: RegressionRegularization,
      targetMetric: MetricSpec,
      eigenSolver: SymmetricEigenSolver,
      tolerance: Double,
      policy: StoragePolicy
  ): Either[MultivarError, CanonicalTask] =
    for
      preparedProblem <- PairedOperatorProblem.fromMatrices(
        prepared.xWorking,
        prepared.yWorking,
        prepared.rowMetric,
        "task-components",
        policy
      )
      problem = preparedProblem.value
      sourceMarginal <- taskSemantic(problem.sourceMarginal.toDense)
      crossDense <- taskSemantic(problem.cross.toDense)
      targetMarginal <- taskSemantic(problem.targetMarginal.toDense)
      ridgeScale = Math.max(1, prepared.moments.sampleCount - 1).toDouble
      ridge = pairedRidgeValue(regularization) * ridgeScale
      penalized = MatrixOps.addRidge(sourceMarginal, ridge)
      sourceRoots <- TaskSpdRoots.factor(penalized, eigenSolver, tolerance, "task source geometry")
      targetDense <- targetMetric.toDense(policy)
      targetRoots <- MetricSqrt.factorDense(targetDense, eigenSolver, tolerance, "task target metric")
      _ <-
        if targetRoots.rank < targetDense.rows then
          Left(
            MultivarError.NonInvertibleValue(
              "task target metric eigenvalue",
              targetDense.rows - 1,
              0.0
            )
          )
        else Right(())
      coordinates <- MetricCoordinates.cross(crossDense, sourceRoots, targetRoots)
      zeroObjective = TaskMath.traceProduct(targetMarginal, targetDense)
    yield
      new CanonicalTask(
        coordinates,
        penalized,
        crossDense,
        zeroObjective,
        TaskAxis("predictor features", crossDense.rows),
        TaskAxis("response features", crossDense.cols)
      )

  /** Table / GMD task \(K = Q^{1/2} X R^{1/2}\) after centering. */
  def fromTable(
      working: DMat,
      rowMetric: MetricSpec,
      featureMetric: MetricSpec,
      eigenSolver: SymmetricEigenSolver,
      tolerance: Double,
      policy: StoragePolicy
  ): Either[MultivarError, CanonicalTask] =
    for
      _ <-
        if rowMetric.dim == working.rows then Right(())
        else Left(MultivarError.MetricShapeMismatch(IndexAxis.Row, working.rows, rowMetric.dim))
      _ <-
        if featureMetric.dim == working.cols then Right(())
        else Left(MultivarError.MetricShapeMismatch(IndexAxis.Feature, working.cols, featureMetric.dim))
      coordinates <- MetricCoordinates.table(
        working,
        rowMetric,
        featureMetric,
        eigenSolver,
        tolerance,
        policy,
        "task row metric",
        "task feature metric"
      )
      available = TaskMath.frobeniusNorm2(coordinates.matrix)
    yield
      new CanonicalTask(
        coordinates,
        DMat.eye(working.rows),
        working,
        available,
        TaskAxis("observations", working.rows),
        TaskAxis("features", working.cols)
      )

  private def taskSemantic[A](result: Either[SemanticError, A]): Either[MultivarError, A] =
    result.left.map:
      case SemanticError.MultivarFailure(error)  => error
      case SemanticError.LinearMapFailure(error) => LinalgErrorAdapter.toMultivarError(error)
      case error                                 => MultivarError.SolverFailed(error.message)

/** SPD metric roots matching the paired certification cutoff (full rank required). */
private[paired] object TaskSpdRoots:
  def factor(
      gram: DMat,
      eigenSolver: SymmetricEigenSolver,
      tolerance: Double,
      role: String
  ): Either[MultivarError, MetricRoots] =
    for
      _ <- MatrixOps.checkFinite(role, gram)
      _ <- MatrixOps.checkSymmetric(gram, tolerance)
      symmetrized = MatrixOps.symmetrize(gram)
      eigen <- LinalgErrorAdapter.adapt(eigenSolver.decompose(symmetrized))
      maximum = eigen.values(0)
      minimum = eigen.values(eigen.values.length - 1)
      scale = TaskMath.frobeniusNorm(gram)
      cutoff = math.max(tolerance, 1e-12) * math.max(1.0, math.max(scale, math.max(math.abs(maximum), math.abs(minimum))))
      _ <-
        if minimum <= cutoff then
          Left(MultivarError.NonInvertibleValue(s"$role eigenvalue", eigen.values.length - 1, minimum))
        else Right(())
      roots <- MetricSqrt.factorDense(symmetrized, eigenSolver, tolerance, role)
      _ <-
        if roots.rank < gram.rows then
          Left(MultivarError.NonInvertibleValue(s"$role rank", roots.rank, minimum))
        else Right(())
    yield roots

private[paired] object TaskMath:
  def frobeniusNorm(matrix: DMat): Double =
    math.sqrt(frobeniusNorm2(matrix))

  def frobeniusNorm2(matrix: DMat): Double =
    var sum = 0.0
    var row = 0
    while row < matrix.rows do
      var col = 0
      while col < matrix.cols do
        val value = matrix(row, col)
        sum += value * value
        col += 1
      row += 1
    sum

  def traceProduct(left: DMat, right: DMat): Double =
    var sum = 0.0
    var row = 0
    while row < left.rows do
      var col = 0
      while col < left.cols do
        sum += left(row, col) * right(row, col)
        col += 1
      row += 1
    sum

  def reconstruction(svd: SvdResult): DMat =
    if svd.singularValues.length == 0 then
      DMat.zeros(svd.u.rows, svd.v.rows)
    else
      val scaled = MetricOperator.scaleColumnsDense(svd.u, svd.singularValues)
      GaleNumerics.multiply(scaled, svd.v.transpose)

  def add(left: DMat, right: DMat): DMat =
    require(left.rows == right.rows && left.cols == right.cols, "matrix add shape mismatch")
    val out = Matrix.newBuilder(left.rows, left.cols)
    var row = 0
    while row < left.rows do
      var col = 0
      while col < left.cols do
        out(row, col) = left(row, col) + right(row, col)
        col += 1
      row += 1
    out.result()

  def takeSvd(svd: SvdResult, count: Int): SvdResult =
    if count <= 0 then
      SvdResult(
        DMat.zeros(svd.u.rows, 0),
        GaleNumerics.vectorFromArray(Array.empty[Double]),
        DMat.zeros(svd.v.rows, 0)
      )
    else if count >= svd.singularValues.length then svd
    else
      SvdResult(
        MatrixOps.takeColumns(svd.u, count),
        MatrixOps.takeVector(svd.singularValues, count),
        MatrixOps.takeColumns(svd.v, count)
      )

  /** Supported rank and count of weak-but-supported singular values of \(K\). */
  def supportSummary(values: DVec, tolerance: Double): (Int, Int) =
    var largest = 0.0
    var i = 0
    while i < values.length do
      val absolute = math.abs(values(i))
      if absolute > largest then largest = absolute
      i += 1
    val supportCut = tolerance * math.max(1.0, largest)
    val weakCut = math.sqrt(tolerance) * math.max(1.0, largest)
    var support = 0
    var weak = 0
    i = 0
    while i < values.length do
      val absolute = math.abs(values(i))
      if absolute > supportCut then
        support += 1
        if absolute <= weakCut then weak += 1
      i += 1
    (support, weak)

private[paired] object LocalTaskSvd:
  def fit(
      block: DMat,
      maxComponents: Int,
      solver: SvdSolver,
      tolerance: Double
  ): Either[MultivarError, SvdResult] =
    val limit = math.min(maxComponents, math.min(block.rows, block.cols))
    if limit <= 0 || TaskMath.frobeniusNorm2(block) <= 0.0 then
      Right(
        SvdResult(
          DMat.zeros(block.rows, 0),
          GaleNumerics.vectorFromArray(Array.empty[Double]),
          DMat.zeros(block.cols, 0)
        )
      )
    else
      for
        request <- ComponentCount(limit)
        svd <- solver.decompose(MatrixView.dense(block), request)
      yield svd

private[paired] object SpectralTransport:
  def dualApply(
      task: DMat,
      target: DMat,
      policy: DualSpectrum,
      solver: SvdSolver
  ): Either[MultivarError, DMat] =
    for
      request <- requestedRank(task, policy)
      svd <- solver.decompose(MatrixView.dense(task), request)
      filter <- spectralFilter(svd.singularValues, policy)
    yield
      val rightCoordinates = GaleNumerics.multiply(svd.v.transpose, target)
      val filtered = MatrixView.scaleRows(rightCoordinates, filter)
      GaleNumerics.multiply(svd.u, filtered)

  private def requestedRank(
      task: DMat,
      policy: DualSpectrum
  ): Either[MultivarError, ComponentCount] =
    val full = math.min(task.rows, task.cols)
    policy match
      case DualSpectrum.MoorePenrose(_) =>
        ComponentCount(math.max(1, full))
      case DualSpectrum.Truncated(components, _) =>
        ComponentCount(components).flatMap: count =>
          if count.value > full then Left(MultivarError.InvalidComponentRequest(count.value, full))
          else Right(count)
      case DualSpectrum.Ridge(_) =>
        ComponentCount(math.max(1, full))

  private def spectralFilter(
      values: DVec,
      policy: DualSpectrum
  ): Either[MultivarError, DVec] =
    policy match
      case DualSpectrum.MoorePenrose(tolerance) =>
        RowGeometryOps.requireTolerance("dual moore-penrose tolerance", tolerance).map: _ =>
          mapValues(values): (d, _) =>
            if supported(d, values, tolerance) then 1.0 / d else 0.0
      case DualSpectrum.Truncated(components, tolerance) =>
        for
          count <- ComponentCount(components)
          _ <- RowGeometryOps.requireTolerance("dual truncated tolerance", tolerance)
        yield
          mapValues(values): (d, index) =>
            if index < count.value && supported(d, values, tolerance) then 1.0 / d else 0.0
      case DualSpectrum.Ridge(lambda) =>
        if !lambda.isFinite || lambda <= 0.0 then
          Left(MultivarError.InvalidRegularization("task dual ridge", lambda, "lambda > 0"))
        else
          Right(
            mapValues(values): (d, _) =>
              d / (d * d + lambda)
          )

  private def supported(value: Double, values: DVec, tolerance: Double): Boolean =
    var largest = 0.0
    var i = 0
    while i < values.length do
      val absolute = math.abs(values(i))
      if absolute > largest then largest = absolute
      i += 1
    math.abs(value) > tolerance * math.max(1.0, largest)

  private def mapValues(values: DVec)(f: (Double, Int) => Double): DVec =
    val out = new Array[Double](values.length)
    var i = 0
    while i < values.length do
      out(i) = f(values(i), i)
      i += 1
    GaleNumerics.vectorFromArray(out)

private[paired] object TaskResolver:
  def split(
      task: CanonicalTask,
      resolution: TaskResolution,
      svdSolver: SvdSolver,
      tolerance: Double
  ): Either[MultivarError, OrthonormalSubspace.Split] =
    resolution match
      case TaskResolution.Whole =>
        Right(
          OrthonormalSubspace.Split(
            focus = task.matrix,
            remainder = DMat.zeros(task.matrix.rows, task.matrix.cols)
          )
        )
      case value: TaskResolution.TargetWeights =>
        for
          _ <- requireTargetRows(value.design, task.targetAxis)
          standardized = task.coordinates.right.pinvHalf.applyLeft(value.design)
          generators = GaleNumerics.multiply(task.matrix, standardized)
          split <- splitAround(task.matrix, generators, svdSolver, tolerance)
        yield split
      case value: TaskResolution.TargetStructure =>
        for
          _ <- requireTargetRows(value.design, task.targetAxis)
          standardized = task.coordinates.right.half.applyLeft(value.design)
          generators <- SpectralTransport.dualApply(
            task.matrix,
            standardized,
            value.dual,
            svdSolver
          )
          split <- splitAround(task.matrix, generators, svdSolver, tolerance)
        yield split
      case value: TaskResolution.SourceCovariance =>
        for
          _ <- requireSourceRows(value.design, task.sourceAxis, "covariance")
          // G* = G_ρ^{1/2} G
          standardized = task.coordinates.left.half.applyLeft(value.design)
          // K (K' G*) without forming K K'
          cross = GaleNumerics.multiply(task.matrix.transpose, standardized)
          generators = GaleNumerics.multiply(task.matrix, cross)
          split <- splitAround(task.matrix, generators, svdSolver, tolerance)
        yield split
      case value: TaskResolution.SourceRegression =>
        for
          _ <- requireSourceRows(value.design, task.sourceAxis, "regression")
          standardized = task.coordinates.left.half.applyLeft(value.design)
          support <- OrthonormalSubspace.span(task.matrix, svdSolver, tolerance)
          generators <- support.project(standardized)
          split <- splitAround(task.matrix, generators, svdSolver, tolerance)
        yield split

  private def splitAround(
      task: DMat,
      generators: DMat,
      solver: SvdSolver,
      tolerance: Double
  ): Either[MultivarError, OrthonormalSubspace.Split] =
    OrthonormalSubspace.span(generators, solver, tolerance).flatMap(_.split(task))

  private def requireTargetRows(design: DMat, axis: TaskAxis): Either[MultivarError, Unit] =
    if design.rows == axis.dim then Right(())
    else
      Left(
        MultivarError.MatrixShapeMismatch(
          s"Target resolution design has ${design.rows} rows, but this task's target axis is the ${axis.dim} ${axis.label}"
        )
      )

  private def requireSourceRows(
      design: DMat,
      axis: TaskAxis,
      kind: String
  ): Either[MultivarError, Unit] =
    if design.rows == axis.dim then Right(())
    else
      Left(
        MultivarError.MatrixShapeMismatch(
          s"Source.$kind design has ${design.rows} rows, but this task's source axis is the ${axis.dim} ${axis.label}"
        )
      )

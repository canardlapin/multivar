package multivar
package family.paired

import multivar.capability.*
import multivar.core.*
import multivar.family.spectral.GpcaCentering

import gale.linalg.DMat
import gale.linalg.DVec

/** Opaque predictive result of a resolved least-squares task.
  *
  * Contributions cancel the response affine baseline so
  * `predict(x) = baseline + focusContribution(x) + remainderContribution(x)`.
  */
final class TaskComponentsFit private[paired] (
    val focus: TaskComponentBlock,
    val remainder: TaskComponentBlock,
    val trainingValue: TaskValue,
    val sourceAxis: String,
    val targetAxis: String,
    val resolution: TaskResolution,
    val rank: RankBudget,
    val dualSpectrum: Option[DualSpectrum],
    val boundaryTied: Boolean,
    /** Numerical rank of the canonical task matrix \(K\). */
    val taskSupportRank: Int,
    /** Singular values of \(K\) below \(\sqrt{\tau}\sigma_{\max}\) but still supported. */
    val weakDirectionCount: Int,
    private val rawCoefficients: DMat,
    private val rawIntercept: DVec,
    private val focusWorkingMap: DMat,
    private val remainderWorkingMap: DMat,
    private val fullTransform: FittedCoefficientTransform,
    private val predictorPreprocessor: FittedPreprocessor,
    private val responsePreprocessor: FittedInvertiblePreprocessor
):
  def value: TaskValue =
    trainingValue

  def coefficients: DMat =
    rawCoefficients

  def intercept: DVec =
    rawIntercept

  def predict(input: DMat): Either[MultivarError, DMat] =
    fullTransform.predict(MatrixView.dense(input))

  def predict(input: MatrixView): Either[MultivarError, DMat] =
    fullTransform.predict(input)

  /** Response baseline (inverse of a zero working prediction). */
  def baseline(rows: Int): Either[MultivarError, DMat] =
    responsePreprocessor.inverseTransformDense(DMat.zeros(rows, responsePreprocessor.inputCols))

  def focusContribution(input: DMat): Either[MultivarError, DMat] =
    contribution(input, focusWorkingMap)

  def remainderContribution(input: DMat): Either[MultivarError, DMat] =
    contribution(input, remainderWorkingMap)

  private def contribution(
      input: DMat,
      coefficient: DMat
  ): Either[MultivarError, DMat] =
    for
      prepared <- predictorPreprocessor.transform(
        MatrixView.dense(input),
        policy = StoragePolicy.AllowDense
      )
      dense <- prepared.toDense(StoragePolicy.AllowDense)
      working = GaleNumerics.multiply(dense, coefficient)
      original <- responsePreprocessor.inverseContributionDense(working)
    yield original

/** Low-rank component analysis of a least-squares task, optionally resolved.
  *
  * Use [[fit]] for prediction \(X\to Y\) and [[decompose]] for a metric table
  * with a nontrivial scientific resolution.
  */
object TaskComponents:
  final case class Options(
      regularization: RegressionRegularization = RegressionRegularization.Ols,
      predictorPreprocessing: PreprocessSpec = PreprocessSpec.Center,
      responsePreprocessing: PreprocessSpec = PreprocessSpec.Center,
      rowMetric: Option[MetricSpec] = None,
      targetMetric: Option[MetricSpec] = None,
      rankTolerance: Double = 1e-12
  )

  def fit(
      predictors: DMat,
      responses: DMat,
      components: Int
  ): Either[MultivarError, TaskComponentsFit] =
    ComponentCount(components).flatMap: _ =>
      fit(
        MatrixView.dense(predictors),
        MatrixView.dense(responses),
        RankBudget.Focus(components),
        TaskResolution.Whole,
        Options(),
        DenseSolvers.svd,
        DenseSolvers.symmetricEigen,
        StoragePolicy.AllowDense
      )

  def fit(
      predictors: DMat,
      responses: DMat,
      components: Int,
      resolution: TaskResolution
  ): Either[MultivarError, TaskComponentsFit] =
    ComponentCount(components).flatMap: _ =>
      fit(
        MatrixView.dense(predictors),
        MatrixView.dense(responses),
        RankBudget.Focus(components),
        resolution,
        Options(),
        DenseSolvers.svd,
        DenseSolvers.symmetricEigen,
        StoragePolicy.AllowDense
      )

  def fit(
      predictors: DMat,
      responses: DMat,
      rank: RankBudget,
      resolution: TaskResolution
  ): Either[MultivarError, TaskComponentsFit] =
    fit(predictors, responses, rank, resolution, Options())

  def fit(
      predictors: DMat,
      responses: DMat,
      rank: RankBudget,
      resolution: TaskResolution,
      options: Options
  ): Either[MultivarError, TaskComponentsFit] =
    fit(
      MatrixView.dense(predictors),
      MatrixView.dense(responses),
      rank,
      resolution,
      options,
      DenseSolvers.svd,
      DenseSolvers.symmetricEigen,
      StoragePolicy.AllowDense
    )

  def fit(
      predictors: MatrixView,
      responses: MatrixView,
      rank: RankBudget,
      resolution: TaskResolution,
      options: Options,
      solver: SvdSolver,
      eigenSolver: SymmetricEigenSolver,
      policy: StoragePolicy
  ): Either[MultivarError, TaskComponentsFit] =
    for
      _ <- validateOptions(options)
      budget <- RankBudget.validate(rank)
      prepared <- PreparedPair.fromPredictive(
        predictors,
        responses,
        options.predictorPreprocessing,
        options.responsePreprocessing,
        options.rowMetric,
        policy
      )
      fittedY <- prepared.invertibleResponse
      targetMetric <- options.targetMetric match
        case Some(metric) => Right(metric)
        case None         => MetricSpec.identity(prepared.moments.yFeatures)
      _ <- validateComponentBudget(budget, prepared.moments)
      task <- CanonicalTask.fromPrepared(
        prepared,
        options.regularization,
        targetMetric,
        eigenSolver,
        options.rankTolerance,
        policy
      )
      (xWorking, _) <- prepared.workingDense(policy)
      core <- ResolvedTaskProblem(task, resolution).fit(budget, solver, options.rankTolerance)
      focusWorking = task.workingMap(core.focusSvd)
      remainderWorking = task.workingMap(core.remainderSvd)
      working = TaskMath.add(focusWorking, remainderWorking)
      rawMap <- PairedCoordinateMap.decode(working, prepared.xPreprocessor, fittedY)
      transform <- FittedCoefficientTransform.from(
        working,
        prepared.xPreprocessor,
        fittedY,
        "task-components"
      )
      focusBlock =
        if core.focusSvd.singularValues.length == 0 then
          TaskComponentBlock(
            TaskPart.Focus,
            GaleNumerics.vectorFromArray(Array.empty[Double]),
            DMat.zeros(task.sourceAxis.dim, 0),
            DMat.zeros(task.targetAxis.dim, 0),
            DMat.zeros(xWorking.rows, 0),
            TaskPartValue(core.focusAvailable, core.focusRetained, task.availableGain)
          )
        else
          blockFrom(
            TaskPart.Focus,
            core.focusSvd,
            task,
            xWorking,
            core.focusAvailable,
            core.focusRetained,
            task.availableGain
          )
      remainderBlock =
        if core.remainderSvd.singularValues.length == 0 then
          TaskComponentBlock(
            TaskPart.Remainder,
            GaleNumerics.vectorFromArray(Array.empty[Double]),
            DMat.zeros(task.sourceAxis.dim, 0),
            DMat.zeros(task.targetAxis.dim, 0),
            DMat.zeros(xWorking.rows, 0),
            TaskPartValue(core.remainderAvailable, core.remainderRetained, task.availableGain)
          )
        else
          blockFrom(
            TaskPart.Remainder,
            core.remainderSvd,
            task,
            xWorking,
            core.remainderAvailable,
            core.remainderRetained,
            task.availableGain
          )
      value = TaskValue(
        baselineLoss = task.zeroObjective,
        fullLoss = task.fullLoss,
        fittedLoss = math.max(0.0, task.zeroObjective - (core.focusRetained + core.remainderRetained)),
        availableGain = task.availableGain,
        retainedGain = core.focusRetained + core.remainderRetained
      )
      dual = resolution match
        case TaskResolution.TargetStructure(_, dualSpectrum) => Some(dualSpectrum)
        case _                                               => None
    yield
      new TaskComponentsFit(
        focusBlock,
        remainderBlock,
        value,
        task.sourceAxis.label,
        task.targetAxis.label,
        resolution,
        budget,
        dual,
        core.boundaryTied,
        core.taskSupportRank,
        core.weakDirectionCount,
        rawMap.coefficients,
        rawMap.intercept,
        focusWorking,
        remainderWorking,
        transform,
        prepared.xPreprocessor,
        fittedY
      )

  private def validateOptions(options: Options): Either[MultivarError, Unit] =
    val ridgeValue =
      options.regularization match
        case RegressionRegularization.Ols          => 0.0
        case RegressionRegularization.Ridge(value) => value.value
    val identityTarget = options.targetMetric.forall(_.isIdentity)
    if ridgeValue > 0.0 && !identityTarget then
      Left(
        MultivarError.InvalidRegularization(
          "task-components target metric with Frobenius ridge",
          ridgeValue,
          "Phase 1 rejects nonidentity target metric R with Frobenius ridge (Sylvester gap); use identity R with ridge, or arbitrary R with OLS"
        )
      )
    else RowGeometryOps.requireTolerance("task rank tolerance", options.rankTolerance)

  private def validateComponentBudget(
      budget: RankBudget,
      moments: PairedMoments
  ): Either[MultivarError, Unit] =
    val limit = math.min(moments.sampleCount, math.min(moments.xFeatures, moments.yFeatures))
    val requested =
      budget match
        case RankBudget.Focus(components)       => components
        case RankBudget.Total(components)       => components
        case RankBudget.Split(focus, remainder) => focus + remainder
    if requested > limit then Left(MultivarError.InvalidComponentRequest(requested, limit))
    else Right(())

  final case class DecompositionOptions(
      rowMetric: Option[MetricSpec] = None,
      featureMetric: Option[MetricSpec] = None,
      centering: GpcaCentering = GpcaCentering.Auto,
      rankTolerance: Double = 1e-12
  )

  /** Resolved metric decomposition. Requires a nontrivial scientific resolution.
    *
    * Unresolved metric SVD remains [[multivar.family.spectral.Gpca]]; do not use
    * `TaskResolution.Whole` on this ordinary entry point.
    */
  def decompose(
      input: DMat,
      resolution: TaskResolution,
      components: Int
  ): Either[MultivarError, TaskDecompositionFit] =
    ComponentCount(components).flatMap: _ =>
      decompose(
        MatrixView.dense(input),
        RankBudget.Focus(components),
        resolution,
        DecompositionOptions(),
        DenseSolvers.svd,
        DenseSolvers.symmetricEigen,
        StoragePolicy.AllowDense
      )

  def decompose(
      input: DMat,
      rank: RankBudget,
      resolution: TaskResolution
  ): Either[MultivarError, TaskDecompositionFit] =
    decompose(input, rank, resolution, DecompositionOptions())

  def decompose(
      input: DMat,
      rank: RankBudget,
      resolution: TaskResolution,
      options: DecompositionOptions
  ): Either[MultivarError, TaskDecompositionFit] =
    decompose(
      MatrixView.dense(input),
      rank,
      resolution,
      options,
      DenseSolvers.svd,
      DenseSolvers.symmetricEigen,
      StoragePolicy.AllowDense
    )

  def decompose(
      input: MatrixView,
      rank: RankBudget,
      resolution: TaskResolution,
      options: DecompositionOptions,
      solver: SvdSolver,
      eigenSolver: SymmetricEigenSolver,
      policy: StoragePolicy
  ): Either[MultivarError, TaskDecompositionFit] =
    TaskDecompositionEngine.run(
      input,
      rank,
      resolution,
      options,
      solver,
      eigenSolver,
      policy,
      allowWhole = false
    )

  private def blockFrom(
      part: TaskPart,
      svd: SvdResult,
      task: CanonicalTask,
      workingX: DMat,
      available: Double,
      retained: Double,
      totalAvailable: Double
  ): TaskComponentBlock =
    TaskComponentBlock(
      part,
      svd.singularValues,
      task.coordinates.sourceWeights(svd.u),
      task.coordinates.targetWeights(svd.v),
      task.sourceScores(workingX, svd),
      TaskPartValue(available, retained, totalAvailable)
    )

private[paired] final case class ResolvedTaskCore(
    focusSvd: SvdResult,
    remainderSvd: SvdResult,
    focusBlock: DMat,
    remainderBlock: DMat,
    focusAvailable: Double,
    remainderAvailable: Double,
    focusRetained: Double,
    remainderRetained: Double,
    boundaryTied: Boolean,
    taskSupportRank: Int,
    weakDirectionCount: Int,
    componentOrder: Vector[TaskComponentIndex]
)

private[paired] final class ResolvedTaskProblem(
    task: CanonicalTask,
    resolution: TaskResolution
):
  def fit(
      budget: RankBudget,
      solver: SvdSolver,
      tolerance: Double
  ): Either[MultivarError, ResolvedTaskCore] =
    for
      split <- TaskResolver.split(task, resolution, solver, tolerance)
      maxFocus =
        budget match
          case RankBudget.Focus(k)        => k
          case RankBudget.Total(k)        => k
          case RankBudget.Split(focus, _) => focus
      maxRemainder =
        budget match
          case RankBudget.Focus(_)            => 0
          case RankBudget.Total(k)            => k
          case RankBudget.Split(_, remainder) => remainder
      focusCandidate <- LocalTaskSvd.fit(split.focus, maxFocus, solver, tolerance)
      remainderCandidate <- LocalTaskSvd.fit(split.remainder, maxRemainder, solver, tolerance)
      allocation = RankAllocator.allocate(
        focusCandidate.singularValues,
        remainderCandidate.singularValues,
        budget
      )
      focusSvd = TaskMath.takeSvd(focusCandidate, allocation.focus)
      remainderSvd = TaskMath.takeSvd(remainderCandidate, allocation.remainder)
      spectrum <- LocalTaskSvd.fit(
        task.matrix,
        math.min(task.matrix.rows, task.matrix.cols),
        solver,
        tolerance
      )
      (supportRank, weakCount) = TaskMath.supportSummary(spectrum.singularValues, tolerance)
    yield
      ResolvedTaskCore(
        focusSvd,
        remainderSvd,
        focusBlock = split.focus,
        remainderBlock = split.remainder,
        focusAvailable = TaskMath.frobeniusNorm2(split.focus),
        remainderAvailable = TaskMath.frobeniusNorm2(split.remainder),
        focusRetained = TaskMath.frobeniusNorm2(TaskMath.reconstruction(focusSvd)),
        remainderRetained = TaskMath.frobeniusNorm2(TaskMath.reconstruction(remainderSvd)),
        boundaryTied = allocation.boundaryTied,
        taskSupportRank = supportRank,
        weakDirectionCount = weakCount,
        componentOrder = allocation.order
      )

private[paired] object RankAllocator:
  def allocate(
      focus: DVec,
      remainder: DVec,
      budget: RankBudget
  ): RankAllocation =
    budget match
      case RankBudget.Focus(components) =>
        val count = math.min(components, focus.length)
        new RankAllocation(
          count,
          0,
          boundaryTied = false,
          order = componentIndex(TaskPart.Focus, focus, count, globalOffset = 0)
        )
      case RankBudget.Split(focusCount, remainderCount) =>
        val kf = math.min(focusCount, focus.length)
        val kr = math.min(remainderCount, remainder.length)
        new RankAllocation(
          kf,
          kr,
          boundaryTied = false,
          order =
            componentIndex(TaskPart.Focus, focus, kf, globalOffset = 0) ++
              componentIndex(TaskPart.Remainder, remainder, kr, globalOffset = kf)
        )
      case RankBudget.Total(components) =>
        allocateTotal(focus, remainder, components)

  private def allocateTotal(
      focus: DVec,
      remainder: DVec,
      total: Int
  ): RankAllocation =
    val candidateBuilder = Vector.newBuilder[TaskComponentIndex]
    var index = 0
    while index < focus.length do
      candidateBuilder += TaskComponentIndex(0, TaskPart.Focus, index, focus(index))
      index += 1
    index = 0
    while index < remainder.length do
      candidateBuilder += TaskComponentIndex(0, TaskPart.Remainder, index, remainder(index))
      index += 1
    val ordered = candidateBuilder.result().sorted(using CandidateOrdering)
    val selected = ordered.take(total)
    val boundaryTied =
      if selected.isEmpty || ordered.length <= total then false
      else
        val cut = selected.last.singularValue
        val next = ordered(total).singularValue
        math.abs(cut - next) <= 1e-14 * math.max(1.0, math.abs(cut))
    var focusCount = 0
    var remainderCount = 0
    val orderBuilder = Vector.newBuilder[TaskComponentIndex]
    index = 0
    while index < selected.length do
      val candidate = selected(index)
      candidate.part match
        case TaskPart.Focus     => focusCount += 1
        case TaskPart.Remainder => remainderCount += 1
      orderBuilder += TaskComponentIndex(
        index,
        candidate.part,
        candidate.localIndex,
        candidate.singularValue
      )
      index += 1
    new RankAllocation(
      focus = focusCount,
      remainder = remainderCount,
      boundaryTied = boundaryTied,
      order = orderBuilder.result()
    )

  private def componentIndex(
      part: TaskPart,
      values: DVec,
      count: Int,
      globalOffset: Int
  ): Vector[TaskComponentIndex] =
    val builder = Vector.newBuilder[TaskComponentIndex]
    var index = 0
    while index < count do
      builder += TaskComponentIndex(
        globalOffset + index,
        part,
        index,
        values(index)
      )
      index += 1
    builder.result()

  private object CandidateOrdering extends Ordering[TaskComponentIndex]:
    def compare(left: TaskComponentIndex, right: TaskComponentIndex): Int =
      val byValue = java.lang.Double.compare(right.singularValue, left.singularValue)
      if byValue != 0 then byValue
      else
        val byPart = Integer.compare(left.part.ordinal, right.part.ordinal)
        if byPart != 0 then byPart
        else Integer.compare(left.localIndex, right.localIndex)

package multivar
package family.paired

import multivar.capability.FittedBidirectionalTransform
import multivar.capability.FittedFrameTransform
import multivar.core.*
import multivar.family.spectral.Gpca

import gale.linalg.DMat
import gale.linalg.DVec

/** Declared held-out evaluation loss for a fitted decomposition.
  *
  * Distinct from [[TaskValue]]: training gains are estimand identities, not
  * out-of-sample utility.
  */
enum DecompositionValidationLoss:
  /** Squared Frobenius residual in original feature coordinates. */
  case ReconstructionFrobenius

/** Held-out reconstruction evaluation (not a training [[TaskValue]]). */
final case class DecompositionValidation(
    baselineLoss: Double,
    fittedLoss: Double,
    residualFrobenius: Double
)

/** Diagnostic split of new observations under the fitted resolution transfer.
  *
  * Parts are original-feature contributions (baseline cancelled). They are not a
  * third scientific branch: `unsupported` is the complement of the training task
  * support \(P_{Z^\top}\), while focus/remainder are the full scientific split
  * \(F_f,F_r\) before rank truncation.
  *
  * In whitened features: \(z = z F_f + z F_r + z F_0\).
  */
final case class ResolvedInput(
    focus: DMat,
    supportedRemainder: DMat,
    unsupported: DMat
)

/** One resolved branch with inductive analysis / contribution synthesis. */
final class TaskDecompositionBranch private[paired] (
    val part: TaskPart,
    private val singular: DVec,
    private val source: DMat,
    private val target: DMat,
    private val trainingSourceScores: DMat,
    val value: TaskPartValue,
    private val synthesis: Option[FittedBidirectionalTransform],
    private val featureDim: Int,
    private val preprocessor: FittedInvertiblePreprocessor
):
  def singularValues: DVec =
    singular

  def sourceWeights: DMat =
    source

  def targetPatterns: DMat =
    target

  def sourceScores: DMat =
    trainingSourceScores

  def rank: Int =
    singular.length

  def componentValues: DVec =
    val out = new Array[Double](singular.length)
    var i = 0
    while i < singular.length do
      out(i) = singular(i) * singular(i)
      i += 1
    GaleNumerics.vectorFromArray(out)

  /** Express new observations in this branch's retained component coordinates. */
  def transform(input: DMat): Either[MultivarError, DMat] =
    transform(MatrixView.dense(input))

  def transform(input: MatrixView): Either[MultivarError, DMat] =
    synthesis match
      case Some(map) => map.analysis.project(input)
      case None =>
        Right(DMat.zeros(input.rows, 0))

  /** Branch reconstruction contribution (no shared baseline). */
  def reconstructContribution(input: DMat): Either[MultivarError, DMat] =
    reconstructContribution(MatrixView.dense(input))

  def reconstructContribution(input: MatrixView): Either[MultivarError, DMat] =
    synthesis match
      case Some(map) => map.reconstructContribution(input)
      case None =>
        Right(DMat.zeros(input.rows, featureDim)).flatMap(preprocessor.inverseContributionDense)

/** Opaque resolved decomposition of a metric table task with inductive transfer.
  *
  * Feature-side maps \(F_b = Z^{+}Z_b\) are folded into original-feature
  * analysis weights \(R^{1/2}A_b\) and synthesis \(V_b^\top R^{-1/2}\) so the
  * reusable bidirectional capability owns project / contribution paths.
  */
final class TaskDecompositionFit private[paired] (
    val focus: TaskDecompositionBranch,
    val remainder: TaskDecompositionBranch,
    val trainingValue: TaskValue,
    val sourceAxis: String,
    val targetAxis: String,
    val resolution: TaskResolution,
    val resolutionTransfer: TaskResolutionTransfer,
    val rank: RankBudget,
    val dualSpectrum: Option[DualSpectrum],
    val boundaryTied: Boolean,
    val taskSupportRank: Int,
    val weakDirectionCount: Int,
    val componentIndex: Vector[TaskComponentIndex],
    private val combined: Option[FittedBidirectionalTransform],
    private val focusWorking: DMat,
    private val remainderWorking: DMat,
    private val preprocessor: FittedInvertiblePreprocessor,
    private val featureDim: Int,
    private val featureRoots: MetricRoots,
    private val focusTransferFull: DMat,
    private val remainderTransferFull: DMat,
    private val unsupportedTransfer: DMat
):
  def value: TaskValue =
    trainingValue

  def reconstructWorking: DMat =
    TaskMath.add(focusWorking, remainderWorking)

  def transform(input: DMat): Either[MultivarError, DMat] =
    transform(MatrixView.dense(input))

  def transform(input: MatrixView): Either[MultivarError, DMat] =
    combined match
      case Some(map) => map.analysis.project(input)
      case None      => Right(DMat.zeros(input.rows, 0))

  def reconstruct(input: DMat): Either[MultivarError, DMat] =
    reconstruct(MatrixView.dense(input))

  def reconstruct(input: MatrixView): Either[MultivarError, DMat] =
    for
      base <- baseline(input.rows)
      focusPart <- focusContribution(input)
      remPart <- remainderContribution(input)
    yield TaskMath.add(base, TaskMath.add(focusPart, remPart))

  def reconstructTraining: Either[MultivarError, DMat] =
    preprocessor.inverseTransformDense(reconstructWorking)

  def baseline(rows: Int): Either[MultivarError, DMat] =
    preprocessor.inverseTransformDense(DMat.zeros(rows, preprocessor.inputCols))

  def focusContribution(input: DMat): Either[MultivarError, DMat] =
    focus.reconstructContribution(input)

  def focusContribution(input: MatrixView): Either[MultivarError, DMat] =
    focus.reconstructContribution(input)

  def remainderContribution(input: DMat): Either[MultivarError, DMat] =
    remainder.reconstructContribution(input)

  def remainderContribution(input: MatrixView): Either[MultivarError, DMat] =
    remainder.reconstructContribution(input)

  def focusContributionTraining: Either[MultivarError, DMat] =
    preprocessor.inverseContributionDense(focusWorking)

  def remainderContributionTraining: Either[MultivarError, DMat] =
    preprocessor.inverseContributionDense(remainderWorking)

  def residual(input: DMat): Either[MultivarError, DMat] =
    residual(MatrixView.dense(input), StoragePolicy.AllowDense)

  def residual(
      input: MatrixView,
      policy: StoragePolicy
  ): Either[MultivarError, DMat] =
    for
      dense <- input.toDense(policy)
      hat <- reconstruct(input)
    yield MatrixOps.subtract(dense, hat)

  /** Split new rows into scientific focus, supported remainder, and task-unsupported mass. */
  def resolve(input: DMat): Either[MultivarError, ResolvedInput] =
    resolve(MatrixView.dense(input), StoragePolicy.AllowDense)

  def resolve(
      input: MatrixView,
      policy: StoragePolicy
  ): Either[MultivarError, ResolvedInput] =
    for
      z <- whitenedFeatures(input, policy)
      focusW = GaleNumerics.multiply(z, focusTransferFull)
      remW = GaleNumerics.multiply(z, remainderTransferFull)
      unsupW = GaleNumerics.multiply(z, unsupportedTransfer)
      focusPart <- decodeWhitenedContribution(focusW)
      remPart <- decodeWhitenedContribution(remW)
      unsupPart <- decodeWhitenedContribution(unsupW)
    yield ResolvedInput(focusPart, remPart, unsupPart)

  /** Fraction of whitened feature energy outside the training task support. */
  def unsupportedFraction(input: DMat): Either[MultivarError, Double] =
    unsupportedFraction(MatrixView.dense(input), StoragePolicy.AllowDense)

  def unsupportedFraction(
      input: MatrixView,
      policy: StoragePolicy
  ): Either[MultivarError, Double] =
    whitenedFeatures(input, policy).map: z =>
      val total = TaskMath.frobeniusNorm2(z)
      if total <= 0.0 then 0.0
      else TaskMath.frobeniusNorm2(GaleNumerics.multiply(z, unsupportedTransfer)) / total

  def validationLoss(
      input: DMat,
      loss: DecompositionValidationLoss = DecompositionValidationLoss.ReconstructionFrobenius
  ): Either[MultivarError, DecompositionValidation] =
    validationLoss(MatrixView.dense(input), loss, StoragePolicy.AllowDense)

  def validationLoss(
      input: MatrixView,
      loss: DecompositionValidationLoss,
      policy: StoragePolicy
  ): Either[MultivarError, DecompositionValidation] =
    loss match
      case DecompositionValidationLoss.ReconstructionFrobenius =>
        for
          dense <- input.toDense(policy)
          base <- baseline(input.rows)
          hat <- reconstruct(input)
          resid = MatrixOps.subtract(dense, hat)
          aboutBaseline = MatrixOps.subtract(dense, base)
        yield
          DecompositionValidation(
            baselineLoss = TaskMath.frobeniusNorm2(aboutBaseline),
            fittedLoss = TaskMath.frobeniusNorm2(resid),
            residualFrobenius = TaskMath.frobeniusNorm(resid)
          )

  private def whitenedFeatures(
      input: MatrixView,
      policy: StoragePolicy
  ): Either[MultivarError, DMat] =
    for
      prepared <- preprocessor.transform(input, policy = policy)
      dense <- prepared.toDense(policy)
      _ <-
        if dense.cols != featureRoots.half.dim then
          Left(
            MultivarError.MatrixShapeMismatch(
              s"task decomposition expected ${featureRoots.half.dim} features, got ${dense.cols}"
            )
          )
        else Right(())
    yield featureRoots.half.applyRight(dense)

  private def decodeWhitenedContribution(whitened: DMat): Either[MultivarError, DMat] =
    val working = featureRoots.pinvHalf.applyRight(whitened)
    preprocessor.inverseContributionDense(working)

/** Package-private decomposition engine. Kept off the ordinary surface list so
  * the Whole/GMD law path does not appear on [[TaskComponents]].
  */
private[paired] object TaskDecompositionEngine:
  def run(
      input: MatrixView,
      rank: RankBudget,
      resolution: TaskResolution,
      options: TaskComponents.DecompositionOptions,
      solver: SvdSolver,
      eigenSolver: SymmetricEigenSolver,
      policy: StoragePolicy,
      allowWhole: Boolean
  ): Either[MultivarError, TaskDecompositionFit] =
    for
      _ <-
        if !allowWhole && resolution == TaskResolution.Whole then
          Left(
            MultivarError.UnsupportedEstimator(
              "TaskComponents.decompose requires a scientific resolution; use Gpca.fit for unresolved metric decomposition"
            )
          )
        else Right(())
      _ <- RowGeometryOps.requireTolerance("task decomposition rank tolerance", options.rankTolerance)
      budget <- RankBudget.validate(rank)
      rowMetric <- options.rowMetric match
        case Some(metric) => Right(metric)
        case None         => MetricSpec.identity(input.rows)
      featureMetric <- options.featureMetric match
        case Some(metric) => Right(metric)
        case None         => MetricSpec.identity(input.cols)
      fittedPrep <- Gpca.fitCentering(options.centering, rowMetric, input)
      invertible <- fittedPrep.requireInvertible
      preparedView <- fittedPrep.transform(input, policy = policy)
      working <- preparedView.toDense(policy)
      _ <- MatrixOps.checkFinite("task decomposition table", working)
      _ <- validateBudget(budget, working.rows, working.cols)
      task <- CanonicalTask.fromTable(
        working,
        rowMetric,
        featureMetric,
        eigenSolver,
        options.rankTolerance,
        policy
      )
      core <- ResolvedTaskProblem(task, resolution).fit(budget, solver, options.rankTolerance)
      zFocus = task.whitenedReconstruction(core.focusSvd)
      zRemainder = task.whitenedReconstruction(core.remainderSvd)
      focusTransfer <- FeatureTransfer.applyPseudoInverse(
        task.matrix,
        zFocus,
        solver,
        options.rankTolerance
      )
      remainderTransfer <- FeatureTransfer.applyPseudoInverse(
        task.matrix,
        zRemainder,
        solver,
        options.rankTolerance
      )
      focusTransferFull <- FeatureTransfer.applyPseudoInverse(
        task.matrix,
        core.focusBlock,
        solver,
        options.rankTolerance
      )
      remainderTransferFull <- FeatureTransfer.applyPseudoInverse(
        task.matrix,
        core.remainderBlock,
        solver,
        options.rankTolerance
      )
      taskSupport <- FeatureTransfer.applyPseudoInverse(
        task.matrix,
        task.matrix,
        solver,
        options.rankTolerance
      )
      unsupportedTransfer = MatrixOps.subtract(DMat.eye(task.matrix.cols), taskSupport)
      focusAnalysis = FeatureTransfer.encoder(focusTransfer, core.focusSvd.v)
      remainderAnalysis = FeatureTransfer.encoder(remainderTransfer, core.remainderSvd.v)
      (combinedAnalysis, combinedSynthesis, componentIndex) = combineComponents(
        core.componentOrder,
        focusAnalysis,
        remainderAnalysis,
        core.focusSvd.v,
        core.remainderSvd.v
      )
      focusWorking = task.decodeBlock(zFocus)
      remainderWorking = task.decodeBlock(zRemainder)
      focusMap <- biorthogonalMap(
        input,
        focusAnalysis,
        core.focusSvd.v,
        task.coordinates.right,
        invertible,
        "task-decomposition-focus",
        core.focusSvd.singularValues
      )
      remainderMap <- biorthogonalMap(
        input,
        remainderAnalysis,
        core.remainderSvd.v,
        task.coordinates.right,
        invertible,
        "task-decomposition-remainder",
        core.remainderSvd.singularValues
      )
      combinedMap <- biorthogonalMap(
        input,
        combinedAnalysis,
        combinedSynthesis,
        task.coordinates.right,
        invertible,
        "task-decomposition",
        combinedSpectrum(componentIndex)
      )
      focusBranch = branch(
        TaskPart.Focus,
        core.focusSvd,
        task,
        FeatureTransfer.canonicalScores(core.focusSvd),
        core.focusAvailable,
        core.focusRetained,
        task.availableGain,
        focusMap,
        invertible
      )
      remainderBranch = branch(
        TaskPart.Remainder,
        core.remainderSvd,
        task,
        FeatureTransfer.canonicalScores(core.remainderSvd),
        core.remainderAvailable,
        core.remainderRetained,
        task.availableGain,
        remainderMap,
        invertible
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
      new TaskDecompositionFit(
        focusBranch,
        remainderBranch,
        value,
        task.sourceAxis.label,
        task.targetAxis.label,
        resolution,
        resolutionTransfer(resolution),
        budget,
        dual,
        core.boundaryTied,
        core.taskSupportRank,
        core.weakDirectionCount,
        componentIndex,
        combinedMap,
        focusWorking,
        remainderWorking,
        invertible,
        working.cols,
        task.coordinates.right,
        focusTransferFull,
        remainderTransferFull,
        unsupportedTransfer
      )

  private def resolutionTransfer(resolution: TaskResolution): TaskResolutionTransfer =
    resolution match
      case TaskResolution.Whole | _: TaskResolution.TargetWeights | _: TaskResolution.TargetStructure =>
        TaskResolutionTransfer.FeatureSpecified
      case _: TaskResolution.SourceCovariance | _: TaskResolution.SourceRegression =>
        TaskResolutionTransfer.TrainingSourceInduced

  def decomposeAllowingWhole(
      input: MatrixView,
      rank: RankBudget,
      options: TaskComponents.DecompositionOptions,
      solver: SvdSolver = DenseSolvers.svd,
      eigenSolver: SymmetricEigenSolver = DenseSolvers.symmetricEigen,
      policy: StoragePolicy = StoragePolicy.AllowDense
  ): Either[MultivarError, TaskDecompositionFit] =
    run(
      input,
      rank,
      TaskResolution.Whole,
      options,
      solver,
      eigenSolver,
      policy,
      allowWhole = true
    )

  /** Fold whitened \(A,V\) into original-feature bidirectional maps via \(R^{\pm 1/2}\). */
  private def biorthogonalMap(
      training: MatrixView,
      whitenedAnalysis: DMat,
      whitenedSynthesis: DMat,
      featureRoots: MetricRoots,
      preprocessor: FittedInvertiblePreprocessor,
      method: String,
      spectrum: DVec
  ): Either[MultivarError, Option[FittedBidirectionalTransform]] =
    if whitenedAnalysis.cols == 0 then Right(None)
    else
      val weights = featureRoots.half.applyLeft(whitenedAnalysis)
      val decoder = featureRoots.pinvHalf.applyLeft(whitenedSynthesis).transpose
      for
        count <- ComponentCount(weights.cols)
        frame <- FittedFrameTransform.fromTraining(
          training,
          weights,
          preprocessor,
          method,
          count,
          Some(spectrum)
        )
        map <- frame.withExplicitSynthesis(
          decoder,
          ValueIdentity.source(ValueId.unsafe(s"$method.synthesis"))
        )
      yield Some(map)

  private def combineComponents(
      order: Vector[TaskComponentIndex],
      focusAnalysis: DMat,
      remainderAnalysis: DMat,
      focusSynthesis: DMat,
      remainderSynthesis: DMat
  ): (DMat, DMat, Vector[TaskComponentIndex]) =
    val cols = order.length
    val rows = math.max(focusAnalysis.rows, remainderAnalysis.rows)
    val analysis =
      if cols == 0 then DMat.zeros(rows, 0)
      else takeColumns(order, focusAnalysis, remainderAnalysis, rows)
    val synthesis =
      if cols == 0 then DMat.zeros(rows, 0)
      else takeColumns(order, focusSynthesis, remainderSynthesis, rows)
    (analysis, synthesis, order)

  private def takeColumns(
      order: Vector[TaskComponentIndex],
      focus: DMat,
      remainder: DMat,
      rows: Int
  ): DMat =
    val builder = gale.linalg.Matrix.newBuilder(rows, order.length)
    var j = 0
    while j < order.length do
      val ref = order(j)
      val source = ref.part match
        case TaskPart.Focus     => focus
        case TaskPart.Remainder => remainder
      var i = 0
      while i < rows do
        builder(i, j) = source(i, ref.localIndex)
        i += 1
      j += 1
    builder.result()

  private def combinedSpectrum(index: Vector[TaskComponentIndex]): DVec =
    GaleNumerics.vectorFromArray(index.map(_.singularValue).toArray)

  private def validateBudget(
      budget: RankBudget,
      rows: Int,
      cols: Int
  ): Either[MultivarError, Unit] =
    val limit = math.min(rows, cols)
    val requested =
      budget match
        case RankBudget.Focus(components)       => components
        case RankBudget.Total(components)       => components
        case RankBudget.Split(focus, remainder) => focus + remainder
    if requested > limit then Left(MultivarError.InvalidComponentRequest(requested, limit))
    else Right(())

  private def branch(
      part: TaskPart,
      svd: SvdResult,
      task: CanonicalTask,
      scores: DMat,
      available: Double,
      retained: Double,
      totalAvailable: Double,
      map: Option[FittedBidirectionalTransform],
      preprocessor: FittedInvertiblePreprocessor
  ): TaskDecompositionBranch =
    if svd.singularValues.length == 0 then
      new TaskDecompositionBranch(
        part,
        GaleNumerics.vectorFromArray(Array.empty[Double]),
        DMat.zeros(task.sourceAxis.dim, 0),
        DMat.zeros(task.targetAxis.dim, 0),
        DMat.zeros(task.sourceAxis.dim, 0),
        TaskPartValue(available, retained, totalAvailable),
        None,
        task.targetAxis.dim,
        preprocessor
      )
    else
      new TaskDecompositionBranch(
        part,
        svd.singularValues,
        task.coordinates.sourceWeights(svd.u),
        task.coordinates.targetWeights(svd.v),
        scores,
        TaskPartValue(available, retained, totalAvailable),
        map,
        task.targetAxis.dim,
        preprocessor
      )

package multivar
package family.paired

import multivar.core.*
import multivar.capability.*

import gale.linalg.DMat
import gale.linalg.DVec

/** Opaque PLSC result. Typed operator and frames open through [[multivar.advanced]]. */
final class PlscFit private[multivar] (
    private val svd: SvdResult,
    private val pairedOperator: PairedOperatorFit[? <: SemanticSpace, ? <: SemanticSpace, ? <: SemanticSpace],
    private val source: FittedFrameTransform,
    private val target: FittedFrameTransform
):
  def xScores: DMat = source.scores
  def yScores: DMat = target.scores
  def xWeights: DMat = source.weights
  def yWeights: DMat = target.weights
  def covariances: DVec = svd.singularValues

  def transformX(input: DMat): Either[MultivarError, DMat] =
    transformX(MatrixView.dense(input))

  def transformX(input: MatrixView): Either[MultivarError, DMat] =
    source.project(input)

  def transformY(input: DMat): Either[MultivarError, DMat] =
    transformY(MatrixView.dense(input))

  def transformY(input: MatrixView): Either[MultivarError, DMat] =
    target.project(input)

object PlscFit:
  private[multivar] def sourceOf(fit: PlscFit): FittedFrameTransform = fit.source
  private[multivar] def targetOf(fit: PlscFit): FittedFrameTransform = fit.target
  private[multivar] def operatorOf(fit: PlscFit): PairedOperatorFit[?, ?, ?] = fit.pairedOperator
  private[multivar] def resultOf(fit: PlscFit): SvdResult = fit.svd

object Plsc:
  def fit(
      x: DMat,
      y: DMat,
      components: Int
  ): Either[MultivarError, PlscFit] =
    fit(x, y, components, PreprocessSpec.Center, PreprocessSpec.Center)

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      xPreproc: PreprocessSpec,
      yPreproc: PreprocessSpec
  ): Either[MultivarError, PlscFit] =
    ComponentCount(components).flatMap: checked =>
      fit(MatrixView.dense(x), MatrixView.dense(y), checked, xPreproc, yPreproc)

  def fit(
      x: MatrixView,
      y: MatrixView,
      components: ComponentCount,
      xPreproc: PreprocessSpec = PreprocessSpec.Center,
      yPreproc: PreprocessSpec = PreprocessSpec.Center,
      solver: SvdSolver = DenseSolvers.svd,
      rowMetric: Option[MetricSpec] = None,
      eigenSolver: SymmetricEigenSolver = DenseSolvers.symmetricEigen,
      policy: StoragePolicy = StoragePolicy.AllowDense
  ): Either[MultivarError, PlscFit] =
    for
      prepared <- PreparedPair.from(x, y, xPreproc, yPreproc, rowMetric, policy)
      problem <- PairedOperatorProblem.fromMatrices(
        prepared.xWorking,
        prepared.yWorking,
        prepared.rowMetric,
        "plsc",
        policy
      )
      operator <- problem.fitPlsc(components, prepared.moments.covarianceScale, solver, eigenSolver)
      xWeights <- operator.sourceWeights
      yWeights <- operator.targetWeights
      xTransform <- FittedFrameTransform.fromTraining(
        prepared.xOriginal,
        xWeights,
        prepared.xPreprocessor,
        "plsc.source",
        components,
        Some(operator.result.singularValues)
      )
      yTransform <- FittedFrameTransform.fromTraining(
        prepared.yOriginal,
        yWeights,
        prepared.yPreprocessor,
        "plsc.target",
        components,
        Some(operator.result.singularValues)
      )
    yield PlscFit(operator.result, operator, xTransform, yTransform)

/** Opaque CCA result. Typed operator and frames open through [[multivar.advanced]]. */
final class CcaFit private[multivar] (
    private val svd: SvdResult,
    private val pairedOperator: PairedOperatorFit[? <: SemanticSpace, ? <: SemanticSpace, ? <: SemanticSpace],
    private val source: FittedFrameTransform,
    private val target: FittedFrameTransform
):
  def xScores: DMat = source.scores
  def yScores: DMat = target.scores
  def xWeights: DMat = source.weights
  def yWeights: DMat = target.weights
  def correlations: DVec = svd.singularValues

  def transformX(input: DMat): Either[MultivarError, DMat] =
    transformX(MatrixView.dense(input))

  def transformX(input: MatrixView): Either[MultivarError, DMat] =
    source.project(input)

  def transformY(input: DMat): Either[MultivarError, DMat] =
    transformY(MatrixView.dense(input))

  def transformY(input: MatrixView): Either[MultivarError, DMat] =
    target.project(input)

object CcaFit:
  private[multivar] def sourceOf(fit: CcaFit): FittedFrameTransform = fit.source
  private[multivar] def targetOf(fit: CcaFit): FittedFrameTransform = fit.target
  private[multivar] def operatorOf(fit: CcaFit): PairedOperatorFit[?, ?, ?] = fit.pairedOperator
  private[multivar] def resultOf(fit: CcaFit): SvdResult = fit.svd

object Cca:
  def fit(
      x: DMat,
      y: DMat,
      components: Int
  ): Either[MultivarError, CcaFit] =
    fit(x, y, components, ridge = 1e-8)

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      ridge: Double
  ): Either[MultivarError, CcaFit] =
    ComponentCount(components).flatMap: checked =>
      fit(MatrixView.dense(x), MatrixView.dense(y), checked, ridge)

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      xRidge: Double,
      yRidge: Double
  ): Either[MultivarError, CcaFit] =
    for
      checked <- ComponentCount(components)
      regularization <- CcaRegularization.asymmetric(xRidge, yRidge)
      fit <- fitRegularized(MatrixView.dense(x), MatrixView.dense(y), checked, regularization)
    yield fit

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      ridge: Double,
      xPreproc: PreprocessSpec,
      yPreproc: PreprocessSpec
  ): Either[MultivarError, CcaFit] =
    ComponentCount(components).flatMap: checked =>
      fit(MatrixView.dense(x), MatrixView.dense(y), checked, ridge, xPreproc, yPreproc)

  def fit(
      x: MatrixView,
      y: MatrixView,
      components: ComponentCount,
      ridge: Double = 1e-8,
      xPreproc: PreprocessSpec = PreprocessSpec.Center,
      yPreproc: PreprocessSpec = PreprocessSpec.Center,
      eigenSolver: SymmetricEigenSolver = DenseSolvers.symmetricEigen,
      solver: SvdSolver = DenseSolvers.svd,
      rowMetric: Option[MetricSpec] = None,
      policy: StoragePolicy = StoragePolicy.AllowDense
  ): Either[MultivarError, CcaFit] =
    CcaRegularization.symmetric(ridge).flatMap: regularization =>
      fitRegularized(x, y, components, regularization, xPreproc, yPreproc, eigenSolver, solver, rowMetric, policy)

  def fitRegularized(
      x: DMat,
      y: DMat,
      components: Int,
      regularization: CcaRegularization
  ): Either[MultivarError, CcaFit] =
    ComponentCount(components).flatMap: checked =>
      fitRegularized(MatrixView.dense(x), MatrixView.dense(y), checked, regularization)

  def fitRegularized(
      x: MatrixView,
      y: MatrixView,
      components: ComponentCount,
      regularization: CcaRegularization,
      xPreproc: PreprocessSpec = PreprocessSpec.Center,
      yPreproc: PreprocessSpec = PreprocessSpec.Center,
      eigenSolver: SymmetricEigenSolver = DenseSolvers.symmetricEigen,
      solver: SvdSolver = DenseSolvers.svd,
      rowMetric: Option[MetricSpec] = None,
      policy: StoragePolicy = StoragePolicy.AllowDense
  ): Either[MultivarError, CcaFit] =
    for
      prepared <- PreparedPair.from(x, y, xPreproc, yPreproc, rowMetric, policy)
      problem <- PairedOperatorProblem.fromMatrices(
        prepared.xWorking,
        prepared.yWorking,
        prepared.rowMetric,
        "cca",
        policy
      )
      operator <- problem.fitCca(
        components,
        regularization,
        prepared.moments.covarianceScale,
        solver,
        eigenSolver
      )
      xWeights <- operator.sourceWeights
      yWeights <- operator.targetWeights
      xTransform <- FittedFrameTransform.fromTraining(
        prepared.xOriginal,
        xWeights,
        prepared.xPreprocessor,
        "cca.source",
        components,
        Some(operator.result.singularValues)
      )
      yTransform <- FittedFrameTransform.fromTraining(
        prepared.yOriginal,
        yWeights,
        prepared.yPreprocessor,
        "cca.target",
        components,
        Some(operator.result.singularValues)
      )
    yield CcaFit(operator.result, operator, xTransform, yTransform)

/** Opaque reduced-rank regression result.
  *
  * [[coefficients]] and [[intercept]] are in original predictor/response units.
  * [[workingCoefficients]] is the map in preprocessed coordinates that
  * [[predictWorking]] applies after predictor preprocessing.
  */
final class ReducedRankRegressionFit private[multivar] (
    private val rawCoefficients: DMat,
    private val rawIntercept: DVec,
    private val workingMap: DMat,
    private val unconstrainedWorking: DMat,
    private val coefficientTransform: FittedCoefficientTransform,
    private val source: FittedFrameTransform,
    private val target: FittedFrameTransform,
    private val pairedOperator: PairedOperatorFit[? <: SemanticSpace, ? <: SemanticSpace, ? <: SemanticSpace]
):
  /** Rank-constrained coefficient map from raw predictors to raw responses. */
  def coefficients: DMat =
    rawCoefficients

  /** Intercept in original response units: `predict(x) ≈ x * coefficients + broadcast(intercept)`. */
  def intercept: DVec =
    rawIntercept

  /** Coefficient map in preprocessed coordinates. */
  def workingCoefficients: DMat =
    workingMap

  def xScores: DMat = source.scores
  def yScores: DMat = target.scores
  def xWeights: DMat = source.weights
  def yLoadings: DMat = target.weights

  def predictWorking(input: DMat): Either[MultivarError, DMat] =
    predictWorking(MatrixView.dense(input))

  def predictWorking(input: MatrixView): Either[MultivarError, DMat] =
    coefficientTransform.predictWorking(input)

  def predict(input: DMat): Either[MultivarError, DMat] =
    predict(MatrixView.dense(input))

  def predict(input: MatrixView): Either[MultivarError, DMat] =
    coefficientTransform.predict(input)

  def transformX(input: DMat): Either[MultivarError, DMat] =
    transformX(MatrixView.dense(input))

  def transformX(input: MatrixView): Either[MultivarError, DMat] =
    source.project(input)

  def transformY(input: DMat): Either[MultivarError, DMat] =
    transformY(MatrixView.dense(input))

  def transformY(input: MatrixView): Either[MultivarError, DMat] =
    target.project(input)

object ReducedRankRegressionFit:
  private[multivar] def sourceOf(fit: ReducedRankRegressionFit): FittedFrameTransform = fit.source
  private[multivar] def targetOf(fit: ReducedRankRegressionFit): FittedFrameTransform = fit.target
  private[multivar] def operatorOf(fit: ReducedRankRegressionFit): PairedOperatorFit[?, ?, ?] =
    fit.pairedOperator
  private[multivar] def coefficientTransformOf(fit: ReducedRankRegressionFit): FittedCoefficientTransform =
    fit.coefficientTransform
  private[multivar] def unconstrainedWorkingCoefficients(fit: ReducedRankRegressionFit): DMat =
    fit.unconstrainedWorking

object ReducedRankRegression:
  def fit(
      x: DMat,
      y: DMat,
      components: Int
  ): Either[MultivarError, ReducedRankRegressionFit] =
    fit(x, y, components, RegressionRegularization.Ols)

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      regularization: RegressionRegularization
  ): Either[MultivarError, ReducedRankRegressionFit] =
    fit(x, y, components, regularization, PreprocessSpec.Center, PreprocessSpec.Center)

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      ridge: Double
  ): Either[MultivarError, ReducedRankRegressionFit] =
    for
      regularization <- RegressionRegularization.ridge(ridge)
      fit <- fit(x, y, components, regularization)
    yield fit

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      regularization: RegressionRegularization,
      xPreproc: PreprocessSpec,
      yPreproc: PreprocessSpec
  ): Either[MultivarError, ReducedRankRegressionFit] =
    ComponentCount(components).flatMap: checked =>
      fit(
        MatrixView.dense(x),
        MatrixView.dense(y),
        checked,
        regularization,
        RegressionDirection.XToY,
        xPreproc,
        yPreproc
      )

  def fit(
      x: MatrixView,
      y: MatrixView,
      components: ComponentCount,
      regularization: RegressionRegularization = RegressionRegularization.Ols,
      direction: RegressionDirection = RegressionDirection.XToY,
      xPreproc: PreprocessSpec = PreprocessSpec.Center,
      yPreproc: PreprocessSpec = PreprocessSpec.Center,
      solver: SvdSolver = DenseSolvers.svd,
      eigenSolver: SymmetricEigenSolver = DenseSolvers.symmetricEigen,
      rowMetric: Option[MetricSpec] = None,
      policy: StoragePolicy = StoragePolicy.AllowDense
  ): Either[MultivarError, ReducedRankRegressionFit] =
    direction match
      case RegressionDirection.XToY =>
        fitXToY(x, y, components, regularization, xPreproc, yPreproc, solver, eigenSolver, rowMetric, policy)

  private def fitXToY(
      x: MatrixView,
      y: MatrixView,
      components: ComponentCount,
      regularization: RegressionRegularization,
      xPreproc: PreprocessSpec,
      yPreproc: PreprocessSpec,
      solver: SvdSolver,
      eigenSolver: SymmetricEigenSolver,
      rowMetric: Option[MetricSpec],
      policy: StoragePolicy
  ): Either[MultivarError, ReducedRankRegressionFit] =
    for
      _ <- validateRrrComponentRequest(components, x.rows, x.cols, y.cols)
      prepared <- PreparedPair.fromPredictive(x, y, xPreproc, yPreproc, rowMetric, policy)
      fittedY <- prepared.invertibleResponse
      problem <- PairedOperatorProblem.fromMatrices(
        prepared.xWorking,
        prepared.yWorking,
        prepared.rowMetric,
        "rrr",
        policy
      )
      operator <- problem.fitReducedRankRegression(
        components,
        regularization,
        Math.max(1, prepared.moments.sampleCount - 1).toDouble,
        solver,
        eigenSolver
      )
      sourceWeights <- operator.sourceWeights
      responseLoadings <- operator.targetWeights
      unconstrained <- operator.coefficient match
        case Some(value) => decompositionSemantic(value.toDense)
        case None        => Left(MultivarError.SolverFailed("RRR operator fit omitted its directed coefficient"))
      encoderWeights = MetricOperator.scaleColumnsDense(sourceWeights, operator.result.singularValues)
      working = GaleNumerics.multiply(encoderWeights, responseLoadings.transpose)
      (raw, intercept) <- rawCoordinateMap(working, prepared.xPreprocessor, fittedY)
      coefficientTransform <- FittedCoefficientTransform.from(working, prepared.xPreprocessor, fittedY, "rrr")
      sourceTransform <- FittedFrameTransform.fromTraining(
        prepared.xOriginal,
        encoderWeights,
        prepared.xPreprocessor,
        "rrr.source",
        components,
        Some(operator.result.singularValues)
      )
      targetTransform <- FittedFrameTransform.fromTraining(
        prepared.yOriginal,
        responseLoadings,
        fittedY,
        "rrr.target",
        components,
        Some(operator.result.singularValues)
      )
    yield ReducedRankRegressionFit(
      raw,
      intercept,
      working,
      unconstrained,
      coefficientTransform,
      sourceTransform,
      targetTransform,
      operator
    )

/** Convert a working-space coefficient map into original predictor/response units.
  *
  * With forward preprocessing `x ↦ x ⊙ D + a`, the identity
  * `y = ((x ⊙ D_x + a_x) B_w − a_y) ⊙ D_y^{-1}` rearranges to
  * `y = x B_raw + intercept` with `B_raw = D_x B_w D_y^{-1}` and
  * `intercept = (a_x B_w − a_y) ⊙ D_y^{-1}`.
  */
private def rawCoordinateMap(
    working: DMat,
    predictor: FittedPreprocessor,
    response: FittedInvertiblePreprocessor
): Either[MultivarError, (DMat, DVec)] =
  for
    predictorAffine <- columnAffine(predictor, "predictor")
    responseAffine <- response match
      case affine: InvertibleColumnAffine => Right(affine)
      case other =>
        Left(
          MultivarError.InvalidMap(
            s"response preprocessor must be a column affine to expose raw coefficients, got ${other.getClass.getName}"
          )
        )
  yield
    val dx = predictorAffine.scale
    val ax = predictorAffine.shift
    val dyInv = responseAffine.summary.scale
    val ay = responseAffine.forward.shift
    val rowScaled = MatrixView.scaleRows(working, dx)
    val raw = MetricOperator.scaleColumnsDense(rowScaled, dyInv)
    val intercept = interceptFrom(ax, working, ay, dyInv)
    (raw, intercept)

private def columnAffine(
    preprocessor: FittedPreprocessor,
    role: String
): Either[MultivarError, FittedColumnAffine] =
  preprocessor match
    case affine: FittedColumnAffine       => Right(affine)
    case affine: InvertibleColumnAffine   => Right(affine.forward)
    case other =>
      Left(
        MultivarError.InvalidMap(
          s"$role preprocessor must be a column affine to expose raw coefficients, got ${other.getClass.getName}"
        )
      )

private def interceptFrom(ax: DVec, working: DMat, ay: DVec, dyInv: DVec): DVec =
  val out = new Array[Double](working.cols)
  var col = 0
  while col < working.cols do
    var sum = 0.0
    var row = 0
    while row < working.rows do
      sum += ax(row) * working(row, col)
      row += 1
    out(col) = (sum - ay(col)) * dyInv(col)
    col += 1
  GaleNumerics.vectorFromArray(out)

private def validateRrrComponentRequest(
    components: ComponentCount,
    rows: Int,
    xCols: Int,
    yCols: Int
): Either[MultivarError, Unit] =
  val limit = Math.min(rows, Math.min(xCols, yCols))
  if components.value > limit then Left(MultivarError.InvalidComponentRequest(components.value, limit))
  else Right(())

private def decompositionSemantic[A](result: Either[SemanticError, A]): Either[MultivarError, A] =
  result.left.map:
    case SemanticError.MultivarFailure(error)  => error
    case SemanticError.LinearMapFailure(error) => LinalgErrorAdapter.toMultivarError(error)
    case error                                 => MultivarError.SolverFailed(error.message)

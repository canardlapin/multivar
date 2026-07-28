package multivar
package family.paired

import multivar.capability.{FittedCoefficientTransform, FittedFrameTransform}
import multivar.core.*

import gale.linalg.DMat
import gale.linalg.DVec
import gale.spectral.EigenOrder
import gale.spectral.EigenSelection

/** Algorithm used to extract PLS latent components. */
enum PlsAlgorithm:
  /** de Jong SIMPLS: successive SVD of the deflated cross-product `X'Y`. */
  case Simpls

  def label: String =
    this match
      case Simpls => "simpls"

final case class PlsRegressionDiagnostics(
    algorithm: PlsAlgorithm,
    requestedComponents: ComponentCount,
    effectiveComponents: Int
):
  require(effectiveComponents > 0, "PLS must retain at least one component")
  require(
    effectiveComponents <= requestedComponents.value,
    "effective PLS components cannot exceed the request"
  )

/** Opaque PLS regression result.
  *
  * [[coefficients]] / [[intercept]] are in original predictor/response units.
  * [[workingCoefficients]] is the map in preprocessed coordinates.
  * [[transform]] returns X-scores for new predictors after the fitted preprocessing.
  */
final class PlsRegressionFit private[multivar] (
    private val rawCoefficients: DMat,
    private val rawIntercept: DVec,
    private val workingMap: DMat,
    private val coefficientTransform: FittedCoefficientTransform,
    private val xScoreFrame: FittedFrameTransform,
    private val trainingYScores: DMat,
    private val xWeightMatrix: DMat,
    private val yWeightMatrix: DMat,
    private val xLoadingMatrix: DMat,
    private val yLoadingMatrix: DMat,
    private val xRotationMatrix: DMat,
    private val fitDiagnostics: PlsRegressionDiagnostics
):
  def coefficients: DMat = rawCoefficients
  def intercept: DVec = rawIntercept
  def workingCoefficients: DMat = workingMap

  def xScores: DMat = xScoreFrame.scores
  def yScores: DMat = trainingYScores
  def xWeights: DMat = xWeightMatrix
  def yWeights: DMat = yWeightMatrix
  def xLoadings: DMat = xLoadingMatrix
  def yLoadings: DMat = yLoadingMatrix

  /** Predictor map such that training scores satisfy `T ≈ X_working * xRotations`. */
  def xRotations: DMat = xRotationMatrix

  def requestedComponents: Int = fitDiagnostics.requestedComponents.value
  def effectiveComponents: Int = fitDiagnostics.effectiveComponents
  def diagnostics: PlsRegressionDiagnostics = fitDiagnostics

  def transform(input: DMat): Either[MultivarError, DMat] =
    transform(MatrixView.dense(input))

  def transform(input: MatrixView): Either[MultivarError, DMat] =
    xScoreFrame.project(input)

  def predict(input: DMat): Either[MultivarError, DMat] =
    predict(MatrixView.dense(input))

  def predict(input: MatrixView): Either[MultivarError, DMat] =
    coefficientTransform.predict(input)

  def predictWorking(input: DMat): Either[MultivarError, DMat] =
    predictWorking(MatrixView.dense(input))

  def predictWorking(input: MatrixView): Either[MultivarError, DMat] =
    coefficientTransform.predictWorking(input)

object PlsRegressionFit:
  private[multivar] def coefficientTransformOf(fit: PlsRegressionFit): FittedCoefficientTransform =
    fit.coefficientTransform

  private[multivar] def frameOf(fit: PlsRegressionFit): FittedFrameTransform =
    fit.xScoreFrame

object PlsRegression:
  def fit(
      x: DMat,
      y: DMat,
      components: Int
  ): Either[MultivarError, PlsRegressionFit] =
    fit(x, y, components, PreprocessSpec.Center, PreprocessSpec.Center, PlsAlgorithm.Simpls)

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      algorithm: PlsAlgorithm
  ): Either[MultivarError, PlsRegressionFit] =
    fit(x, y, components, PreprocessSpec.Center, PreprocessSpec.Center, algorithm)

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      xPreproc: PreprocessSpec,
      yPreproc: PreprocessSpec
  ): Either[MultivarError, PlsRegressionFit] =
    fit(x, y, components, xPreproc, yPreproc, PlsAlgorithm.Simpls)

  def fit(
      x: DMat,
      y: DMat,
      components: Int,
      xPreproc: PreprocessSpec,
      yPreproc: PreprocessSpec,
      algorithm: PlsAlgorithm
  ): Either[MultivarError, PlsRegressionFit] =
    ComponentCount(components).flatMap: checked =>
      fit(
        MatrixView.dense(x),
        MatrixView.dense(y),
        checked,
        xPreproc,
        yPreproc,
        algorithm
      )

  def fit(
      x: MatrixView,
      y: MatrixView,
      components: ComponentCount,
      xPreproc: PreprocessSpec = PreprocessSpec.Center,
      yPreproc: PreprocessSpec = PreprocessSpec.Center,
      algorithm: PlsAlgorithm = PlsAlgorithm.Simpls,
      rowMetric: Option[MetricSpec] = None,
      policy: StoragePolicy = StoragePolicy.AllowDense,
      eigenSolver: SymmetricEigenSolver = DenseSolvers.symmetricEigen
  ): Either[MultivarError, PlsRegressionFit] =
    for
      prepared <- PreparedPair.fromPredictive(x, y, xPreproc, yPreproc, rowMetric, policy)
      _ <- validateComponentRequest(components, prepared.moments)
      fit <-
        algorithm match
          case PlsAlgorithm.Simpls =>
            fitSimpls(prepared, components, eigenSolver, policy)
    yield fit

  private def validateComponentRequest(
      components: ComponentCount,
      moments: PairedMoments
  ): Either[MultivarError, Unit] =
    // SIMPLS deflates in predictor space; response width does not cap the
    // number of successive components (univariate multi-component PLS is valid).
    val limit = Math.min(moments.sampleCount, moments.xFeatures)
    if components.value > limit then Left(MultivarError.InvalidComponentRequest(components.value, limit))
    else Right(())

  private def fitSimpls(
      prepared: PreparedPair,
      components: ComponentCount,
      eigenSolver: SymmetricEigenSolver,
      policy: StoragePolicy
  ): Either[MultivarError, PlsRegressionFit] =
    for
      (x, y) <- prepared.workingDense(policy)
      fit <- fitSimplsDense(prepared, x, y, components, eigenSolver)
    yield fit

  /** SIMPLS matching R `pls::simpls.fit` on already-prepared working matrices. */
  private def fitSimplsDense(
      prepared: PreparedPair,
      x: DMat,
      y: DMat,
      components: ComponentCount,
      eigenSolver: SymmetricEigenSolver
  ): Either[MultivarError, PlsRegressionFit] =
    val n = x.rows
    val p = x.cols
    val m = y.cols
    val aMax = components.value
    val R = Array.ofDim[Double](p * aMax)
    val P = Array.ofDim[Double](p * aMax)
    val V = Array.ofDim[Double](p * aMax)
    val T = Array.ofDim[Double](n * aMax)
    val U = Array.ofDim[Double](n * aMax)
    val Qt = Array.ofDim[Double](aMax * m)
    var S = GaleNumerics.transposeMultiply(x, y)
    val sBuf = S.copyData
    var a = 0
    var failure = Option.empty[MultivarError]
    while a < aMax && failure.isEmpty do
      // Rebuild the live cross-product view from the reused buffer.
      S = GaleNumerics.matrixFromRowMajor(p, m, sBuf)
      computeQa(S, m, p, eigenSolver) match
        case Left(error) => failure = Some(error)
        case Right(qInit) =>
          var ra = matVec(S, qInit)
          var ta = matVec(x, ra)
          // Working matrices may already be centered; keep the R mean-correction
          // so a non-centered request stays consistent with pls::simpls.fit.
          centerVectorInPlace(ta)
          val tnorm = euclideanNorm(ta)
          if !(tnorm > 0.0 && tnorm.isFinite) then
            failure = Some(MultivarError.SolverFailed(s"SIMPLS component ${a + 1} produced a zero score"))
          else
            scaleVectorInPlace(ta, 1.0 / tnorm)
            scaleVectorInPlace(ra, 1.0 / tnorm)
            orientWeight(ra) match
              case (orientedR, flip) =>
                ra = orientedR
                if flip then scaleVectorInPlace(ta, -1.0)
            val pa = transposeMatVec(x, ta)
            val qa = transposeMatVec(y, ta)
            var va = pa.clone()
            if a > 0 then va = subtractProjection(va, V, p, a)
            val vnorm = euclideanNorm(va)
            if !(vnorm > 0.0 && vnorm.isFinite) then
              failure = Some(MultivarError.SolverFailed(s"SIMPLS component ${a + 1} produced a singular loading basis"))
            else
              scaleVectorInPlace(va, 1.0 / vnorm)
              deflateCrossInPlace(sBuf, p, m, va)
              copyColumn(R, p, a, ra)
              copyColumn(P, p, a, pa)
              copyColumn(V, p, a, va)
              copyColumn(T, n, a, ta)
              val ua = matVec(y, qa)
              val uOrtho =
                if a == 0 then ua
                else subtractProjection(ua, T, n, a)
              copyColumn(U, n, a, uOrtho)
              copyRow(Qt, m, a, qa)
              a += 1
    failure match
      case Some(error) => Left(error)
      case None =>
        val keep = a
        val rMat = matrixFromColumns(R, p, keep)
        val pMat = matrixFromColumns(P, p, keep)
        val uMat = matrixFromColumns(U, n, keep)
        val qMat = matrixFromRows(Qt, keep, m)
        val working = GaleNumerics.multiply(rMat, qMat)
        for
          response <- prepared.invertibleResponse
          (raw, intercept) <- rawCoordinateMap(working, prepared.xPreprocessor, response)
          coefficientTransform <- FittedCoefficientTransform.from(
            working,
            prepared.xPreprocessor,
            response,
            "pls"
          )
          xFrame <- FittedFrameTransform.fromTraining(
            prepared.xOriginal,
            rMat,
            prepared.xPreprocessor,
            "pls.source",
            components,
            None
          )
        yield
          new PlsRegressionFit(
            raw,
            intercept,
            working,
            coefficientTransform,
            xFrame,
            uMat,
            rMat,
            qMat.transpose,
            pMat,
            qMat.transpose,
            rMat,
            PlsRegressionDiagnostics(PlsAlgorithm.Simpls, components, keep)
          )

  private def computeQa(
      s: DMat,
      responses: Int,
      predictors: Int,
      eigenSolver: SymmetricEigenSolver
  ): Either[MultivarError, Array[Double]] =
    val top1 = EigenSelection.Count(1, EigenOrder.LargestAlgebraic)
    if responses == 1 then Right(Array(1.0))
    else if responses < predictors then
      // eigen(S'S): GaleNumerics.crossProduct(S) = S' S
      val sts = GaleNumerics.crossProduct(s)
      LinalgErrorAdapter.adapt(eigenSolver.decompose(MatrixOps.symmetrize(sts), top1)).map: eigen =>
        columnOf(eigen.vectors, 0)
    else
      // eigen(SS'): crossProduct(S') = S S'
      val sst = GaleNumerics.crossProduct(s.transpose)
      LinalgErrorAdapter.adapt(eigenSolver.decompose(MatrixOps.symmetrize(sst), top1)).flatMap: eigen =>
        val left = columnOf(eigen.vectors, 0)
        val qa = transposeMatVec(s, left)
        val norm = euclideanNorm(qa)
        if !(norm > 0.0 && norm.isFinite) then
          Left(MultivarError.SolverFailed("SIMPLS Y-weight is zero"))
        else Right(scaleVector(qa, 1.0 / norm))

  /** `S <- S - v (v' S)` into a reused row-major buffer. */
  private def deflateCrossInPlace(sBuf: Array[Double], rows: Int, cols: Int, v: Array[Double]): Unit =
    val vs = new Array[Double](cols)
    var col = 0
    while col < cols do
      var acc = 0.0
      var row = 0
      while row < rows do
        acc += v(row) * sBuf(row * cols + col)
        row += 1
      vs(col) = acc
      col += 1
    col = 0
    while col < cols do
      val scale = vs(col)
      var row = 0
      while row < rows do
        sBuf(row * cols + col) -= v(row) * scale
        row += 1
      col += 1

  private def subtractProjection(vector: Array[Double], basis: Array[Double], rows: Int, cols: Int): Array[Double] =
    val out = vector.clone()
    var col = 0
    while col < cols do
      var dot = 0.0
      var row = 0
      while row < rows do
        dot += basis(row + col * rows) * vector(row)
        row += 1
      row = 0
      while row < rows do
        out(row) -= basis(row + col * rows) * dot
        row += 1
      col += 1
    out

  private def orientWeight(weight: Array[Double]): (Array[Double], Boolean) =
    var best = 0
    var bestAbs = 0.0
    var i = 0
    while i < weight.length do
      val value = Math.abs(weight(i))
      if value > bestAbs then
        bestAbs = value
        best = i
      i += 1
    if weight(best) < 0.0 then (scaleVector(weight, -1.0), true) else (weight, false)

  private def matVec(matrix: DMat, vector: Array[Double]): Array[Double] =
    val out = new Array[Double](matrix.rows)
    var row = 0
    while row < matrix.rows do
      var acc = 0.0
      var col = 0
      while col < matrix.cols do
        acc += matrix(row, col) * vector(col)
        col += 1
      out(row) = acc
      row += 1
    out

  private def transposeMatVec(matrix: DMat, vector: Array[Double]): Array[Double] =
    val out = new Array[Double](matrix.cols)
    var col = 0
    while col < matrix.cols do
      var acc = 0.0
      var row = 0
      while row < matrix.rows do
        acc += matrix(row, col) * vector(row)
        row += 1
      out(col) = acc
      col += 1
    out

  private def centerVectorInPlace(values: Array[Double]): Unit =
    var sum = 0.0
    var i = 0
    while i < values.length do
      sum += values(i)
      i += 1
    val mean = sum / values.length.toDouble
    i = 0
    while i < values.length do
      values(i) -= mean
      i += 1

  private def scaleVectorInPlace(values: Array[Double], factor: Double): Unit =
    var i = 0
    while i < values.length do
      values(i) *= factor
      i += 1

  private def scaleVector(values: Array[Double], factor: Double): Array[Double] =
    val out = values.clone()
    scaleVectorInPlace(out, factor)
    out

  private def euclideanNorm(values: Array[Double]): Double =
    var acc = 0.0
    var i = 0
    while i < values.length do
      acc += values(i) * values(i)
      i += 1
    Math.sqrt(acc)

  private def columnOf(matrix: DMat, col: Int): Array[Double] =
    val out = new Array[Double](matrix.rows)
    var row = 0
    while row < matrix.rows do
      out(row) = matrix(row, col)
      row += 1
    out

  private def copyColumn(target: Array[Double], rows: Int, col: Int, values: Array[Double]): Unit =
    var row = 0
    while row < rows do
      target(row + col * rows) = values(row)
      row += 1

  private def copyRow(target: Array[Double], cols: Int, row: Int, values: Array[Double]): Unit =
    var col = 0
    while col < cols do
      target(row * cols + col) = values(col)
      col += 1

  private def matrixFromColumns(data: Array[Double], rows: Int, cols: Int): DMat =
    val out = new Array[Double](rows * cols)
    var col = 0
    while col < cols do
      var row = 0
      while row < rows do
        out(row * cols + col) = data(row + col * rows)
        row += 1
      col += 1
    GaleNumerics.matrixFromRowMajor(rows, cols, out)

  private def matrixFromRows(data: Array[Double], rows: Int, cols: Int): DMat =
    GaleNumerics.matrixFromRowMajor(rows, cols, java.util.Arrays.copyOf(data, rows * cols))

package multivar
package family.spectral

import gale.backend.Backend.given
import gale.linalg.{DMat, Matrix, QROptions, QRPivoting}
import gale.spectral.SpectralBackend.given
import multivar.core.{LinalgErrorAdapter, MatrixOps, MultivarError}

/** Row scaling used while optimizing a rotation criterion.  Kaiser scaling makes
  * each variable contribute equally to the varimax criterion; it does not change
  * the returned loading scale.
  */
enum RotationNormalization:
  case None
  case Kaiser

/** Bounded execution policy for a dense factor rotation.
  *
  * `maxElements` bounds the conservative sum of the full dense arrays owned by
  * this implementation at an allocation point (working matrices, retained
  * matrices and explicit small factors).  It does not make a claim about
  * Gale's private kernel workspace.
  */
final case class RotationOptions(
    normalization: RotationNormalization = RotationNormalization.Kaiser,
    tolerance: Double = 1e-7,
    maxIterations: Int = 256,
    maxElements: Long = 4_000_000L,
    maxCondition: Double = 1e8
)

/** Evidence from the local stationary-point search.  A converged result is not a
  * claim that the non-convex criterion's global optimum was found.
  */
final case class RotationReceipt(
    method: String,
    normalization: RotationNormalization,
    iterations: Int,
    initialCriterion: Double,
    finalCriterion: Double,
    criterionChange: Double,
    polarSum: Double,
    stoppingTolerance: Double,
    converged: Boolean
)

/** Promax keeps its least-squares target residual separate from varimax's
  * fourth-power criterion; they have different units and are not comparable.
  */
final case class PromaxReceipt(
    baseline: RotationReceipt,
    leastSquaresResidual: Double,
    normalization: RotationNormalization,
    conditionNumber: Double
)

/** An orthogonal factor rotation.  `loadings = input * transform` and
  * `inverseTransform` is retained as an actual matrix identity receipt.
  */
final case class VarimaxRotation(
    loadings: DMat,
    transform: DMat,
    inverseTransform: DMat,
    receipt: RotationReceipt
)

/** An oblique promax rotation.  Scores are assumed to have covariance I before
  * this rotation.  `scoreCorrelation` is therefore computed from the retained
  * inverse transform, not estimated from a data sample.
  */
final case class PromaxRotation(
    loadings: DMat,
    transform: DMat,
    inverseTransform: DMat,
    scoreCorrelation: DMat,
    conditionNumber: Double,
    varimax: VarimaxRotation,
    receipt: PromaxReceipt
)

object FactorRotation:
  def varimax(input: DMat, options: RotationOptions = RotationOptions()): Either[MultivarError, VarimaxRotation] =
    for
      _ <- validate(input, options, "varimax", varimaxBudget(input))
      scales <- rowScales(input, options.normalization)
      working = scaleRows(input, scales, inverse = true)
      result <- iterateVarimax(working, options)
      output = input * result._1
    yield VarimaxRotation(output, result._1, result._1.t, result._2)

  /** Promax is deliberately fixed to the common power-four target
    * `sign(V) * abs(V)^4`, where V is the Kaiser-normalized varimax solution.
    * It solves V Q approximately equal to that target by Gale's pivoted QR.
    */
  def promax(input: DMat, options: RotationOptions = RotationOptions()): Either[MultivarError, PromaxRotation] =
    for
      _ <- validate(input, options, "promax", promaxBudget(input))
      baseline <- varimax(input, options)
      scales <- rowScales(input, options.normalization)
      normalized = scaleRows(baseline.loadings, scales, inverse = true)
      target = powerFour(normalized)
      q <- LinalgErrorAdapter.adapt(
        normalized.qr(QROptions(pivoting = QRPivoting.Column)).solveLeastSquares(target)
      )
      rawTransform = baseline.transform * q
      normalizedTransform <- normalizeObliqueTransform(rawTransform, options.maxCondition)
      transform = normalizedTransform._1
      inverse = normalizedTransform._2
      correlation = normalizedTransform._3
      condition = normalizedTransform._4
      output = input * transform
      residual = promaxResidual(normalized, q, target)
      receipt = PromaxReceipt(baseline.receipt, residual, options.normalization, condition)
    yield PromaxRotation(output, transform, inverse, correlation, condition, baseline, receipt)

  private def iterateVarimax(working: DMat, options: RotationOptions): Either[MultivarError, (DMat, RotationReceipt)] =
    val p = working.rows
    val k = working.cols
    var rotation = DMat.eye(k)
    var criterion = varimaxCriterion(working * rotation)
    val initial = criterion
    var change = Double.PositiveInfinity
    var previousPolarSum = 0.0
    var continue = true
    var iteration = 0
    var failure = Option.empty[MultivarError]
    while iteration < options.maxIterations && continue && failure.isEmpty do
      val rotated = working * rotation
      val adjustment = varimaxAdjustment(rotated, p)
      // R is the global polar update from the original normalized loadings:
      // L'[(LR)^3 - (1/p) LR diag((LR)'LR)].  Replacing L' with (LR)'
      // incorrectly discards the accumulated rotation.
      fullSvd(working.t * adjustment) match
        case Left(error) => failure = Some(LinalgErrorAdapter.toMultivarError(error))
        case Right(svd) =>
          var next = svd.u * svd.vt
          var nextCriterion = varimaxCriterion(working * next)
          val criterionTolerance = options.tolerance * Math.max(1.0, Math.abs(criterion))
          val polarSum = singularSum(svd.singularValues)
          var fixedPoint = maximumCoordinateChange(next, rotation) <= options.tolerance
          // Full polar steps can alternate between equal-criterion rotations
          // (including ordinary two-row Kaiser inputs). Test the orthogonal
          // midpoint before declaring a plateau; equal criterion alone does
          // not establish a stationary representation.
          if !fixedPoint && nextCriterion <= criterion + criterionTolerance then
            fullSvd(rotation + next) match
              case Left(error) => failure = Some(LinalgErrorAdapter.toMultivarError(error))
              case Right(midpointSvd) =>
                val midpoint = midpointSvd.u * midpointSvd.vt
                val midpointCriterion = varimaxCriterion(working * midpoint)
                if midpointCriterion > nextCriterion then
                  next = midpoint
                  nextCriterion = midpointCriterion
                  fixedPoint = maximumCoordinateChange(next, rotation) <= options.tolerance
          if nextCriterion < criterion - criterionTolerance then
            failure = Some(MultivarError.NumericalResidualExceeded("varimax monotonic polar update", criterion - nextCriterion, 0.0))
          else
            // Objective rounding within the declared tolerance does not stop
            // the search. Only the actual polar fixed point does; accepting
            // such a step preserves monotonicity to that numerical tolerance.
            change = Math.abs(nextCriterion - criterion)
            rotation = next
            criterion = nextCriterion
            // An exactly zero polar update is stationary. In particular,
            // Kaiser-normalized rank-one rows are +/-1 and have zero update;
            // the relative test alone would keep 0 >= 0 true indefinitely.
            continue = !fixedPoint && polarSum > 0.0
            previousPolarSum = polarSum
          iteration += 1
    failure match
      case Some(error) => Left(error)
      case None if continue => Left(MultivarError.IterationLimitExceeded("varimax", options.maxIterations, change))
      case None =>
        Right(rotation -> RotationReceipt("varimax(gamma=1, svd-polar; plateau midpoint and coordinate fixed-point)", options.normalization, iteration, initial, criterion, change, previousPolarSum, options.tolerance, converged = true))

  private def maximumCoordinateChange(left: DMat, right: DMat): Double =
    var maximum = 0.0
    var row = 0
    while row < left.rows do
      var column = 0
      while column < left.cols do
        maximum = math.max(maximum, math.abs(left(row, column) - right(row, column)))
        column += 1
      row += 1
    maximum

  private def validate(input: DMat, options: RotationOptions, method: String, ownedElements: BigInt): Either[MultivarError, Unit] =
    if input.rows < 2 || input.cols < 1 then Left(MultivarError.MatrixShapeMismatch(s"$method requires at least two rows and one component, got ${input.rows}x${input.cols}"))
    else if input.rows < input.cols then Left(MultivarError.MatrixShapeMismatch(s"$method requires rows >= components for a full-column-rank criterion, got ${input.rows}x${input.cols}"))
    else if ownedElements > BigInt(options.maxElements) || ownedElements > BigInt(Int.MaxValue) then Left(MultivarError.DimensionOverflow(input.rows, input.cols))
    else if !options.tolerance.isFinite || options.tolerance <= 0.0 then Left(MultivarError.InvalidTolerance("rotation tolerance", options.tolerance))
    else if options.maxIterations <= 0 then Left(MultivarError.InvalidDimension("rotation iteration budget", options.maxIterations))
    else if !options.maxCondition.isFinite || options.maxCondition < 1.0 then Left(MultivarError.InvalidRegularization("rotation condition limit", options.maxCondition, "must be finite and at least one"))
    else
      for
        _ <- MatrixOps.checkFinite(s"$method input", input)
        svd <- fullSvd(input).left.map(LinalgErrorAdapter.toMultivarError)
        _ <- if svd.rank == input.cols && svd.singularValues.length == input.cols && svd.singularValues(input.cols - 1) > 0.0 then Right(()) else Left(MultivarError.SolverFailed(s"$method requires full column rank ${input.cols}"))
      yield ()

  private def rowScales(input: DMat, normalization: RotationNormalization): Either[MultivarError, Array[Double]] =
    val out = Array.fill(input.rows)(1.0)
    normalization match
      case RotationNormalization.None => Right(out)
      case RotationNormalization.Kaiser =>
        var row = 0
        var error = Option.empty[MultivarError]
        while row < input.rows && error.isEmpty do
          var sum = 0.0
          var col = 0
          while col < input.cols do
            sum += input(row, col) * input(row, col)
            col += 1
          if !sum.isFinite || sum <= 0.0 then error = Some(MultivarError.NonInvertibleValue("Kaiser row norm", row, sum))
          else out(row) = Math.sqrt(sum)
          row += 1
        error.toLeft(out)

  private def scaleRows(input: DMat, scales: Array[Double], inverse: Boolean): DMat =
    val out = Matrix.newBuilder(input.rows, input.cols)
    var row = 0
    while row < input.rows do
      val factor = if inverse then 1.0 / scales(row) else scales(row)
      var col = 0
      while col < input.cols do
        out(row, col) = input(row, col) * factor
        col += 1
      row += 1
    out.result()

  private def varimaxAdjustment(rotated: DMat, p: Int): DMat =
    val norms = new Array[Double](rotated.cols)
    var col = 0
    while col < rotated.cols do
      var row = 0
      while row < rotated.rows do
        val value = rotated(row, col)
        norms(col) += value * value
        row += 1
      col += 1
    val out = Matrix.newBuilder(rotated.rows, rotated.cols)
    var row = 0
    while row < rotated.rows do
      col = 0
      while col < rotated.cols do
        val value = rotated(row, col)
        out(row, col) = value * value * value - value * norms(col) / p.toDouble
        col += 1
      row += 1
    out.result()

  private def varimaxCriterion(rotated: DMat): Double =
    val p = rotated.rows.toDouble
    var total = 0.0
    var col = 0
    while col < rotated.cols do
      var sum2 = 0.0
      var sum4 = 0.0
      var row = 0
      while row < rotated.rows do
        val square = rotated(row, col) * rotated(row, col)
        sum2 += square
        sum4 += square * square
        row += 1
      total += sum4 - sum2 * sum2 / p
      col += 1
    total / p

  private def singularSum(values: gale.linalg.DVec): Double =
    var total = 0.0
    var index = 0
    while index < values.length do
      total += values(index)
      index += 1
    total

  private def powerFour(input: DMat): DMat =
    val out = Matrix.newBuilder(input.rows, input.cols)
    var row = 0
    while row < input.rows do
      var col = 0
      while col < input.cols do
        val value = input(row, col)
        out(row, col) = Math.copySign(value * value * value * value, value)
        col += 1
      row += 1
    out.result()

  private def promaxResidual(v: DMat, q: DMat, target: DMat): Double =
    val residual = v * q
    var sum = 0.0
    var row = 0
    while row < residual.rows do
      var col = 0
      while col < residual.cols do
        val d = residual(row, col) - target(row, col)
        sum += d * d
        col += 1
      row += 1
    Math.sqrt(sum)

  private def svdCondition(matrix: DMat, role: String): Either[MultivarError, Double] =
    fullSvd(matrix).left.map(LinalgErrorAdapter.toMultivarError).flatMap: svd =>
      val values = svd.singularValues
      if values.length == 0 || values(values.length - 1) <= 0.0 then Left(MultivarError.NonInvertibleValue(role, values.length - 1, if values.length == 0 then 0.0 else values(values.length - 1)))
      else Right(values(0) / values(values.length - 1))

  private def inverseBySvd(matrix: DMat, limit: Double, role: String): Either[MultivarError, DMat] =
    fullSvd(matrix).left.map(LinalgErrorAdapter.toMultivarError).flatMap: svd =>
      val values = svd.singularValues
      if values.length != matrix.rows || values.length != matrix.cols || values(values.length - 1) <= 0.0 then Left(MultivarError.NonInvertibleValue(role, values.length - 1, if values.length == 0 then 0.0 else values(values.length - 1)))
      else
        val condition = values(0) / values(values.length - 1)
        if !condition.isFinite || condition > limit then Left(MultivarError.SolverFailed(s"$role condition number $condition exceeds limit $limit"))
        else
          val inverseValues = new Array[Double](values.length)
          var i = 0
          while i < values.length do
            inverseValues(i) = 1.0 / values(i)
            i += 1
          Right(svd.vt.t * diagonal(inverseValues) * svd.u.t)

  /** Standardize the scores `z' = T^-1 z` under the explicit convention
    * cov(z) = I.  For raw T, cov(z') = T^-1 (T^-1)', so multiplying T on the
    * right by the raw score standard deviations makes the final covariance
    * have unit diagonal.
    */
  private[spectral] def normalizeObliqueTransform(
      rawTransform: DMat,
      maxCondition: Double
  ): Either[MultivarError, (DMat, DMat, DMat, Double)] =
    for
      rawInverse <- inverseBySvd(rawTransform, Double.PositiveInfinity, "promax raw transform")
      rawCovariance = rawInverse * rawInverse.t
      scoreStandardDeviations <- diagonalRoots(rawCovariance, "promax raw score covariance")
      transform = rawTransform * diagonal(scoreStandardDeviations)
      inverse <- inverseBySvd(transform, maxCondition, "promax normalized transform")
      covariance = inverse * inverse.t
      _ <- unitDiagonal(covariance, "promax normalized score covariance")
      condition <- svdCondition(transform, "promax transform")
      _ <- if condition <= maxCondition then Right(()) else Left(MultivarError.SolverFailed(s"promax transform condition number $condition exceeds limit $maxCondition"))
    yield (transform, inverse, covariance, condition)

  private def diagonalRoots(matrix: DMat, role: String): Either[MultivarError, Array[Double]] =
    val out = new Array[Double](matrix.rows)
    var i = 0
    var error = Option.empty[MultivarError]
    while i < matrix.rows && error.isEmpty do
      val value = matrix(i, i)
      if !value.isFinite || value <= 0.0 then error = Some(MultivarError.NonInvertibleValue(role, i, value))
      else out(i) = Math.sqrt(value)
      i += 1
    error.toLeft(out)

  private def unitDiagonal(matrix: DMat, role: String): Either[MultivarError, Unit] =
    var index = 0
    var error = Option.empty[MultivarError]
    while index < matrix.rows && error.isEmpty do
      if Math.abs(matrix(index, index) - 1.0) > 1e-8 then
        error = Some(MultivarError.NumericalResidualExceeded(role, Math.abs(matrix(index, index) - 1.0), 1e-8))
      index += 1
    error.toLeft(())

  private def fullSvd(matrix: DMat) =
    matrix.svd.flatMap(_.requireConverged)

  private def varimaxBudget(input: DMat): BigInt =
    val p = BigInt(input.rows)
    val k = BigInt(input.cols)
    // Working/rotated/adjustment/output matrices, normalization vectors and
    // validation factors; full and midpoint polar factors can coexist.
    // Returned Gale factors count here; its private solver workspace does not.
    12 * p * k + 24 * k * k

  private def promaxBudget(input: DMat): BigInt =
    val p = BigInt(input.rows)
    val k = BigInt(input.cols)
    // Includes the complete varimax baseline plus normalized target/QR
    // factors, residual output and raw/normalized inverse-transform factors.
    20 * p * k + 48 * k * k

  private def diagonal(values: Array[Double]): DMat =
    val out = Matrix.newBuilder(values.length, values.length)
    var i = 0
    while i < values.length do
      out(i, i) = values(i)
      i += 1
    out.result()

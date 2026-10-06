package multivar.inference

import gale.backend.Backend.given
import gale.linalg.{DMat, DVec, DenseDecompositions, Matrix, QROptions}
import gale.spectral.{SingularSelection, Svds}
import resample4s.kernel.Draw

/** Multiplicities of complete, equally weighted sampling units. */
final class ProfileBootstrapCounts private (val multiplicities: Vector[Int]):
  def units: Int = multiplicities.size

  /** Delete one complete copy without weakening the N-of-N bootstrap contract. */
  def deleteOne(unit: Int): Either[InferenceError, ProfileDeletionCounts] =
    ProfileDeletionCounts.fromBootstrap(this, unit)

  /** Stable mean of complete unit-major profiles under these exact counts.
    * Useful for draw-dependent linear readouts without an additional fit.
    */
  def meanProfiles(profiles: DMat, coordinatesPerUnit: Int,
      maxWorkingBytes: Long = 256L * 1024 * 1024,
      cancelled: () => Boolean = () => false): Either[InferenceError, DMat] =
    if coordinatesPerUnit < 1 || profiles.cols < 1 ||
        profiles.rows.toLong != units.toLong * coordinatesPerUnit then
      Left(InferenceError.InvalidReplicatePlan("readout profiles must contain every complete sampling unit"))
    else for
      _ <- ProfileBootstrapMath.allocation(BigInt(32) * coordinatesPerUnit * profiles.cols, maxWorkingBytes)
      _ <- ProfileBootstrapMath.finite(profiles, "bootstrap mean readout", cancelled)
      mean <- ProfileBootstrapMath.aggregate(profiles, this, coordinatesPerUnit, cancelled)
    yield mean

object ProfileBootstrapCounts:
  def fromIndices(units: Int, indices: IndexedSeq[Int]): Either[InferenceError, ProfileBootstrapCounts] =
    if units < 2 || indices.size != units then
      Left(InferenceError.InvalidReplicatePlan("a profile bootstrap draws N complete units from N >= 2 units"))
    else if indices.exists(i => i < 0 || i >= units) then
      Left(InferenceError.InvalidReplicatePlan("profile bootstrap index outside the sampling-unit domain"))
    else
      val counts = new Array[Int](units)
      indices.foreach(i => counts(i) += 1)
      Right(new ProfileBootstrapCounts(counts.toVector))

  def fromDraw(draw: Draw): Either[InferenceError, ProfileBootstrapCounts] =
    fromIndices(draw.domain, draw.toVector)

  def identity(units: Int): Either[InferenceError, ProfileBootstrapCounts] =
    fromIndices(units, Vector.range(0, Math.max(0, units)))

enum ProfileBootstrapMode:
  case Auto, Direct, Reduced

enum ProfileBootstrapPath:
  case Direct, QrReduced, QrDirectRefinement

final case class ProfileBootstrapOptions(
    mode: ProfileBootstrapMode = ProfileBootstrapMode.Auto,
    expectedReplicates: Int = 1000,
    relativeRankTolerance: Double = 1e-6,
    absoluteZeroTolerance: Double = 1e-12,
    accuracyTolerance: Double = 1e-10,
    maxWorkingBytes: Long = 1024L * 1024 * 1024
)

final case class ProfileBootstrapPreparation(
    path: ProfileBootstrapPath,
    estimatedWorkingBytes: Long,
    coordinateDimension: Int,
    maxProfileReconstructionError: Double,
    maxOrthogonalityError: Double
)

/** A fit in the complete profile span of its preparing core. */
sealed abstract class MeanProfileFit private[inference] (
    private[inference] val owner: MeanProfileBootstrapCore,
    val singularValues: DVec,
    val leftDirections: DMat,
    val loadingCoordinates: DMat,
    val path: ProfileBootstrapPath
):
  def rank: Int = singularValues.length

  def loadingBlock(from: Int, until: Int, scaled: Boolean = false): Either[InferenceError, DMat] =
    owner.liftBlock(coordinates(scaled), from, until)

  def coordinates(scaled: Boolean = false): DMat =
    if !scaled then loadingCoordinates
    else
      val out = Matrix.newBuilder(loadingCoordinates.rows, rank)
      for r <- 0 until loadingCoordinates.rows; c <- 0 until rank do
        out(r, c) = loadingCoordinates(r, c) * singularValues(c)
      out.result()

/** N-of-N fit. Its original count type and methods remain source compatible. */
final class MeanProfileBootstrapFit private[inference] (
    owner: MeanProfileBootstrapCore,
    val counts: ProfileBootstrapCounts,
    singularValues: DVec,
    leftDirections: DMat,
    loadingCoordinates: DMat,
    path: ProfileBootstrapPath
) extends MeanProfileFit(owner, singularValues, leftDirections, loadingCoordinates, path)

/** N-1 fit over the same source-unit domain, deleting exactly one copy. */
final class MeanProfileDeletionFit private[inference] (
    owner: MeanProfileBootstrapCore,
    val counts: ProfileDeletionCounts,
    singularValues: DVec,
    leftDirections: DMat,
    loadingCoordinates: DMat,
    path: ProfileBootstrapPath
) extends MeanProfileFit(owner, singularValues, leftDirections, loadingCoordinates, path)

/** A fixed linear readout in the same feature geometry as a prepared core. */
final class ProfileBootstrapProjection private[inference] (
    private val owner: MeanProfileBootstrapCore,
    val coordinates: DMat
):
  def scores(fit: MeanProfileBootstrapFit, scaled: Boolean = false): Either[InferenceError, DMat] =
    scoreFit(fit, scaled)

  def scores(fit: MeanProfileDeletionFit): Either[InferenceError, DMat] = scoreFit(fit, false)

  def scores(fit: MeanProfileDeletionFit, scaled: Boolean): Either[InferenceError, DMat] = scoreFit(fit, scaled)

  private def scoreFit(fit: MeanProfileFit, scaled: Boolean): Either[InferenceError, DMat] =
    if !(fit.owner eq owner) then Left(InferenceError.InvalidReplicatePlan("projection and fit use different bootstrap cores"))
    else
      val values = coordinates * fit.coordinates(scaled)
      ProfileBootstrapMath.finite(values, "bootstrap projected scores").map(_ => values)

/** Bootstrap refits of the mean of N fixed d-by-P profiles, stored unit-major.
  * General linear algebra is delegated to Gale. No observed-rank truncation or
  * profile Gram cache is used in preparation.
  */
final class MeanProfileBootstrapCore private[inference] (
    private val profiles: DMat,
    val coordinatesPerUnit: Int,
    val featureCount: Int,
    private val reduced: DMat,
    private val basis: Option[DMat],
    private val sourceScale: Double,
    private val rowErrors: Vector[Double],
    private val rowNorms: Vector[Double],
    val options: ProfileBootstrapOptions,
    val preparation: ProfileBootstrapPreparation
):
  import ProfileBootstrapMath.*
  val unitCount: Int = profiles.rows / coordinatesPerUnit
  val coordinateDimension: Int = reduced.cols

  def observed(cancelled: () => Boolean = () => false): Either[InferenceError, MeanProfileBootstrapFit] =
    ProfileBootstrapCounts.identity(unitCount).flatMap(refit(_, cancelled))

  def refit(
      counts: ProfileBootstrapCounts,
      cancelled: () => Boolean = () => false
  ): Either[InferenceError, MeanProfileBootstrapFit] =
    refitCounts(counts.multiplicities, counts.units, cancelled).map { (s, u, v, path) =>
      new MeanProfileBootstrapFit(this, counts, s, u, v, path)
    }

  def refitDeletion(counts: ProfileDeletionCounts,
      cancelled: () => Boolean = () => false): Either[InferenceError, MeanProfileDeletionFit] =
    refitCounts(counts.multiplicities, counts.sampleSize, cancelled).map { (s, u, v, path) =>
      new MeanProfileDeletionFit(this, counts, s, u, v, path)
    }

  private def refitCounts(multiplicities: Vector[Int], sampleSize: Int, cancelled: () => Boolean)
      : Either[InferenceError, (DVec, DMat, DMat, ProfileBootstrapPath)] =
    if multiplicities.size != unitCount then
      Left(InferenceError.RowCountMismatch("bootstrap sampling units", unitCount, multiplicities.size))
    else if cancelled() then cancellation
    else aggregate(reduced, multiplicities, sampleSize, coordinatesPerUnit, cancelled).flatMap { mean =>
      val meanNorm = norm(mean)
      var errorBound = 0.0
      var coordinate = 0
      while coordinate < coordinatesPerUnit do
        var rowBound = 0.0
        var unit = 0
        while unit < unitCount do
          val row = unit * coordinatesPerUnit + coordinate
          rowBound += multiplicities(unit).toDouble / sampleSize *
            (rowErrors(row) + 8.0 * Math.ulp(1.0) * unitCount * rowNorms(row))
          unit += 1
        errorBound = Math.hypot(errorBound, rowBound)
        coordinate += 1
      val guarded = basis match
        case Some(q) if errorBound > options.accuracyTolerance * meanNorm =>
          // Cancellation can expose details below the QR's global residual scale.
          // Recompute the mean from the original profile entries, then verify its
          // representation at the *mean's* scale before using the small solve.
          for
            rawMean <- aggregate(profiles, multiplicities, sampleSize, coordinatesPerUnit, cancelled)
            direct <- checkedDivision(rawMean, sourceScale, "profile QR refinement")
            corrected = direct * q
            error <- reconstructionError(direct, corrected, q, cancelled)
            allowance = options.accuracyTolerance * norm(direct)
            result <-
              if error > allowance then
                Left(InferenceError.NumericalFailure("profile QR refinement", s"mean residual $error exceeds $allowance"))
              else Right((corrected, ProfileBootstrapPath.QrDirectRefinement))
          yield result
        case _ => Right((mean, preparation.path))
      guarded.flatMap { (matrix, path) =>
        if cancelled() then cancellation
        else decompose(matrix, sourceScale, options).map { (s, u, v) =>
          (s, u, v, path)
        }
      }
    }

  /** Cache a full or contiguous-feature linear readout. Lag-specific readouts
    * can be prepared separately; no temporal interpretation is imposed here.
    */
  def project(
      readout: DMat,
      featureFrom: Int = 0,
      cancelled: () => Boolean = () => false
  ): Either[InferenceError, ProfileBootstrapProjection] =
    if featureFrom < 0 || readout.cols < 1 || featureFrom.toLong + readout.cols > featureCount then
      Left(InferenceError.InvalidReplicatePlan("linear readout is outside the profile feature domain"))
    else if cancelled() then cancellation
    else for
      _ <- finite(readout, "bootstrap readout")
      _ <- allocation(BigInt(readout.rows) * coordinateDimension * 8, options.maxWorkingBytes)
      coordinates = basis match
        case Some(q) => readout * q.slice(featureFrom, featureFrom + readout.cols, 0, q.cols)
        case None =>
          val out = Matrix.newBuilder(readout.rows, featureCount)
          for r <- 0 until readout.rows; c <- 0 until readout.cols do out(r, featureFrom + c) = readout(r, c)
          out.result()
      _ <- finite(coordinates, "bootstrap cached readout", cancelled)
    yield new ProfileBootstrapProjection(this, coordinates)

  /** Lift only the requested feature rows. Coordinates must belong to this core. */
  def liftBlock(values: DMat, from: Int, until: Int): Either[InferenceError, DMat] =
    if values.rows != coordinateDimension then
      Left(InferenceError.RowCountMismatch("bootstrap loading coordinates", coordinateDimension, values.rows))
    else if from < 0 || until <= from || until > featureCount then
      Left(InferenceError.InvalidReplicatePlan("invalid bootstrap feature block"))
    else for
      _ <- finite(values, "bootstrap loading coordinates")
      _ <- allocation(BigInt(until - from) * values.cols * 8, options.maxWorkingBytes)
      lifted = basis match
        case Some(q) => q.slice(from, until, 0, q.cols) * values
        case None => values.slice(from, until, 0, values.cols)
      _ <- finite(lifted, "bootstrap lifted loadings")
    yield lifted

  /** Basis rows for exact reduced covariance or retained-draw recovery. */
  def basisBlock(from: Int, until: Int): Either[InferenceError, DMat] =
    if from < 0 || until <= from || until > featureCount then
      Left(InferenceError.InvalidReplicatePlan("invalid bootstrap basis block"))
    else allocation(BigInt(until - from) * coordinateDimension * 8, options.maxWorkingBytes).map { _ =>
      basis match
        case Some(q) => q.slice(from, until, 0, q.cols)
        case None =>
          val out = Matrix.newBuilder(until - from, coordinateDimension)
          for row <- from until until do out(row - from, row) = 1.0
          out.result()
    }

object MeanProfileBootstrap:
  import ProfileBootstrapMath.*

  def prepare(
      profiles: DMat,
      coordinatesPerUnit: Int,
      options: ProfileBootstrapOptions = ProfileBootstrapOptions(),
      cancelled: () => Boolean = () => false
  ): Either[InferenceError, MeanProfileBootstrapCore] =
    if coordinatesPerUnit < 1 || profiles.cols < 1 || profiles.rows % coordinatesPerUnit != 0 ||
        profiles.rows / coordinatesPerUnit < 2 then
      return Left(InferenceError.InvalidReplicatePlan("profiles require N >= 2 complete units of d >= 1 rows and P >= 1 columns"))
    if options.expectedReplicates < 1 then return Left(InferenceError.InvalidCount("expected bootstrap replicates", options.expectedReplicates))
    if !options.relativeRankTolerance.isFinite || options.relativeRankTolerance <= 0 || options.relativeRankTolerance >= 1 then
      return Left(InferenceError.InvalidTolerance("profile relative singular-value cutoff", options.relativeRankTolerance))
    if !options.absoluteZeroTolerance.isFinite || options.absoluteZeroTolerance < 0 then
      return Left(InferenceError.InvalidTolerance("profile absolute zero", options.absoluteZeroTolerance))
    if !options.accuracyTolerance.isFinite || options.accuracyTolerance <= 0 || options.accuracyTolerance >= options.relativeRankTolerance then
      return Left(InferenceError.InvalidTolerance("profile numerical accuracy", options.accuracyTolerance))
    if cancelled() then return cancellation
    val m = profiles.rows
    val p = profiles.cols
    val d = coordinatesPerUnit
    val r = Math.min(m, p)
    // Conservative live-array estimate, including caller input and QR scratch.
    val reducedBytes = 8 * (BigInt(9) * m * p + BigInt(4) * r * r + BigInt(5) * d * r)
    val directBytes = 8 * (BigInt(2) * m * p + BigInt(8) * d * p)
    val directWork = BigInt(options.expectedReplicates) * p * (m + d * d)
    val reducedWork = BigInt(p) * m * m * (2 + d) + BigInt(options.expectedReplicates) * r * (m + d * d)
    val useReduced = options.mode match
      case ProfileBootstrapMode.Direct => false
      case ProfileBootstrapMode.Reduced => true
      case ProfileBootstrapMode.Auto => p > m && reducedWork < directWork && reducedBytes <= options.maxWorkingBytes
    val bytes = if useReduced then reducedBytes else directBytes
    for
      _ <- allocation(bytes, options.maxWorkingBytes)
      _ <- finite(profiles, "bootstrap profile", cancelled)
      core <-
        var maximum = 0.0
        for row <- 0 until m; col <- 0 until p do maximum = Math.max(maximum, Math.abs(profiles(row, col)))
        val scale = if maximum == 0.0 then 1.0 else maximum
        if useReduced then
          checkedDivision(profiles, scale, "profile QR normalization").flatMap { normalized =>
            val qr = DenseDecompositions.qr(normalized.t, QROptions())
            // QR.q is P-by-P. Applying Q to the first r coordinate vectors builds
            // only its thin P-by-r part through Gale's public API.
            val firstColumns = Matrix.newBuilder(p, r)
            for i <- 0 until r do firstColumns(i, i) = 1.0
            qr.applyQ(firstColumns.result()).left.map(e => InferenceError.NumericalFailure("profile thin QR", e.getMessage)).flatMap { q =>
              val coreBuilder = Matrix.newBuilder(m, r)
              for row <- 0 until m; col <- 0 until r do coreBuilder(row, col) = qr.r(col, row)
              val c = coreBuilder.result()
              val orthogonal = q.t * q
              var orthogonality = 0.0
              for i <- 0 until r; j <- 0 until r do
                orthogonality = Math.max(orthogonality, Math.abs(orthogonal(i, j) - (if i == j then 1.0 else 0.0)))
              if orthogonality > options.accuracyTolerance then
                Left(InferenceError.NumericalFailure("profile thin QR", s"orthogonality error $orthogonality exceeds tolerance"))
              else
                rowReconstructionErrors(normalized, c, q, cancelled).map { rowErrors =>
                    new MeanProfileBootstrapCore(profiles, d, p, c, Some(q), scale, rowErrors,
                      Vector.tabulate(m)(row => norm(c.slice(row, row + 1, 0, r))), options,
                      ProfileBootstrapPreparation(ProfileBootstrapPath.QrReduced, bytes.toLong, r, rowErrors.max, orthogonality))
                }
            }
          }
        else Right(new MeanProfileBootstrapCore(profiles, d, p, profiles, None, 1.0,
          Vector.fill(m)(0.0), Vector.fill(m)(0.0), options,
          ProfileBootstrapPreparation(ProfileBootstrapPath.Direct, bytes.toLong, p, 0.0, 0.0)))
    yield core

private[inference] object ProfileBootstrapMath:
  def cancellation[A]: Either[InferenceError, A] = Left(InferenceError.NumericalFailure("profile bootstrap", "cancelled"))

  def allocation(bytes: BigInt, limit: Long): Either[InferenceError, Unit] =
    if limit <= 0 || bytes > limit || bytes / 8 > Int.MaxValue then
      Left(InferenceError.UnsupportedProblem(s"bootstrap working allocation $bytes bytes exceeds limit $limit or array capacity"))
    else Right(())

  def finite(matrix: DMat, role: String, cancelled: () => Boolean = () => false): Either[InferenceError, Unit] =
    var row = 0
    while row < matrix.rows do
      if cancelled() then return cancellation
      var col = 0
      while col < matrix.cols do
        if !matrix(row, col).isFinite then return Left(InferenceError.NonFiniteStatistic(role, matrix(row, col)))
        col += 1
      row += 1
    Right(())

  def norm(matrix: DMat): Double =
    var result = 0.0
    for row <- 0 until matrix.rows; col <- 0 until matrix.cols do result = Math.hypot(result, matrix(row, col))
    result

  def scaled(matrix: DMat, scale: Double): DMat =
    val out = Matrix.newBuilder(matrix.rows, matrix.cols)
    for row <- 0 until matrix.rows; col <- 0 until matrix.cols do out(row, col) = matrix(row, col) * scale
    out.result()

  def divided(matrix: DMat, scale: Double): DMat =
    val out = Matrix.newBuilder(matrix.rows, matrix.cols)
    for row <- 0 until matrix.rows; col <- 0 until matrix.cols do out(row, col) = matrix(row, col) / scale
    out.result()

  def checkedDivision(matrix: DMat, scale: Double, role: String): Either[InferenceError, DMat] =
    val out = divided(matrix, scale)
    var row = 0
    while row < matrix.rows do
      var col = 0
      while col < matrix.cols do
        if matrix(row, col) != 0.0 && out(row, col) == 0.0 then
          return Left(InferenceError.NumericalFailure(role, "scaling would erase a nonzero profile entry"))
        col += 1
      row += 1
    Right(out)

  def aggregate(matrix: DMat, counts: ProfileBootstrapCounts, d: Int,
      cancelled: () => Boolean = () => false): Either[InferenceError, DMat] =
    aggregate(matrix, counts.multiplicities, counts.units, d, cancelled)

  def aggregate(matrix: DMat, multiplicities: Vector[Int], sampleSize: Int, d: Int,
      cancelled: () => Boolean): Either[InferenceError, DMat] =
    val out = Matrix.newBuilder(d, matrix.cols)
    var coordinate = 0
    while coordinate < d do
      var col = 0
      while col < matrix.cols do
        if cancelled() then return cancellation
        var sum = 0.0
        var correction = 0.0
        var unit = 0
        while unit < multiplicities.size do
          val term = matrix(unit * d + coordinate, col)
          var copy = 0
          // Integer multiplicities describe complete copies. Dividing each term
          // first introduces product rounding that compensated summation cannot
          // recover when large profiles cancel (for example, thirds of 1e8).
          while copy < multiplicities(unit) do
            val next = sum + term
            correction += (if Math.abs(sum) >= Math.abs(term) then (sum - next) + term else (term - next) + sum)
            sum = next
            copy += 1
          unit += 1
        val total = sum + correction
        val average = total / sampleSize
        if !average.isFinite || (total != 0.0 && average == 0.0) then
          return Left(InferenceError.NumericalFailure("profile aggregation", "mean accumulation exceeded floating-point range"))
        out(coordinate, col) = average
        col += 1
      coordinate += 1
    Right(out.result())

  def reconstructionError(original: DMat, coordinates: DMat, q: DMat, cancelled: () => Boolean): Either[InferenceError, Double] =
    rowReconstructionErrors(original, coordinates, q, cancelled).map(_.foldLeft(0.0)(Math.hypot))

  def rowReconstructionErrors(original: DMat, coordinates: DMat, q: DMat, cancelled: () => Boolean): Either[InferenceError, Vector[Double]] =
    val errors = new Array[Double](original.rows)
    var from = 0
    while from < original.cols do
      if cancelled() then return cancellation
      val until = Math.min(original.cols, from + 256)
      val block = coordinates * q.slice(from, until, 0, q.cols).t
      for row <- 0 until original.rows; col <- from until until do
        errors(row) = Math.hypot(errors(row), original(row, col) - block(row, col - from))
      from = until
    errors.find(e => !e.isFinite) match
      case Some(error) => Left(InferenceError.NonFiniteStatistic("QR reconstruction error", error))
      case None => Right(errors.toVector)

  def decompose(matrix: DMat, scale: Double, options: ProfileBootstrapOptions): Either[InferenceError, (DVec, DMat, DMat)] =
    val size = norm(matrix)
    val magnitude = size * scale
    if !magnitude.isFinite then Left(InferenceError.NonFiniteStatistic("bootstrap mean norm", magnitude))
    else if size > 0.0 && magnitude == 0.0 && options.absoluteZeroTolerance == 0.0 then
      Left(InferenceError.NumericalFailure("bootstrap mean norm", "nonzero magnitude is below floating-point range"))
    else if magnitude <= options.absoluteZeroTolerance then
      Right((InferenceNumerics.vectorFromSeq(Vector.empty), Matrix.zeros(matrix.rows, 0), Matrix.zeros(matrix.cols, 0)))
    else if matrix.rows == 1 then
      Right((InferenceNumerics.vectorFromSeq(Vector(magnitude)), Matrix(1, 1)(1.0), divided(matrix.t, size)))
    else
      Svds.svd(divided(matrix, size), SingularSelection.All)
        .flatMap(_.requireConverged).left.map(e => InferenceError.NumericalFailure("bootstrap reduced SVD", e.getMessage))
        .flatMap { svd =>
          val largest = svd.singularValues(0)
          val ambiguous = (1 until svd.size).exists { i =>
            Math.abs(svd.singularValues(i) / largest - options.relativeRankTolerance) <= options.accuracyTolerance
          }
          if ambiguous then Left(InferenceError.NumericalFailure("bootstrap rank", "singular value lies within numerical uncertainty of the rank cutoff"))
          else
            val rank = (0 until svd.size).count(i => svd.singularValues(i) > options.relativeRankTolerance * largest)
            Right((InferenceNumerics.vectorFromSeq(Vector.tabulate(rank)(i => svd.singularValues(i) * magnitude)),
              svd.u.slice(0, matrix.rows, 0, rank), svd.vt.slice(0, rank, 0, matrix.cols).t))
        }

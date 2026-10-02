package multivar.inference

import gale.backend.Backend.given
import gale.linalg.{DMat, QROptions, QRPivoting}
import gale.spectral.SpectralBackend.given
import gale.spectral.SpectralDiagnostics
import multivar.core.{ComponentCount, MatrixOps, MatrixView, PreprocessSpec}
import multivar.family.paired.{Cca, CcaFit}
import resample4s.kernel.{Permutation, Seed}

/** A compact, full-rank classical CCA problem in already prepared coordinates.
  * Lower observed modes are removed before each randomized higher-rank fit.
  * Complements of the coefficient matrices preserve unequal candidate spaces.
  * This numerical procedure does not itself qualify a scientific error law. */
final class StepwiseCanonicalRank private (
    val leftVariables: DMat, val rightVariables: DMat,
    val correlations: Vector[Double], val observedWilks: Vector[Double],
    val plannedElements: Long, val observedCanonicalDiagnostics: SpectralDiagnostics,
    private val stepBases: Vector[(DMat, DMat)]
):
  def candidateRank: Int = correlations.size
  val qrRankTolerance: Double = StepwiseCanonicalRank.qrRankTolerance

  def nullStatistics(permutation: Permutation): Either[InferenceError, Vector[Double]] =
    nullStatisticsWithDiagnostics(permutation).map(_.map(_._1))

  private[inference] def nullStatisticsWithDiagnostics(permutation: Permutation): Either[InferenceError, Vector[(Double, SpectralDiagnostics)]] =
    if permutation.domain != leftVariables.rows then Left(InferenceError.RowCountMismatch("canonical residual action", leftVariables.rows, permutation.domain))
    else
      val rows = permutation.toIArray
      val output = Vector.newBuilder[(Double, SpectralDiagnostics)]
      var step = 0
      while step < candidateRank do
        val (left, right) = stepBases(step)
        StepwiseCanonicalRank.crossFit(left, right, i => rows(i)) match
          case Left(error) => return Left(error)
          case Right(value) => output += value
        step += 1
      Right(output.result())

object StepwiseCanonicalRank:
  val qrRankTolerance: Double = 1e-12
  def from(left: DMat, right: DMat, maximumElements: Long = 4_000_000L): Either[InferenceError, StepwiseCanonicalRank] =
    val n = BigInt(left.rows); val p = BigInt(left.cols); val q = BigInt(right.cols)
    val k = p.min(q)
    val retainedBasisColumns = k * (p + q) - k * (k - 1)
    val required = 12 * n * (p + q) + 4 * n * retainedBasisColumns + 40 * (p * p + q * q + p * q)
    if left.rows != right.rows then Left(InferenceError.RowCountMismatch("canonical paired coordinates", left.rows, right.rows))
    else if left.cols < 1 || right.cols < 1 || left.rows <= left.cols.toLong + right.cols then Left(InferenceError.InvalidPartition("classical rank confirmation requires rows > left columns + right columns"))
    else if maximumElements < 0L || required > maximumElements || Vector(n * p, n * q, p * p, q * q).exists(_ > Int.MaxValue) then Left(InferenceError.UnsupportedProblem(s"stepwise CCA requires $required owned elements, allowed $maximumElements"))
    else
      for
        _ <- finite("left canonical input", left)
        _ <- finite("right canonical input", right)
        _ <- fullRank(left)
        _ <- fullRank(right)
        canonical <- completeClassical(left, right)
        (roots, leftWeights, rightWeights, observedDiagnostics) = canonical
        _ <- wilks(roots).map(_ => ())
        leftCoefficients <- complete(leftWeights)
        rightCoefficients <- complete(rightWeights)
        u = left * leftCoefficients
        v = right * rightCoefficients
        _ <- finite("complete left canonical variables", u)
        _ <- finite("complete right canonical variables", v)
        bases <-
          val out = Vector.newBuilder[(DMat, DMat)]
          var step = 0
          var failure = Option.empty[InferenceError]
          while step < roots.size && failure.isEmpty do
            val x = DMat.tabulate(u.rows, u.cols - step)((i, j) => u(i, j + step))
            val y = DMat.tabulate(v.rows, v.cols - step)((i, j) => v(i, j + step))
            (for
              leftBasis <- orthonormal(x)
              rightBasis <- orthonormal(y)
            yield (leftBasis, rightBasis)) match
              case Left(error) => failure = Some(error)
              case Right(value) => out += value
            step += 1
          failure.toLeft(out.result())
        observed <-
          val out = Vector.newBuilder[Double]
          var index = 0
          var failure = Option.empty[InferenceError]
          while index < bases.size && failure.isEmpty do
            crossWilks(bases(index)._1, bases(index)._2, identity) match
              case Left(error) => failure = Some(error)
              case Right(value) => out += value
            index += 1
          failure.toLeft(out.result())
      yield new StepwiseCanonicalRank(u, v, roots, observed, required.toLong, observedDiagnostics, bases)

  /** A row permutation is orthogonal, so each side's Gram matrix is unchanged.
    * CCA on permuted tail variables is therefore exactly the converged SVD of
    * (P Q_left)^T Q_right. Cache the step-specific QR bases, not a randomized
    * statistic or a representation selected from confirmation permutations.
    * Every draw still fits a new compact cross-basis SVD at every rank. */
  private[inference] def crossWilks(left: DMat, right: DMat, sourceRow: Int => Int): Either[InferenceError, Double] =
    crossFit(left, right, sourceRow).map(_._1)

  private[inference] def crossFit(left: DMat, right: DMat, sourceRow: Int => Int): Either[InferenceError, (Double, SpectralDiagnostics)] =
    val cross = DMat.tabulate(left.cols, right.cols): (l, r) =>
      var total = 0.0
      var row = 0
      while row < left.rows do
        total += left(sourceRow(row), l) * right(row, r)
        row += 1
      total
    cross.svd.flatMap(_.requireConverged).left.map(error => InferenceError.NumericalFailure("stepwise compact CCA SVD", error.toString)).flatMap(result => wilks(result.singularValues.toVector).map(value => (value, result.diagnostics)))

  private def orthonormal(input: DMat): Either[InferenceError, DMat] =
    val qr = input.qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(qrRankTolerance)))
    val rank = qr.diagnostics.rank.getOrElse(0)
    if rank != input.cols then Left(InferenceError.RankLoss(input.cols, rank))
    else
      val selectors = DMat.tabulate(input.rows, input.cols)((i, j) => if i == j then 1.0 else 0.0)
      qr.applyQ(selectors).left.map(error => InferenceError.NumericalFailure("stepwise canonical QR basis", error.toString))

  /** The fitted latent-model CCA surface truncates numerical zero roots.
    * Confirmation must retain all K hypotheses, including zero roots. The
    * full thin Gale SVD supplies their coordinate vectors without a local
    * eigensolver or covariance inverse. */
  private def completeClassical(left: DMat, right: DMat): Either[InferenceError, (Vector[Double], DMat, DMat, SpectralDiagnostics)] =
    val leftQr = left.qr(QROptions(QRPivoting.Column, Some(qrRankTolerance)))
    val rightQr = right.qr(QROptions(QRPivoting.Column, Some(qrRankTolerance)))
    def numerical[A](role: String, value: Either[gale.linalg.LinAlgError, A]): Either[InferenceError, A] =
      value.left.map(error => InferenceError.NumericalFailure(role, error.toString))
    for
      leftBasis <- numerical("observed left canonical basis", leftQr.applyQ(DMat.tabulate(left.rows, left.cols)((i, j) => if i == j then 1.0 else 0.0)))
      rightBasis <- numerical("observed right canonical basis", rightQr.applyQ(DMat.tabulate(right.rows, right.cols)((i, j) => if i == j then 1.0 else 0.0)))
      spectral <- numerical("complete observed canonical SVD", (leftBasis.t * rightBasis).svd.flatMap(_.requireConverged))
      _ <- if spectral.singularValues.length == math.min(left.cols, right.cols) then Right(()) else Left(InferenceError.InvalidSpectrum("complete observed CCA omitted candidate roots"))
      leftWeights <- numerical("observed left canonical weights", leftQr.solveLeastSquares(leftBasis * spectral.u))
      rightWeights <- numerical("observed right canonical weights", rightQr.solveLeastSquares(rightBasis * spectral.vt.t))
    yield (spectral.singularValues.toVector, leftWeights, rightWeights, spectral.diagnostics)

  private[inference] def classical(left: DMat, right: DMat): Either[InferenceError, CcaFit] =
    for
      count <- ComponentCount(math.min(left.cols, right.cols)).left.map(error => InferenceError.NumericalFailure("classical CCA count", error.message))
      fit <- Cca.fit(MatrixView.dense(left), MatrixView.dense(right), count, ridge = 0.0,
        xPreproc = PreprocessSpec.Pass, yPreproc = PreprocessSpec.Pass).left.map(error => InferenceError.NumericalFailure("classical CCA", error.message))
    yield fit

  private def complete(coefficients: DMat): Either[InferenceError, DMat] =
    if coefficients.rows == coefficients.cols then Right(coefficients)
    else
      val qr = coefficients.qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(qrRankTolerance)))
      val rank = qr.diagnostics.rank.getOrElse(0)
      if rank != coefficients.cols then Left(InferenceError.RankLoss(coefficients.cols, rank))
      else
        val count = coefficients.rows - coefficients.cols
        val selectors = DMat.tabulate(coefficients.rows, count)((i, j) => if i == j + coefficients.cols then 1.0 else 0.0)
        qr.applyQ(selectors).left.map(error => InferenceError.NumericalFailure("canonical coefficient complement", error.toString)).map: complement =>
          DMat.tabulate(coefficients.rows, coefficients.rows)((i, j) => if j < coefficients.cols then coefficients(i, j) else complement(i, j - coefficients.cols))

  private[inference] def wilks(values: Vector[Double]): Either[InferenceError, Double] =
    if values.isEmpty || values.exists(x => !x.isFinite || x < 0.0 || x >= 1.0 || 1.0 - x * x <= 1e-12) then Left(InferenceError.InvalidSpectrum("classical canonical correlations require finite [0,1) values with 1-rho^2 > 1e-12; unit or numerically unit correlations are unavailable"))
    else
      val result = -values.map(x => math.log1p(-x * x)).sum
      if result.isFinite then Right(result) else Left(InferenceError.NonFiniteStatistic("stepwise Wilks", result))

  private def finite(role: String, input: DMat): Either[InferenceError, Unit] =
    MatrixOps.checkFinite(role, input).left.map(error => InferenceError.NumericalFailure(role, error.message))

  private def fullRank(input: DMat): Either[InferenceError, Unit] =
    val rank = input.qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(qrRankTolerance))).diagnostics.rank.getOrElse(0)
    if rank == input.cols then Right(()) else Left(InferenceError.RankLoss(input.cols, rank))

enum CanonicalRankRefitMethod:
  case ExactPermutedTailQrCrossSvd

enum CanonicalRankSampling:
  case WithReplacement
  /** Unrestricted finite-group sampling, without evaluating rejected draws. */
  case DistinctNonIdentity(maximumCandidates: Int)

final class FixedCanonicalRankResult private[inference] (
    val correlations: Vector[Double], val receipts: Vector[MonteCarloReceipt],
    val adjustedPValues: Vector[PValue], val detectableRank: Int,
    val refitMethod: CanonicalRankRefitMethod, val completedCompactFits: Long, val identityDraws: Int,
    val sampling: CanonicalRankSampling, val action: PermutationAction,
    val candidateDraws: Int, val duplicateDraws: Int, val sampledCandidateIds: Vector[ReplicateId],
    val observedCanonicalDiagnostics: SpectralDiagnostics, val qrRankTolerance: Double,
    val maximumNullSvdResidual: Double, val maximumNullOrthogonalityError: Double,
    val allNullFitsExtremeCertified: Boolean
)

object FixedCanonicalRank:
  private val distinctAlgorithm: resample4s.spi.AlgorithmId =
    resample4s.spi.AlgorithmId.of("canonical-rank-distinct-nonidentity/v1").fold(
      error => throw IllegalStateException(error.message),
      (value: resample4s.spi.AlgorithmId) => value
    )
  /** Fixed B, inclusive ties and the existing plus-one reducer. A failed
    * replicate refuses the complete result; it is never treated as planned B.
    * The same legal action is used across every rank hypothesis. The null-value
    * ceiling bounds receipt cardinality, not collection overhead or peak heap. */
  def run(problem: StepwiseCanonicalRank, action: PermutationAction, seed: Seed,
      draws: MonteCarloDraws, alpha: Alpha, maximumNullValues: Long = 1_000_000L,
      sampling: CanonicalRankSampling = CanonicalRankSampling.WithReplacement,
      maximumTransformElements: Long = 1_000_000L): Either[InferenceError, FixedCanonicalRankResult] =
    val required = BigInt(draws.value) * problem.candidateRank
    val distinct = sampling match
      case CanonicalRankSampling.WithReplacement => false
      case CanonicalRankSampling.DistinctNonIdentity(_) => true
    val candidateLimit = sampling match
      case CanonicalRankSampling.WithReplacement => draws.value
      case CanonicalRankSampling.DistinctNonIdentity(limit) => limit
    val transformElements = if distinct then BigInt(draws.value) * problem.leftVariables.rows else BigInt(0)
    val enoughTransforms = if !distinct then true else action match
      case PermutationAction.Unrestricted(rows) =>
        // Capped factorial: answer feasibility without constructing n!.
        val cap = BigInt(draws.value) + 1
        var count = BigInt(1)
        var index = 2
        while index <= rows.value && count < cap do
          count = (count * index).min(cap)
          index += 1
        count >= cap
      case _ => false
    if action.rowCount.value != problem.leftVariables.rows then Left(InferenceError.RowCountMismatch("canonical action", problem.leftVariables.rows, action.rowCount.value))
    else if maximumNullValues < 0L || required > maximumNullValues || required > Int.MaxValue then Left(InferenceError.UnsupportedProblem(s"canonical null receipt requires $required values, allowed $maximumNullValues"))
    else if candidateLimit < draws.value then Left(InferenceError.CanonicalDrawBudgetExhausted(0, 0, draws.value))
    else if !enoughTransforms then Left(InferenceError.UnsupportedProblem("distinct non-identity sampling requires an unrestricted group containing at least B+1 transformations"))
    else if maximumTransformElements < 0L || transformElements > maximumTransformElements || transformElements > Int.MaxValue then Left(InferenceError.UnsupportedProblem(s"distinct transform keys require $transformElements retained row indices, allowed $maximumTransformElements; collection overhead is separate"))
    else
      val values = Array.fill(problem.candidateRank)(Vector.newBuilder[ReplicateStatistic])
      val sampledIds = Vector.newBuilder[ReplicateId]
      val seen = scala.collection.mutable.HashSet.empty[Vector[Int]]
      var index = 0
      var candidates = 0
      var identityDraws = 0
      var duplicates = 0
      var maximumResidual = 0.0
      var maximumOrthogonalityError = 0.0
      var allExtremeCertified = true
      while index < draws.value && candidates < candidateLimit do
        val candidateId = ReplicateId(candidates).toOption.get
        val permutation = action.draw(seed, candidateId) match
          case Left(error) => return Left(InferenceError.ReplicateFailure(candidateId, error))
          case Right(value) => value
        candidates += 1
        val rows = permutation.toIArray
        val isIdentity = (0 until rows.length).forall(i => rows(i) == i)
        if isIdentity then identityDraws += 1
        val accepted = if !distinct then true else if isIdentity then false else
          val novel = seen.add(rows.toVector)
          if !novel then duplicates += 1
          novel
        if accepted then
          val replicate = ReplicateId(index).toOption.get
          val statistics = problem.nullStatisticsWithDiagnostics(permutation) match
            case Left(error) => return Left(InferenceError.ReplicateFailure(candidateId, error))
            case Right(value) => value
          var step = 0
          while step < statistics.size do
            val (statistic, diagnostics) = statistics(step)
            if !diagnostics.orthogonalityError.isFinite || (0 until diagnostics.residuals.length).exists(i => !diagnostics.residuals(i).isFinite) then return Left(InferenceError.ReplicateFailure(candidateId, InferenceError.NumericalFailure("canonical SVD diagnostics", "nonfinite diagnostics")))
            maximumResidual = math.max(maximumResidual, diagnostics.worstResidual)
            maximumOrthogonalityError = math.max(maximumOrthogonalityError, diagnostics.orthogonalityError)
            allExtremeCertified = allExtremeCertified && diagnostics.extremalityCertified
            ReplicateStatistic.from(replicate, statistic) match
              case Left(error) => return Left(InferenceError.ReplicateFailure(candidateId, error))
              case Right(value) => values(step) += value
            step += 1
          sampledIds += candidateId
          index += 1
      if index != draws.value then return Left(InferenceError.CanonicalDrawBudgetExhausted(candidates, index, draws.value))
      val algorithm: resample4s.spi.AlgorithmId = if distinct then distinctAlgorithm else ResamplingPlans.permutationAlgorithm
      val receipts = Vector.newBuilder[MonteCarloReceipt]
      val adjusted = Vector.newBuilder[PValue]
      var maximum = 0.0
      var detectable = 0
      var step = 0
      while step < problem.candidateRank do
        MonteCarlo.fixed(problem.observedWilks(step), Alternative.Greater, values(step).result(),
          ReplicateProvenance.Deterministic(seed, algorithm, Seed.derivationAlgorithm)) match
          case Left(error) => return Left(error)
          case Right(receipt) =>
            receipts += receipt
            maximum = math.max(maximum, receipt.pValue.value)
            adjusted += PValue(maximum).toOption.get
            if detectable == step && maximum <= alpha.value then detectable += 1
        step += 1
      Right(new FixedCanonicalRankResult(problem.correlations, receipts.result(), adjusted.result(), detectable,
        CanonicalRankRefitMethod.ExactPermutedTailQrCrossSvd, draws.value.toLong * problem.candidateRank, identityDraws,
        sampling, action, candidates, duplicates, sampledIds.result(), problem.observedCanonicalDiagnostics, problem.qrRankTolerance,
        maximumResidual, maximumOrthogonalityError, allExtremeCertified))

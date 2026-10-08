package multivar.inference

import gale.linalg.DMat
import resample4s.core.Rand
import resample4s.kernel.{Seed, StreamDomain, StreamPath}
import resample4s.spi.AlgorithmId

/** Complete-null Gaussian comparison for one rank null. Only the smaller
  * block loses variables: dimensions (a-s,b), with m zero-mean rows.
  * The constructor cannot represent the unjustified (a-s,b-s) reference. */
final class GaussianCanonicalReference private[inference] (
    val rows: Int, val smallerColumns: Int, val largerColumns: Int, val nullRank: Int
):
  /** The seed is local to this hypothesis. Both blocks and every replicate
    * use distinct children; there is no fitted covariance or recentering. */
  def draw(hypothesisSeed: Seed, replicate: ReplicateId, maximumElements: Long = 4_000_000L): Either[InferenceError, CanonicalRankSpectrum] =
    val n = BigInt(rows); val p = BigInt(smallerColumns); val q = BigInt(largerColumns)
    val required = 14 * n * (p + q) + 40 * (p * p + q * q + p * q)
    if maximumElements < 0L || required > maximumElements || Vector(n * p, n * q).exists(_ > Int.MaxValue) then
      Left(InferenceError.UnsupportedProblem(s"Gaussian canonical reference requires $required elements, allowed $maximumElements"))
    else
      val seed = GaussianCanonicalRank.child(hypothesisSeed, GaussianCanonicalRank.replicateDomain, replicate.value)
      val left = GaussianCanonicalRank.gaussian(rows, smallerColumns, GaussianCanonicalRank.child(seed, GaussianCanonicalRank.blockDomain, 0))
      val right = GaussianCanonicalRank.gaussian(rows, largerColumns, GaussianCanonicalRank.child(seed, GaussianCanonicalRank.blockDomain, 1))
      CanonicalRankSpectrum.from(left, right, maximumElements)

object GaussianCanonicalReference:
  def forHypothesis(problem: CanonicalRankSpectrum, step: ComponentIx): Either[InferenceError, GaussianCanonicalReference] =
    if step.value >= problem.candidateRank then Left(InferenceError.ComponentOutOfRange(step.value, problem.candidateRank))
    else Right(new GaussianCanonicalReference(problem.rows, problem.candidateRank - step.value,
      math.max(problem.leftColumns, problem.rightColumns), step.value))

object GaussianCanonicalRank:
  val referenceAlgorithm: AlgorithmId = AlgorithmId.of("gaussian-canonical-interlacing-splitmix64-boxmuller52/v1").toOption.get
  private val hypothesisDomain = StreamDomain.custom(1011).toOption.get
  private[inference] val replicateDomain = StreamDomain.custom(1012).toOption.get
  private[inference] val blockDomain = StreamDomain.custom(1013).toOption.get

  def hypothesisSeed(root: Seed, step: ComponentIx): Seed = child(root, hypothesisDomain, step.value)

  /** Under rank(SigmaXY)<=s, rootwise interlacing bounds the unscaled observed
    * Wilks tail by full Wilks for independent blocks (a-s,b). Independent
    * Gaussian Monte Carlo draws and inclusive plus-one counting preserve the
    * stochastic bound. Prefix maxima control false rank overestimation.
    *
    * Assumptions: fixed candidate spaces, iid zero-mean joint Gaussian rows,
    * nonsingular marginal covariances, and rows > p+q. This routine does not
    * infer those assumptions from the input or admit downstream study claims.
    * No failed draw is dropped, retried, or counted as completed B. */
  def run(problem: CanonicalRankSpectrum, seed: Seed, draws: MonteCarloDraws, alpha: Alpha,
      maximumNullValues: Long = 1_000_000L, maximumElements: Long = 4_000_000L
  ): Either[InferenceError, GaussianCanonicalRankResult] =
    val n = BigInt(problem.rows); val p = BigInt(problem.leftColumns); val q = BigInt(problem.rightColumns)
    val nullValues = BigInt(draws.value) * problem.candidateRank
    val required = 14 * n * (p + q) + 40 * (p * p + q * q + p * q) + 4 * nullValues
    if maximumNullValues < 0L || nullValues > maximumNullValues || nullValues > Int.MaxValue then
      Left(InferenceError.UnsupportedProblem(s"Gaussian canonical receipt requires $nullValues values, allowed $maximumNullValues"))
    else if maximumElements < 0L || required > maximumElements then
      Left(InferenceError.UnsupportedProblem(s"Gaussian canonical run requires $required elements, allowed $maximumElements"))
    else
      val receipts = Vector.newBuilder[MonteCarloReceipt]
      val adjusted = Vector.newBuilder[PValue]
      val references = Vector.newBuilder[GaussianCanonicalReference]
      val seeds = Vector.newBuilder[Seed]
      var maximum = 0.0
      var detectable = 0
      var maximumResidual = 0.0
      var maximumOrthogonality = 0.0
      var step = 0
      while step < problem.candidateRank do
        val component = ComponentIx(step).toOption.get
        val reference = GaussianCanonicalReference.forHypothesis(problem, component).toOption.get
        val currentSeed = hypothesisSeed(seed, component)
        val statistics = Vector.newBuilder[ReplicateStatistic]
        var index = 0
        while index < draws.value do
          val replicate = ReplicateId(index).toOption.get
          reference.draw(currentSeed, replicate, maximumElements) match
            case Left(error) => return Left(InferenceError.ReplicateFailure(replicate, error))
            case Right(sample) =>
              val diagnostics = sample.observedCanonicalDiagnostics
              if !diagnostics.orthogonalityError.isFinite || !diagnostics.worstResidual.isFinite then
                return Left(InferenceError.ReplicateFailure(replicate, InferenceError.NumericalFailure("Gaussian reference diagnostics", "nonfinite diagnostics")))
              maximumResidual = math.max(maximumResidual, diagnostics.worstResidual)
              maximumOrthogonality = math.max(maximumOrthogonality, diagnostics.orthogonalityError)
              statistics += ReplicateStatistic.from(replicate, sample.observedWilks.head).toOption.get
          index += 1
        MonteCarlo.fixed(problem.observedWilks(step), Alternative.Greater, statistics.result(),
          ReplicateProvenance.Deterministic(currentSeed, referenceAlgorithm, Seed.derivationAlgorithm)) match
          case Left(error) => return Left(error)
          case Right(receipt) =>
            receipts += receipt
            maximum = math.max(maximum, receipt.pValue.value)
            adjusted += PValue(maximum).toOption.get
            if detectable == step && maximum <= alpha.value then detectable += 1
        references += reference
        seeds += currentSeed
        step += 1
      Right(new GaussianCanonicalRankResult(problem.correlations, receipts.result(), adjusted.result(), detectable,
        references.result(), seeds.result(), draws.value.toLong * problem.candidateRank, required.toLong,
        problem.observedCanonicalDiagnostics, problem.qrRankTolerance, maximumResidual, maximumOrthogonality))

  private[inference] def child(seed: Seed, domain: StreamDomain, ordinal: Int): Seed =
    seed.derive(StreamPath.of(domain, ordinal).toOption.get)

  /** Box-Muller from open, midpoint 52-bit uniforms. Rand owns SplitMix64;
    * only the Gaussian reference transformation is implemented here.
    * DMat.tabulate visits row-major entries; each pair consumes two words. */
  private[inference] def gaussian(rows: Int, columns: Int, seed: Seed): DMat =
    var random = Rand.fromSeed(seed)
    var spare = 0.0
    var hasSpare = false
    def uniform(): Double =
      val (next, word) = random.nextLong
      random = next
      ((word >>> 12).toDouble + 0.5) / 4503599627370496.0
    DMat.tabulate(rows, columns): (_, _) =>
      if hasSpare then
        hasSpare = false
        spare
      else
        val radius = math.sqrt(-2.0 * math.log(uniform()))
        val angle = 2.0 * math.Pi * uniform()
        spare = radius * math.sin(angle)
        hasSpare = true
        radius * math.cos(angle)

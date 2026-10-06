package multivar.inference

import gale.backend.Backend.given
import gale.linalg.{DMat, Matrix}
import resample4s.kernel.Seed
import scala.util.control.NonFatal

/** Sampling is part of the test, not inferred from the number of supplied draws. */
enum MeanProfileSampling:
  case Exact
  case MonteCarlo(draws: Int, seed: Long)
  /** Caller attests that supplied draws were uniform, independent, with replacement. */
  case ReplayMonteCarlo(signs: Vector[Vector[Int]], source: String)

final case class MeanProfileLimits(maxDraws: Int = 100000, maxStoredSigns: Long = 5000000L,
    maxProfileElements: Long = 50000000L, maxGramElements: Long = 1000000L):
  require(maxDraws > 0 && maxStoredSigns > 0 && maxProfileElements > 0 && maxGramElements > 0)

final case class MeanProfileDraw(signs: Vector[Int], statistic: Double)
final case class MeanProfileSignFailure(reason: String, requested: Int,
    completed: Vector[MeanProfileDraw], failed: Int, cancelled: Boolean)

/** Complete successful evidence only. No scientific symmetry claim is inferred. */
final case class MeanProfileSignReceipt private[inference] (
    sampling: String, provenance: String, participants: Int, activeProfiles: Int,
    groupSize: BigInt, dataOrbitSize: BigInt, requested: Int,
    observed: Double, draws: Vector[MeanProfileDraw], exceedances: Int,
    identityOccurrences: Int, pValue: PValue, minimumP: Double,
    relativeTolerance: Double, absoluteTolerance: Double, comparisonTolerance: Double
)

/** A complete-profile mean-energy statistic with a participant-sized Gram cache.
  * Coordinates are never selected using observed signal or numerical rank.
  */
final class PreparedMeanProfileSignTest private[inference] (
    private val gram: DMat, private val scaleSquared: Double,
    val activeProfiles: Int, val referenceEnergy: Double,
    val relativeTolerance: Double, val absoluteTolerance: Double,
    private val limits: MeanProfileLimits
):
  val participants: Int = gram.rows
  val comparisonTolerance: Double = relativeTolerance * referenceEnergy + absoluteTolerance
  private val gramRoundoff = 64 * Math.ulp(1.0) * participants *
    (0 until participants).map(i => gram(i,i)).sum / participants

  def statistic(signs: Vector[Int]): Either[InferenceError, Double] =
    if signs.size != participants || signs.exists(s => s != -1 && s != 1) then
      Left(InferenceError.UnsupportedProblem("One -1/+1 sign is required per complete participant profile"))
    else
      var sum = 0.0
      var i = 0
      while i < participants do
        sum += gram(i,i)
        var j = 0
        while j < i do
          sum += 2 * signs(i) * signs(j) * gram(i,j)
          j += 1
        i += 1
      val value = sum / (participants.toDouble * participants)
      if !value.isFinite || value < -gramRoundoff then
        Left(InferenceError.NumericalFailure("mean-profile sign test","Mean-profile Gram energy is non-finite or materially negative"))
      else
        val physical = math.max(0.0,value) * scaleSquared
        if physical.isFinite then Right(physical)
        else Left(InferenceError.NonFiniteStatistic("mean-profile energy",physical))

  def run(sampling: MeanProfileSampling, cancelled: () => Boolean = () => false): Either[MeanProfileSignFailure, MeanProfileSignReceipt] =
    val groupSize = BigInt(1) << participants
    val count = sampling match
      case MeanProfileSampling.Exact => if groupSize.isValidInt then groupSize.toInt else -1
      case MeanProfileSampling.MonteCarlo(draws,_) => draws
      case MeanProfileSampling.ReplayMonteCarlo(signs,_) => signs.size
    val completed = scala.collection.mutable.ArrayBuffer.empty[MeanProfileDraw]
    def fail(reason: String, failed: Int = 0, stopped: Boolean = false) =
      Left(MeanProfileSignFailure(reason,count,completed.toVector,failed,stopped))
    if count <= 0 || count > limits.maxDraws || count.toLong * participants > limits.maxStoredSigns then
      return fail("Sign distribution exceeds the declared draw or stored-sign budget")
    val (mode,provenance) = sampling match
      case MeanProfileSampling.Exact => ("exact","All independent signs; identity first")
      case MeanProfileSampling.MonteCarlo(_,seed) => ("monte-carlo",s"SignFlipAction/independent/domain-1003/seed=$seed")
      case MeanProfileSampling.ReplayMonteCarlo(_,source) =>
        if source.trim.isEmpty then return fail("Replay requires explicit sampling provenance")
        ("monte-carlo",s"External uniform independent replacement draws: $source")
    try
      if cancelled() then return fail("Cancelled before observation",stopped = true)
      val observed = statistic(Vector.fill(participants)(1)) match
        case Left(error) => return fail(error.toString,failed = 1)
        case Right(value) => value
      val action = RowCount(participants).flatMap(SignFlipAction.independent) match
        case Left(error) => return fail(error.toString,failed = 1)
        case Right(value) => value
      var index = 0
      while index < count do
        if cancelled() then return fail("Cancelled during sign distribution",stopped = true)
        val signs = sampling match
          case MeanProfileSampling.Exact => Right(Vector.tabulate(participants)(p => if ((index >>> p) & 1) == 0 then 1 else -1))
          case MeanProfileSampling.MonteCarlo(_,seed) =>
            ReplicateId(index).flatMap(action.draw(Seed.fromLong(seed),_)).map(_.signs)
          case MeanProfileSampling.ReplayMonteCarlo(supplied,_) => Right(supplied(index))
        signs.flatMap(s => statistic(s).map(t => MeanProfileDraw(s,t))) match
          case Left(error) => return fail(error.toString,failed = 1)
          case Right(draw) => completed += draw
        index += 1
      val draws = completed.toVector
      val exceedances = draws.count(d => d.statistic >= observed || observed - d.statistic <= comparisonTolerance)
      val p = if mode == "exact" then exceedances.toDouble / count else (1.0 + exceedances) / (count + 1.0)
      PValue(p) match
        case Left(error) => fail(error.toString,failed = 1)
        case Right(value) => Right(MeanProfileSignReceipt(mode,provenance,participants,activeProfiles,
          groupSize,BigInt(1) << activeProfiles,count,observed,draws,exceedances,
          draws.count(_.signs.forall(_ == 1)),value,
          if mode != "exact" then 1.0 / (count + 1.0) else if activeProfiles == 0 then 1.0 else math.pow(2.0,1-activeProfiles),
          relativeTolerance,absoluteTolerance,comparisonTolerance))
    catch case NonFatal(error) => fail(Option(error.getMessage).getOrElse(error.toString),failed = 1)

object MeanProfileSignTest:
  /** Input rows are independent units; columns jointly constitute each complete profile.
    * The optional absolute tolerance is in squared response units and fixed for the orbit.
    */
  def prepare(profiles: DMat, relativeTolerance: Double = 1e-12, absoluteTolerance: Double = 0.0,
      limits: MeanProfileLimits = MeanProfileLimits()): Either[InferenceError,PreparedMeanProfileSignTest] =
    if profiles.rows < 2 || profiles.cols < 1 || profiles.rows.toLong * profiles.cols > limits.maxProfileElements ||
        profiles.rows.toLong * profiles.rows > limits.maxGramElements then
      return Left(InferenceError.UnsupportedProblem("At least two profiles and one coordinate required within profile/Gram budgets"))
    if !relativeTolerance.isFinite || relativeTolerance < 0 || relativeTolerance > 1e-6 ||
        !absoluteTolerance.isFinite || absoluteTolerance < 0 then
      return Left(InferenceError.UnsupportedProblem("Invalid mean-profile comparison tolerances"))
    var scale = 0.0
    var active = 0
    var r = 0
    while r < profiles.rows do
      var nonzero = false
      var c = 0
      while c < profiles.cols do
        val value = profiles(r,c)
        if !value.isFinite then return Left(InferenceError.NonFiniteStatistic("mean profile",value))
        scale = math.max(scale,math.abs(value))
        nonzero ||= value != 0.0
        c += 1
      if nonzero then active += 1
      r += 1
    val scaleSquared = if scale == 0 then 1.0 else scale * scale
    if !scaleSquared.isFinite || scaleSquared == 0 then
      return Left(InferenceError.NumericalFailure("mean-profile sign test","Squared profile scale is outside finite floating-point range"))
    try
      val builder = Matrix.newBuilder(profiles.rows,profiles.cols)
      for r <- 0 until profiles.rows; c <- 0 until profiles.cols do
        builder(r,c) = if scale == 0 then 0.0 else profiles(r,c) / scale
      val normalized = builder.result()
      val gram = normalized * normalized.t
      val energy = (0 until profiles.rows).map(i => gram(i,i)).sum / profiles.rows * scaleSquared
      if !energy.isFinite || !(relativeTolerance * energy + absoluteTolerance).isFinite then
        Left(InferenceError.NonFiniteStatistic("mean profile reference energy",energy))
      else Right(new PreparedMeanProfileSignTest(gram,scaleSquared,active,energy,relativeTolerance,absoluteTolerance,limits))
    catch case NonFatal(error) => Left(InferenceError.NumericalFailure("mean-profile sign test",Option(error.getMessage).getOrElse(error.toString)))

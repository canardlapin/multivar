package multivar.inference

/** A scalar estimate and a strictly positive, finite standard error, not a
  * variance. Construction does not establish scientific interval admission.
  */
final class StudentizedValue private (val estimate: Double, val standardError: Double)

object StudentizedValue:
  def apply(estimate: Double, standardError: Double): Either[InferenceError, StudentizedValue] =
    if !estimate.isFinite then Left(InferenceError.NonFiniteStatistic("studentized estimate", estimate))
    else if !standardError.isFinite || standardError <= 0 then
      Left(InferenceError.NumericalFailure("studentized standard error", "requires a finite positive SE, not variance"))
    else Right(new StudentizedValue(estimate, standardError))

/** Every requested deletion appears once, including identical copies. */
object JackknifeStandardError:
  def compute(deletions: IndexedSeq[Double], expectedDeletions: Int,
      cancelled: () => Boolean = () => false): Either[InferenceError, Double] =
    if expectedDeletions < 2 || deletions.size != expectedDeletions then
      return Left(InferenceError.InvalidReplicatePlan("jackknife requires exactly N deletions for N >= 2 observations"))
    if cancelled() then return ProfileBootstrapMath.cancellation
    val origin = deletions.head
    var absoluteScale = 0.0
    var differenceScale = 0.0
    var finiteDifferences = true
    var i = 0
    while i < deletions.size do
      if cancelled() then return ProfileBootstrapMath.cancellation
      val value = deletions(i)
      if !value.isFinite then return Left(InferenceError.NonFiniteStatistic("jackknife deletion", value))
      absoluteScale = Math.max(absoluteScale, Math.abs(value))
      val difference = value - origin
      finiteDifferences = finiteDifferences && difference.isFinite
      differenceScale = Math.max(differenceScale, Math.abs(difference))
      i += 1
    // Anchor first to preserve tiny variation on a large location. Scale the
    // raw values only when subtraction itself would overflow.
    val scale = if finiteDifferences then differenceScale else absoluteScale
    if scale == 0.0 then return Right(0.0)
    def normalized(i: Int): Double =
      if finiteDifferences then (deletions(i) - origin) / scale
      else deletions(i) / scale - origin / scale
    var sum = 0.0
    var correction = 0.0
    i = 0
    while i < deletions.size do
      if cancelled() then return ProfileBootstrapMath.cancellation
      val term = normalized(i)
      val next = sum + term
      correction += (if Math.abs(sum) >= Math.abs(term) then (sum - next) + term else (term - next) + sum)
      sum = next
      i += 1
    val mean = (sum + correction) / expectedDeletions
    var norm = 0.0
    i = 0
    while i < deletions.size do
      if cancelled() then return ProfileBootstrapMath.cancellation
      norm = Math.hypot(norm, normalized(i) - mean)
      i += 1
    val se = scale * (norm * Math.sqrt((expectedDeletions - 1.0) / expectedDeletions))
    if !se.isFinite || (norm > 0 && se == 0) then
      Left(InferenceError.NumericalFailure("jackknife SE", "standard error exceeds floating-point range"))
    else Right(se)

final case class StudentizedDraw(id: ReplicateId, value: Either[InferenceError, StudentizedValue])

final case class StudentizedIntervalOptions(
    confidenceLevel: Double,
    minimumStandardError: Double = 0.0,
    maxWorkingBytes: Long = 64L * 1024 * 1024
)

enum IntervalQuantileConvention:
  case Type7

/** Complete-draw pointwise bootstrap-t arithmetic; no scientific admission is
  * inferred from successful numerical construction. No transformation or extra
  * bias correction is applied. Caller owns the meaning/alignment of each scalar.
  */
final class StudentizedInterval private (
    val observed: StudentizedValue,
    val replicateIds: Vector[ReplicateId],
    val confidenceLevel: Double,
    val lowerPivotQuantile: Double,
    val upperPivotQuantile: Double,
    val lower: Double,
    val upper: Double
):
  val quantileConvention: IntervalQuantileConvention = IntervalQuantileConvention.Type7

object StudentizedInterval:
  def compute(observed: StudentizedValue, expectedIds: IndexedSeq[ReplicateId],
      draws: Iterator[StudentizedDraw], options: StudentizedIntervalOptions,
      cancelled: () => Boolean = () => false): Either[InferenceError, StudentizedInterval] =
    if !options.confidenceLevel.isFinite || options.confidenceLevel <= 0 || options.confidenceLevel >= 1 then
      return Left(InferenceError.InvalidProbability("interval confidence", options.confidenceLevel, false))
    if !options.minimumStandardError.isFinite || options.minimumStandardError < 0 then
      return Left(InferenceError.InvalidTolerance("minimum interval SE", options.minimumStandardError))
    if expectedIds.size < 2 then return Left(InferenceError.InvalidCount("studentized draws", expectedIds.size))
    if cancelled() then return ProfileBootstrapMath.cancellation
    // Conservative bound includes requested IDs, uniqueness validation, retained
    // pivots, sorting and returned IDs. The caller budgets its fit/readout state.
    ProfileBootstrapMath.allocation(BigInt(128) * expectedIds.size + 1024, options.maxWorkingBytes) match
      case Left(error) => return Left(error)
      case Right(_) => ()
    if expectedIds.distinct.size != expectedIds.size then
      return Left(InferenceError.InvalidReplicatePlan("studentized replicate identities must be unique"))
    if observed.standardError <= options.minimumStandardError then
      return Left(InferenceError.NumericalFailure("observed interval SE", "standard error does not exceed the declared floor"))
    val lowProbability = (1 - options.confidenceLevel) / 2
    val highProbability = 1 - lowProbability
    if lowProbability >= highProbability then
      return Left(InferenceError.NumericalFailure("interval tails", "tail probabilities are not numerically distinct"))
    val pivots = new Array[Double](expectedIds.size)
    var index = 0
    while index < expectedIds.size do
      if cancelled() then return ProfileBootstrapMath.cancellation
      if !draws.hasNext then return Left(InferenceError.InvalidReplicatePlan("studentized stream ended before every requested draw"))
      val draw = draws.next()
      if draw.id != expectedIds(index) then
        return Left(InferenceError.InvalidReplicatePlan("studentized draw identity/order differs from the frozen plan"))
      draw.value match
        case Left(error) => return Left(InferenceError.ReplicateFailure(draw.id, error))
        case Right(value) =>
          if value.standardError <= options.minimumStandardError then
            return Left(InferenceError.ReplicateFailure(draw.id,
              InferenceError.NumericalFailure("draw interval SE", "standard error does not exceed the declared floor")))
          val difference = value.estimate - observed.estimate
          val pivot = if difference.isFinite then difference / value.standardError
            else value.estimate / value.standardError - observed.estimate / value.standardError
          if !pivot.isFinite then return Left(InferenceError.ReplicateFailure(draw.id,
            InferenceError.NonFiniteStatistic("studentized pivot", pivot)))
          pivots(index) = pivot
      index += 1
    if draws.hasNext then return Left(InferenceError.InvalidReplicatePlan("studentized stream contains unrequested draws"))
    if cancelled() then return ProfileBootstrapMath.cancellation
    scala.util.Sorting.quickSort(pivots)
    def quantile(probability: Double): Double =
      val position = (pivots.length - 1.0) * probability
      val from = Math.floor(position).toInt
      val fraction = position - from
      if from == pivots.length - 1 || fraction == 0 then pivots(from)
      else
        val a = pivots(from)
        val b = pivots(from + 1)
        val difference = b - a
        if difference.isFinite then a + fraction * difference
        else (1 - fraction) * a + fraction * b
    val low = quantile(lowProbability)
    val high = quantile(highProbability)
    val lower = observed.estimate - high * observed.standardError
    val upper = observed.estimate - low * observed.standardError
    if cancelled() then ProfileBootstrapMath.cancellation
    else if !lower.isFinite || !upper.isFinite || lower > upper then
      Left(InferenceError.NumericalFailure("studentized interval", "endpoints exceed finite range or are reversed"))
    else Right(new StudentizedInterval(observed, expectedIds.toVector, options.confidenceLevel, low, high, lower, upper))

package multivar.inference

import gale.linalg.{DMat, LinAlgError, Matrix}
import gale.spectral.{MatrixEnclosure, RealInterval, SingularSelection, Svds, ValidatedSvd}
import resample4s.core.Rand
import resample4s.kernel.Seed

/** Fixed coordinate provenance, supplied by the scientific adapter. */
final case class TaskRankCoordinates(version: String, transform: String, source: String, features: Vector[String])
final case class TaskRankParticipant(id: String, profile: DMat)

/** A fixed row map and positive weight, including their formation uncertainty.
  * All participants in a group use the same map. A matrix interval encloses one
  * fixed scientific map; it is not permission to learn a map from these data.
  */
final case class TaskRankGroup(
    id: String,
    participants: Vector[TaskRankParticipant],
    rowMap: MatrixEnclosure,
    weight: RealInterval
)
final case class TaskRankSigns(id: String, groups: Vector[Vector[Int]])
enum TaskRankMethod:
  case FrobeniusFeasibility, OperatorFirstThenBound
  def version: String = this match
    case FrobeniusFeasibility => "task-mean-rank/frobenius-feasibility/v1"
    case OperatorFirstThenBound => "task-mean-rank/operator-first-bound/v1"
enum TaskRankSampling:
  case ExhaustiveFullSignGroup
  case FixedUniform(seed: Long, algorithm: String)
  /** Replay is conditional on the caller's fixed-budget iid uniform full-group
    * sampling declaration. Supplied arrays cannot prove that sampling law. */
  case SuppliedFixedUniformDeclaration

/** Private construction prevents an incomplete set being labeled exhaustive. */
final class TaskRankPlan private (
    val groupSizes: Vector[Int], val draws: Vector[TaskRankSigns], val sampling: TaskRankSampling
)
object TaskRankPlan:
  private def invalid(message: String) = Left(InferenceError.InvalidReplicatePlan(message))
  private def validSizes(sizes: Vector[Int]): Boolean =
    sizes.nonEmpty && sizes.forall(_ > 0) && sizes.map(_.toLong).sum <= Int.MaxValue

  /** Recheck all representation invariants at execution time. Scala-private
    * constructors remain callable through the JVM ABI, so construction alone
    * cannot be the only boundary for a plan that affects a denominator. */
  private[inference] def validate(plan: TaskRankPlan): Either[InferenceError, Unit] =
    val draws = plan.draws
    if !validSizes(plan.groupSizes) || draws.isEmpty then
      invalid("task-rank requires positive group sizes and a nonempty fixed budget")
    else if draws.map(_.id).distinct.size != draws.size || draws.exists(_.id.trim.isEmpty) then
      invalid("task-rank draw ids must be nonempty and unique")
    else if draws.exists(draw =>
        draw.groups.map(_.size) != plan.groupSizes ||
          draw.groups.exists(_.exists(sign => sign != -1 && sign != 1))) then
      invalid("one exact -1/+1 sign is required per identified participant")
    else
      plan.sampling match
        case TaskRankSampling.ExhaustiveFullSignGroup =>
          val dimensions = plan.groupSizes.map(_.toLong).sum
          if dimensions >= 31 || draws.size != (1 << dimensions.toInt) then
            invalid("exhaustive task-rank plans must contain every full-group sign pattern")
          else if draws.map(_.groups).distinct.size != draws.size then
            invalid("exhaustive task-rank plans must contain unique full-group sign patterns")
          else Right(())
        case _ => Right(())

  def replayFixed(sizes: Vector[Int], draws: Vector[TaskRankSigns]): Either[InferenceError, TaskRankPlan] =
    if !validSizes(sizes) || draws.isEmpty then invalid("task-rank requires positive group sizes and a nonempty fixed budget")
    else if draws.map(_.id).distinct.size != draws.size || draws.exists(_.id.trim.isEmpty) then invalid("task-rank draw ids must be nonempty and unique")
    else if draws.exists(d => d.groups.map(_.size) != sizes || d.groups.exists(_.exists(s => s != -1 && s != 1))) then
      invalid("one exact -1/+1 sign is required per identified participant")
    else Right(new TaskRankPlan(sizes, draws, TaskRankSampling.SuppliedFixedUniformDeclaration))

  /** Uses the provider's integer SplitMix64 bounded draws. Choose seed and B
    * independently of the response data; do not retry until significance. */
  def fixedUniform(sizes: Vector[Int], budget: Int, seed: Seed): Either[InferenceError, TaskRankPlan] =
    if !validSizes(sizes) || budget <= 0 then invalid("task-rank requires positive group sizes and fixed budget")
    else
      var rng = Rand.fromSeed(seed)
      val draws = Vector.newBuilder[TaskRankSigns]
      var b = 0
      while b < budget do
        val groups = Vector.newBuilder[Vector[Int]]
        var j = 0
        while j < sizes.size do
          val signs = Vector.newBuilder[Int]
          var i = 0
          while i < sizes(j) do
            rng.nextIntBounded(2) match
              case Left(error) => return Left(InferenceError.Resampling(error))
              case Right((next, value)) =>
                rng = next
                signs += (if value == 0 then -1 else 1)
            i += 1
          groups += signs.result()
          j += 1
        draws += TaskRankSigns(b.toString, groups.result())
        b += 1
      Right(new TaskRankPlan(sizes, draws.result(), TaskRankSampling.FixedUniform(seed.value, "resample4s/SplitMix64-bounded2/v1")))

  /** Materialized enumeration is deliberately budget bounded. No sign pattern
    * is discarded, including identity and constant-within-group patterns. */
  def exhaustive(sizes: Vector[Int], maximumDraws: Int): Either[InferenceError, TaskRankPlan] =
    if !validSizes(sizes) || maximumDraws <= 0 then invalid("invalid exhaustive task-rank budget")
    else
      val n = sizes.sum
      if n >= 31 || (1L << n) > maximumDraws.toLong then invalid("full sign group exceeds declared materialization budget")
      else
        val offsets = sizes.scanLeft(0)(_ + _)
        val draws = Vector.tabulate(1 << n) { pattern =>
          TaskRankSigns(pattern.toString, sizes.indices.map { j =>
            Vector.tabulate(sizes(j))(i => if ((pattern >>> (offsets(j) + i)) & 1) == 0 then -1 else 1)
          }.toVector)
        }
        Right(new TaskRankPlan(sizes, draws, TaskRankSampling.ExhaustiveFullSignGroup))

enum TaskRankDecision:
  case Counted, Excluded, Ambiguous, ConservativeCarry
final case class TaskRankComparison(rank: Int, decision: TaskRankDecision,
    statisticLower: Double, statisticUpper: Double, thresholdLower: Double, thresholdUpper: Double)
final case class TaskRankDrawReceipt(id: String, signs: Vector[Vector[Int]], comparisons: Vector[TaskRankComparison])
final case class TaskRankGroupReceipt(id: String, participantIds: Vector[String], rowMap: MatrixEnclosure, weight: RealInterval)
final case class TaskRankReceipt(
    method: TaskRankMethod, coordinates: TaskRankCoordinates, sampling: TaskRankSampling,
    groups: Vector[TaskRankGroupReceipt], draws: Vector[TaskRankDrawReceipt],
    counts: Vector[BigInt], denominator: BigInt, alphaNumerator: BigInt, alphaDenominator: BigInt,
    testedRanks: Vector[Int], lowerRankBound: Int
):
  /** Display only. Rejection/stopping uses the exact integer comparison. */
  def pValues: Vector[Double] = counts.map(_.toDouble / denominator.toDouble)
  val transformationGroup: String = "independent whole-participant full sign group"
  val estimand: String = "rank of stacked positive-weighted population group profile means after fixed row maps"

object TaskMeanRank:
  private type Calc[A] = Either[InferenceError, A]
  private type Entries = Vector[RealInterval]
  private final case class Prepared(group: TaskRankGroup, profiles: Vector[Entries], observed: Entries, rows: Int, columns: Int)
  private final case class Transformed(signed: Entries, centered: Entries, z: Entries, rho: RealInterval,
      oneMinusAbsA: RealInterval, constant: Boolean)
  private val zero = RealInterval.exact(0.0).fold(e => throw IllegalStateException(e.toString), identity)
  private def invalid(message: String): Calc[Nothing] = Left(InferenceError.InvalidReplicatePlan(message))
  private def arithmetic[A](value: Either[LinAlgError, A]): Calc[A] =
    value.left.map(e => InferenceError.NumericalFailure("certified task-rank arithmetic", e.toString))
  private def exact(value: Double): Calc[RealInterval] = arithmetic(RealInterval.exact(value))
  private def sequence[A, B](values: Vector[A])(f: A => Calc[B]): Calc[Vector[B]] =
    val out = Vector.newBuilder[B]
    var i = 0
    while i < values.size do
      f(values(i)) match
        case Left(error) => return Left(error)
        case Right(value) => out += value
      i += 1
    Right(out.result())
  private def sum(values: Entries): Calc[RealInterval] =
    values.foldLeft[Calc[RealInterval]](Right(zero))((previous, value) => previous.flatMap(x => arithmetic(x.add(value))))
  private def nonnegative(value: RealInterval): Calc[RealInterval] =
    arithmetic(RealInterval.checked(math.max(0.0, value.lower), value.upper))
  private def squaredNorm(values: Entries): Calc[RealInterval] =
    sequence(values)(x => arithmetic(x.square)).flatMap(sum).flatMap(nonnegative)

  /** Every draw is evaluated before any receipt is returned. Arithmetic or SVD
    * certificate failure and cancellation return Left, never partial evidence.
    * Finite-sample validity remains conditional on the input symmetry model and
    * the declared fixed sampling law; this API cannot diagnose population laws.
    */
  def run(groups: Vector[TaskRankGroup], coordinates: TaskRankCoordinates, plan: TaskRankPlan,
      method: TaskRankMethod, alphaNumerator: BigInt, alphaDenominator: BigInt,
      isCancelled: () => Boolean = () => false): Either[InferenceError, TaskRankReceipt] =
    if alphaNumerator <= 0 || alphaDenominator <= alphaNumerator then
      return invalid("alpha must be an exact rational strictly between zero and one")
    TaskRankPlan.validate(plan) match
      case Left(error) => return Left(error)
      case Right(_)    => ()
    if groups.isEmpty || groups.map(_.participants.size) != plan.groupSizes then invalid("task-rank plan does not match input groups")
    else if Vector(coordinates.version, coordinates.transform, coordinates.source).exists(_.trim.isEmpty) then invalid("coordinate and source provenance must be nonempty")
    else if coordinates.features.isEmpty || coordinates.features.exists(_.trim.isEmpty) || coordinates.features.distinct.size != coordinates.features.size then invalid("feature coordinates must be nonempty, ordered and unique")
    else if groups.exists(_.id.trim.isEmpty) || groups.map(_.id).distinct.size != groups.size then invalid("group ids must be nonempty and unique")
    else if groups.flatMap(_.participants.map(_.id)).distinct.size != groups.map(_.participants.size).sum || groups.exists(_.participants.exists(_.id.trim.isEmpty)) then
      invalid("whole participant identities must be nonempty and globally unique")
    else if isCancelled() then invalid("task-rank execution cancelled before preparation")
    else
      for
        prepared <- sequence(groups)(g => prepare(g, coordinates.features.size))
        observed <- spectrum(prepared.flatMap(_.observed), prepared.map(_.rows).sum, coordinates.features.size)
        rows <- sequence(plan.draws) { draw =>
          if isCancelled() then invalid("task-rank execution cancelled; no partial receipt")
          else evaluate(prepared, draw, method, observed)
        }
        receipt <-
          if isCancelled() then invalid("task-rank execution cancelled; no partial receipt")
          else
            val plusOne = if plan.sampling == TaskRankSampling.ExhaustiveFullSignGroup then BigInt(0) else BigInt(1)
            val denominator = BigInt(rows.size) + plusOne
            val counts = Vector.tabulate(observed.size + 1) { r =>
              plusOne + BigInt(rows.count(_.comparisons(r).decision != TaskRankDecision.Excluded))
            }
            // The terminal rank-q null always counts and alpha<1, so this exists.
            val stop = counts.indices.find(r => counts(r) * alphaDenominator > alphaNumerator * denominator).getOrElse(observed.size)
            Right(TaskRankReceipt(method, coordinates, plan.sampling,
              groups.map(g => TaskRankGroupReceipt(g.id, g.participants.map(_.id), g.rowMap, g.weight)),
              rows, counts, denominator, alphaNumerator, alphaDenominator, (0 to stop).toVector, stop))
      yield receipt

  private def prepare(group: TaskRankGroup, features: Int): Calc[Prepared] =
    val map = group.rowMap
    if group.weight.lower <= 0.0 || map.lower.rows <= 0 || map.lower.cols <= 0 ||
        group.participants.exists(p => p.profile.rows != map.lower.cols || p.profile.cols != features) then
      invalid("positive weights and complete shape-compatible participant profiles are required")
    else
      for
        profiles <- sequence(group.participants) { participant =>
          sequence((0 until map.lower.rows * features).toVector) { index =>
            val row = index / features
            val column = index % features
            sequence((0 until map.lower.cols).toVector) { k =>
              for
                coefficient <- arithmetic(RealInterval.checked(map.lower(row, k), map.upper(row, k)))
                value <- exact(participant.profile(k, column))
                product <- arithmetic(coefficient.multiply(value))
              yield product
            }.flatMap(sum)
          }
        }
        n <- exact(profiles.size.toDouble)
        scale <- arithmetic(group.weight.divide(n))
        observed <- sequence(profiles.head.indices.toVector) { index =>
          sum(profiles.map(_(index))).flatMap(x => arithmetic(x.multiply(scale)))
        }
      yield Prepared(group, profiles, observed, map.lower.rows, features)

  private def transform(group: Prepared, signs: Vector[Int]): Calc[Transformed] =
    val positive = signs.count(_ == 1)
    val negative = signs.size - positive
    val constant = positive == 0 || negative == 0
    for
      n <- exact(signs.size.toDouble)
      scale <- arithmetic(group.group.weight.divide(n))
      signed <- sequence(group.observed.indices.toVector) { index =>
        sequence(signs.indices.toVector)(i => exact(signs(i).toDouble).flatMap(s => arithmetic(group.profiles(i)(index).multiply(s))))
          .flatMap(sum).flatMap(x => arithmetic(x.multiply(scale)))
      }
      result <- if constant then
        Right(Transformed(signed, Vector.fill(signed.size)(zero), Vector.fill(signed.size)(zero), zero, zero, true))
      else
        for
          difference <- exact((positive - negative).toDouble)
          a <- arithmetic(difference.divide(n))
          p <- exact(positive.toDouble).flatMap(x => arithmetic(x.divide(n)))
          m <- exact(negative.toDouble).flatMap(x => arithmetic(x.divide(n)))
          four <- exact(4.0)
          h <- arithmetic(p.multiply(m)).flatMap(x => arithmetic(x.multiply(four)))
          rootH <- arithmetic(h.sqrt)
          aOverRootH <- arithmetic(a.divide(rootH))
          centered <- sequence(signed.indices.toVector) { i =>
            arithmetic(a.multiply(group.observed(i))).flatMap(x => arithmetic(signed(i).subtract(x)))
          }
          z <- sequence(signed.indices.toVector) { i =>
            for
              left <- arithmetic(rootH.multiply(group.observed(i)))
              right <- arithmetic(aOverRootH.multiply(centered(i)))
              value <- arithmetic(left.subtract(right))
            yield value
          }
          rho <- squaredNorm(centered).flatMap(x => arithmetic(x.divide(h))).flatMap(nonnegative)
          twiceMinority <- exact(2.0 * math.min(positive, negative).toDouble)
          oneMinus <- arithmetic(twiceMinority.divide(n))
        yield Transformed(signed, centered, z, rho, oneMinus, false)
    yield result

  private def spectrum(entries: Entries, rows: Int, columns: Int): Calc[Entries] =
    if entries.forall(x => x.lower == 0.0 && x.upper == 0.0) then
      return Right(Vector.fill(math.min(rows, columns))(zero))
    val lower = Matrix.tabulate(rows, columns)((r, c) => entries(r * columns + c).lower)
    val upper = Matrix.tabulate(rows, columns)((r, c) => entries(r * columns + c).upper)
    // Midpoints are solver candidates only; certification includes input formation.
    val middle = Matrix.tabulate(rows, columns)((r, c) => .5 * lower(r, c) + .5 * upper(r, c))
    for
      enclosure <- arithmetic(MatrixEnclosure.checked(lower, upper))
      candidate <- arithmetic(Svds.svd(middle, SingularSelection.All))
      certified <- arithmetic(ValidatedSvd.enclosure(enclosure, candidate))
      result <- sequence((0 until certified.lower.length).toVector)(i => arithmetic(RealInterval.checked(certified.lower(i), certified.upper(i))))
    yield result

  private def compare(rank: Int, statistic: RealInterval, threshold: RealInterval, lowerTail: Boolean): TaskRankComparison =
    val decision =
      if lowerTail then
        if statistic.upper <= threshold.lower then TaskRankDecision.Counted
        else if statistic.lower > threshold.upper then TaskRankDecision.Excluded
        else TaskRankDecision.Ambiguous
      else
        if statistic.lower >= threshold.upper then TaskRankDecision.Counted
        else if statistic.upper < threshold.lower then TaskRankDecision.Excluded
        else TaskRankDecision.Ambiguous
    TaskRankComparison(rank, decision, statistic.lower, statistic.upper, threshold.lower, threshold.upper)

  private def evaluate(groups: Vector[Prepared], draw: TaskRankSigns, method: TaskRankMethod, observed: Entries): Calc[TaskRankDrawReceipt] =
    for
      transformed <- sequence(groups.indices.toVector)(j => transform(groups(j), draw.groups(j)))
      comparisons <- method match
        case TaskRankMethod.FrobeniusFeasibility =>
          for
            singular <- spectrum(transformed.flatMap(_.z), groups.map(_.rows).sum, groups.head.columns)
            rho <- sum(transformed.map(_.rho)).flatMap(nonnegative)
            comparisons <- sequence((0 to observed.size).toVector) { r =>
              squaredNorm(singular.drop(r)).map(tail => compare(r, tail, rho, lowerTail = true))
            }
          yield comparisons
        case TaskRankMethod.OperatorFirstThenBound =>
          for
            signed <- spectrum(transformed.flatMap(_.signed), groups.map(_.rows).sum, groups.head.columns)
            centered <- spectrum(transformed.flatMap(_.centered), groups.map(_.rows).sum, groups.head.columns)
            minimum <- arithmetic(RealInterval.checked(transformed.map(_.oneMinusAbsA.lower).min, transformed.map(_.oneMinusAbsA.upper).min))
            later <- sequence((1 until observed.size).toVector) { r =>
              if transformed.exists(_.constant) then Right(compare(r, zero, zero, lowerTail = false))
              else arithmetic(minimum.multiply(observed(r))).map(threshold => compare(r, centered.head, threshold, lowerTail = false))
            }
          yield compare(0, signed.head, observed.head, lowerTail = false) +: later :+ compare(observed.size, zero, zero, lowerTail = false)
    yield
      // Enclose the exact monotone indicators even when separately computed
      // numerical bounds made an earlier rank ambiguous. Record each carry.
      var previouslyCounted = false
      val monotone = comparisons.map { c =>
        val result = if previouslyCounted && c.decision == TaskRankDecision.Excluded then c.copy(decision = TaskRankDecision.ConservativeCarry) else c
        previouslyCounted = result.decision != TaskRankDecision.Excluded
        result
      }
      TaskRankDrawReceipt(draw.id, draw.groups, monotone)

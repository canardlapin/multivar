package multivar.inference

import gale.linalg.{DMat, DVec}
import gale.spectral.{MatrixEnclosure, SVD, ValidatedSvd}

/** Versioned, formed-input boundary for certified task-mean rank comparisons.
  * The caller must form Z and rho with directed interval arithmetic and record
  * its fixed coordinate map; this generic API refuses unbounded comparisons.
  */
final case class TaskRankCoordinates(version: String, transform: String, source: String):
  require(version.nonEmpty && transform.nonEmpty && source.nonEmpty)
final case class TaskRankDraw(id: String, signs: Vector[Vector[Int]], z: MatrixEnclosure, rhoLower: Double, rhoUpper: Double):
  require(id.nonEmpty && rhoLower.isFinite && rhoUpper.isFinite && rhoLower >= 0.0 && rhoLower <= rhoUpper)
enum TaskRankMethod:
  case FrobeniusFeasibility, OperatorFirstThenBound
enum TaskRankDecision:
  case Counted, Excluded, Ambiguous
final case class TaskRankComparison(rank: Int, decision: TaskRankDecision, tailLower: Double, tailUpper: Double, comparatorLower: Double, comparatorUpper: Double)
final case class TaskRankReceipt(
    method: TaskRankMethod, coordinates: TaskRankCoordinates, exact: Boolean,
    drawIds: Vector[String], signs: Vector[Vector[Vector[Int]]], comparisons: Vector[Vector[TaskRankComparison]],
    counts: Vector[BigInt], denominator: BigInt, stoppedAt: Option[Int]
)

object TaskMeanRank:
  /** Evaluates all supplied signs before returning any receipt.  A failed SVD
    * certificate therefore cannot leak partial inferential evidence.
    */
  def fixedFrobenius(coordinates: TaskRankCoordinates, draws: Vector[TaskRankDraw], alphaNumerator: BigInt, alphaDenominator: BigInt): Either[InferenceError, TaskRankReceipt] =
    if draws.isEmpty || alphaNumerator < 0 || alphaDenominator <= 0 || alphaNumerator > alphaDenominator then
      Left(InferenceError.InvalidReplicatePlan("task-rank requires nonempty fixed draws and an alpha rational in [0,1]"))
    else
      val all = Vector.newBuilder[Vector[TaskRankComparison]]
      var failure: Option[InferenceError] = None
      var i = 0
      while i < draws.length && failure.isEmpty do
        evaluateFrobenius(draws(i)) match
          case Left(error) => failure = Some(error)
          case Right(value) => all += value
        i += 1
      failure match
        case Some(error) => Left(error)
        case None =>
          val rows = all.result()
          val q = rows.head.length
          // Ambiguity is conservative counting.  Cumulative max makes finite
          // certificates monotone even if independently rounded draw bounds vary.
          var running = BigInt(0)
          val counts = Vector.tabulate(q) { r =>
            val raw = BigInt(rows.count(_(r).decision != TaskRankDecision.Excluded))
            running = running.max(raw); running
          }
          val denominator = BigInt(draws.size + 1) // identity-plus-B convention
          val stopped = counts.indices.find(r => counts(r) * alphaDenominator > alphaNumerator * denominator)
          Right(TaskRankReceipt(TaskRankMethod.FrobeniusFeasibility, coordinates, exact = false,
            draws.map(_.id), draws.map(_.signs), rows, counts, denominator, stopped))

  private def evaluateFrobenius(draw: TaskRankDraw): Either[InferenceError, Vector[TaskRankComparison]] =
    if draw.signs.exists(_.exists(s => s != -1 && s != 1)) then Left(InferenceError.InvalidReplicatePlan("task-rank signs must be exact -1/+1 values"))
    else
      // A candidate SVD is never trusted directly: Gale checks factor residual,
      // orthogonality and every singular-value enclosure.
      import gale.backend.Backend.given
      gale.spectral.Svds.svd(draw.z.lower, gale.spectral.SingularSelection.All) match
        case Left(error) => Left(InferenceError.NumericalFailure("task-rank candidate SVD", error.toString))
        case Right(candidate) => ValidatedSvd.enclosure(draw.z, candidate).left.map(e => InferenceError.NumericalFailure("task-rank SVD certificate", e.toString)).map { spectrum =>
          Vector.tabulate(spectrum.lower.length) { r =>
            val tail = tailBounds(spectrum.lower, spectrum.upper, r)
            val decision = if tail._2 <= draw.rhoLower then TaskRankDecision.Counted
              else if tail._1 > draw.rhoUpper then TaskRankDecision.Excluded else TaskRankDecision.Ambiguous
            TaskRankComparison(r, decision, tail._1, tail._2, draw.rhoLower, draw.rhoUpper)
          }
        }

  private def tailBounds(lower: DVec, upper: DVec, rank: Int): (Double, Double) =
    var lo = 0.0; var hi = 0.0; var i = rank
    while i < lower.length do
      lo += lower(i) * lower(i); hi += upper(i) * upper(i); i += 1
    lo -> hi

package multivar
package family.paired

import multivar.core.*

/** How retained components are allocated across focus and remainder. */
sealed trait RankBudget:
  def label: String

object RankBudget:
  /** Retain components only from the scientific focus. Ordinary `components = k` means this. */
  final case class Focus(components: Int) extends RankBudget:
    val label: String = "focus"

  /** One total budget allocated by singular-value magnitude across both parts. */
  final case class Total(components: Int) extends RankBudget:
    val label: String = "total"

  /** Explicit ranks for the two orthogonal parts. */
  final case class Split(focus: Int, remainder: Int) extends RankBudget:
    val label: String = "split"

  private[paired] def validate(budget: RankBudget): Either[MultivarError, RankBudget] =
    budget match
      case Focus(components) =>
        ComponentCount(components).map(_ => budget)
      case Total(components) =>
        ComponentCount(components).map(_ => budget)
      case Split(focus, remainder) =>
        for
          _ <-
            if focus < 0 || remainder < 0 then
              Left(MultivarError.InvalidDimension("rank split", focus + remainder))
            else if focus + remainder <= 0 then
              Left(MultivarError.InvalidDimension("rank split total", focus + remainder))
            else Right(())
          _ <- if focus > 0 then ComponentCount(focus).map(_ => ()) else Right(())
          _ <- if remainder > 0 then ComponentCount(remainder).map(_ => ()) else Right(())
        yield budget

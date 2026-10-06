package multivar.inference

import gale.linalg.DMat
import resample4s.kernel.{Seed, StreamDomain}

/** Exact row multipliers, independent of any null hypothesis or validity claim. */
final class RowSignFlip private (
    val rowCount: RowCount,
    val signs: Vector[Int]
):
  // Scala's private constructor is visible in JVM bytecode; keep Java callers checked too.
  require(rowCount.value > 0 && signs.size == rowCount.value && signs.forall(sign => sign == -1 || sign == 1),
    "Row sign multipliers must match a positive row count and contain only -1 or +1")
  def applyTo(input: DMat): Either[InferenceError, DMat] =
    if input.rows != rowCount.value then
      Left(InferenceError.RowCountMismatch("row sign-flip input", rowCount.value, input.rows))
    else Right(InferenceNumerics.matrixFromRows(Vector.tabulate(input.rows) { row =>
      Vector.tabulate(input.cols)(column => signs(row).toDouble * input(row, column))
    }))

  def compose(other: RowSignFlip): Either[InferenceError, RowSignFlip] =
    if other.rowCount != rowCount then
      Left(InferenceError.RowCountMismatch("row sign-flip composition", rowCount.value, other.rowCount.value))
    else RowSignFlip.from(signs.zip(other.signs).map(_ * _))

object RowSignFlip:
  def from(signs: Vector[Int]): Either[InferenceError, RowSignFlip] =
    for
      rows <- RowCount(signs.size)
      _ <- Either.cond(signs.forall(sign => sign == -1 || sign == 1), (),
        InferenceError.UnsupportedProblem("Row sign multipliers must be -1 or +1"))
    yield new RowSignFlip(rows, signs)

/** One independent sign per partition group, shared by every row in that group.
  * Binding the action establishes row structure, not distributional symmetry.
  */
final class SignFlipAction private (val partition: RowPartition):
  def rowCount: RowCount = partition.rowCount
  def draw(seed: Seed, replicate: ReplicateId): Either[InferenceError, RowSignFlip] =
    for
      binaryRows <- RowCount(2)
      domain <- StreamDomain.custom(1003).left.map(error => InferenceError.UnsupportedProblem(error.message))
      permutations <- ResamplingPlans.independent(PermutationAction.unrestricted(binaryRows),
        seed, replicate, partition.groups.size, domain)
      groupSigns = permutations.map(p => if p.toIArray(0) == 0 then 1 else -1)
      pattern <- RowSignFlip.from(partition.groupIndexByRow.toVector.map(groupSigns))
    yield pattern

object SignFlipAction:
  def forPartition(partition: RowPartition): SignFlipAction = new SignFlipAction(partition)
  def independent(rows: RowCount): Either[InferenceError, SignFlipAction] =
    RowPartition.from(rows, (0 until rows.value).map(row => Vector(row))).map(forPartition)

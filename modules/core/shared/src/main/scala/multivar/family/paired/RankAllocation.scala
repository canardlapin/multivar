package multivar
package family.paired

/** Internal rank-allocation result.
  *
  * The class is public only because Scala 3.7 ScalaDoc cannot render the
  * package-private constructor signature. It is not exported from
  * `multivar.analysis`.
  */
final class RankAllocation private[paired] (
    val focus: Int,
    val remainder: Int,
    val boundaryTied: Boolean,
    val order: Vector[TaskComponentIndex]
):
  def total: Int =
    focus + remainder

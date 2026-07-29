package multivar
package family.paired

/** Stable index of one retained component across focus and remainder. */
final case class TaskComponentIndex(
    globalIndex: Int,
    part: TaskPart,
    localIndex: Int,
    singularValue: Double
)

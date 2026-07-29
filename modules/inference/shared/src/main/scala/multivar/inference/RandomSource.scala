package multivar.inference

final case class ReplicateAssignment private[inference] (
    step: ComponentIx,
    withinStep: Int,
    replicate: ReplicateId
)

object ReplicatePlan:
  def forStep(
      step: ComponentIx,
      capacity: MonteCarloDraws,
      draws: MonteCarloDraws
  ): Either[InferenceError, Vector[ReplicateAssignment]] =
    if draws.value > capacity.value then
      Left(InferenceError.InvalidReplicatePlan(
        s"requested ${draws.value} draws from step capacity ${capacity.value}"
      ))
    else
      val start = step.value.toLong * capacity.value.toLong
      val last = start + draws.value.toLong - 1L
      if last > Int.MaxValue then
        Left(InferenceError.InvalidReplicatePlan(
          s"replicate ids overflow Int at step ${step.value} with capacity ${capacity.value}"
        ))
      else
        Right(Vector.tabulate(draws.value) { within =>
          ReplicateAssignment(step, within, acceptedReplicate((start + within).toInt))
        })

  private def acceptedReplicate(value: Int): ReplicateId =
    ReplicateId(value) match
      case Right(replicate) => replicate
      case Left(error)      => throw IllegalStateException(error.message)

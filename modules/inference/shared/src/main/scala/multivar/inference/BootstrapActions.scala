package multivar.inference

import resample4s.Bootstrap
import resample4s.kernel.{Draw, Seed}

sealed trait BootstrapAction:
  def rowCount: RowCount
  def draw(seed: Seed, replicate: ReplicateId): Either[InferenceError, Draw]

object BootstrapAction:
  final case class Rows private[inference] (rowCount: RowCount) extends BootstrapAction:
    override def draw(
        seed: Seed,
        replicate: ReplicateId
    ): Either[InferenceError, Draw] =
      for
        replicateSeed <- ResamplingPlans.replicateSeed(seed, replicate)
        plan <- adapt(Bootstrap
          .unconditional(1)
          .plan(rowCount.value, replicateSeed.value))
      yield plan.compiled.plan.first.analysis

  final case class Clusters private[inference] (
      partition: ClusterPartition
  ) extends BootstrapAction:
    override def rowCount: RowCount = partition.rowCount

    override def draw(
        seed: Seed,
        replicate: ReplicateId
    ): Either[InferenceError, Draw] =
      for
        replicateSeed <- ResamplingPlans.replicateSeed(seed, replicate)
        plan <- adapt(Bootstrap
          .grouped(1, partition.groups)
          .plan(rowCount.value, replicateSeed.value))
      yield plan.compiled.plan.first.analysis

  def rows(count: RowCount): BootstrapAction =
    Rows(count)

  def clusters(partition: ClusterPartition): BootstrapAction =
    Clusters(partition)

  private def adapt[A](
      value: Either[resample4s.spi.DesignError, A]
  ): Either[InferenceError, A] =
    value.left.map(InferenceError.Resampling.apply)

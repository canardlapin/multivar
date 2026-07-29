package multivar.inference

import resample4s.PermutationTest
import resample4s.kernel.*
import resample4s.spi.AlgorithmId

private[inference] object ResamplingPlans:
  val permutationAlgorithm: AlgorithmId =
    AlgorithmId
      .of("permutation/v1")
      .fold(
        error => throw IllegalStateException(error.message),
        (value: AlgorithmId) => value
      )

  private val ReplicateDomain =
    StreamDomain
      .custom(1000)
      .fold(
        error => throw IllegalStateException(error.message),
        (value: StreamDomain) => value
      )

  private[inference] val PcaColumnDomain =
    StreamDomain
      .custom(1001)
      .fold(
        error => throw IllegalStateException(error.message),
        (value: StreamDomain) => value
      )

  private[inference] val MultiblockDomain =
    StreamDomain
      .custom(1002)
      .fold(
        error => throw IllegalStateException(error.message),
        (value: StreamDomain) => value
      )

  def permutation(
      root: Seed,
      replicate: ReplicateId,
      size: Int
  ): Either[InferenceError, Permutation] =
    for
      seed <- replicateSeed(root, replicate)
      permutation <- compilePermutation(seed, size)
    yield permutation

  def replicateSeed(
      root: Seed,
      replicate: ReplicateId
  ): Either[InferenceError, Seed] =
    adapt(StreamPath.of(ReplicateDomain, replicate.value))
      .map(root.derive)

  def pcaColumns(
      root: Seed,
      replicate: ReplicateId,
      rows: Int,
      columns: Int
  ): Either[InferenceError, Vector[Permutation]] =
    independent(root, replicate, rows, columns, PcaColumnDomain)

  def multiblockLanes(
      root: Seed,
      replicate: ReplicateId,
      rows: Int,
      lanes: Int
  ): Either[InferenceError, Vector[Permutation]] =
    independent(root, replicate, rows, lanes, MultiblockDomain)

  def independent(
      action: PermutationAction,
      root: Seed,
      replicate: ReplicateId,
      lanes: Int,
      laneDomain: StreamDomain
  ): Either[InferenceError, Vector[Permutation]] =
    val out = Vector.newBuilder[Permutation]
    var lane = 0
    var failure = Option.empty[InferenceError]
    while lane < lanes && failure.isEmpty do
      val current =
        for
          path <- adapt(StreamPath.of(laneDomain, lane))
          permutation <- action.draw(root.derive(path), replicate)
        yield permutation
      current match
        case Left(error)  => failure = Some(error)
        case Right(value) => out += value
      lane += 1
    failure.toLeft(out.result())

  private def independent(
      root: Seed,
      replicate: ReplicateId,
      size: Int,
      lanes: Int,
      laneDomain: StreamDomain
  ): Either[InferenceError, Vector[Permutation]] =
    for
      replicatePath <- adapt(StreamPath.of(ReplicateDomain, replicate.value))
      permutations <-
        val out = Vector.newBuilder[Permutation]
        var lane = 0
        var failure = Option.empty[InferenceError]
        while lane < lanes && failure.isEmpty do
          val current =
            for
              path <- adapt(replicatePath.append(laneDomain, lane))
              permutation <- compilePermutation(root.derive(path), size)
            yield permutation
          current match
            case Left(error)  => failure = Some(error)
            case Right(value) => out += value
          lane += 1
        failure.toLeft(out.result())
    yield permutations

  private def compilePermutation(
      seed: Seed,
      size: Int
  ): Either[InferenceError, Permutation] =
    adapt(PermutationTest(1).plan(size, seed.value))
      .map(_.compiled.plan.first)

  private def adapt[A](
      value: Either[resample4s.spi.DesignError, A]
  ): Either[InferenceError, A] =
    value.left.map(InferenceError.Resampling.apply)

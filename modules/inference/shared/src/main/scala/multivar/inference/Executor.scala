package multivar.inference

import resample4s.kernel.Seed
import resample4s.spi.AlgorithmId

trait ExactLadderProtocol[
    S,
    F,
    K <: TargetKind,
    N <: NullKind
]:
  type Action

  def fit: FitDescriptor[F]
  def target: TargetSpec[K]
  def nullHypothesis: NullSpec[N]
  def validity: ValidityClaim
  def randomizationAlgorithm: AlgorithmId =
    ResamplingPlans.permutationAlgorithm
  def validateDesign(
      state: S,
      design: ResamplingDesign[?]
  ): Either[InferenceError, Unit]
  def roots(initial: S): Either[InferenceError, Vector[Double]]
  def observed(state: S): Either[InferenceError, Double]
  def action(
      state: S,
      replicate: ReplicateId,
      seed: Seed
  ): Either[InferenceError, Action]
  def actionForDesign(
      state: S,
      design: ResamplingDesign[?],
      replicate: ReplicateId,
      seed: Seed
  ): Either[InferenceError, Action]
  def nullStatistic(
      state: S,
      step: ComponentIx,
      replicate: ReplicateId,
      action: Action
  ): Either[InferenceError, Double]
  def remove(state: S): Either[InferenceError, S]

final case class LadderExecutionConfig(
    perRung: MonteCarloDraws,
    totalBudget: MonteCarloDraws,
    alpha: Alpha,
    batchSize: BatchSize,
    maxSteps: LadderSteps,
    seed: Seed,
    unitPolicy: UnitPolicy = UnitPolicy.SingleAxes,
    boundary: SequentialBoundary = SequentialBoundary.FromAlpha
)

final case class ExactLadderRun[S, K <: TargetKind](
    target: TargetSpec[K],
    roots: Vector[Double],
    units: Vector[LatentUnitResult],
    ladder: LadderResult,
    finalState: S,
    validity: ValidityClaim
)

object InferenceExecutor:
  def runLadder[S, F, K <: TargetKind, N <: NullKind](
      initial: S,
      protocol: ExactLadderProtocol[S, F, K, N],
      config: LadderExecutionConfig
  ): Either[InferenceError, ExactLadderRun[S, K]] =
    protocol.roots(initial).flatMap { roots =>
      val limit = Math.min(config.maxSteps.value, roots.length)
      if limit <= 0 then Left(InferenceError.InvalidCount("available ladder roots", limit))
      else
        requireAxisStatistics(roots, config.unitPolicy).flatMap: _ =>
          execute(
            initial,
            protocol,
            roots,
            limit,
            None,
            config,
            protocol.validity
          )
    }

  def runProgram[
      F,
      S,
      K <: TargetKind,
      N <: NullKind,
      D <: DesignKind
  ](
      initial: S,
      program: InferenceProgram[F, K, N, D],
      protocol: ExactLadderProtocol[S, F, K, N],
      config: LadderExecutionConfig
  ): Either[InferenceError, ExactLadderRun[S, K]] =
    for
      _ <-
        if program.spec.fit.label == protocol.fit.label then Right(())
        else Left(InferenceError.FitProtocolMismatch(
          program.spec.fit,
          protocol.fit
        ))
      _ <-
        if program.spec.target.label == protocol.target.label then Right(())
        else Left(InferenceError.TargetProtocolMismatch(
          program.spec.target.label,
          protocol.target.label
        ))
      _ <-
        if program.spec.nullHypothesis.label == protocol.nullHypothesis.label then Right(())
        else Left(InferenceError.NullProtocolMismatch(
          program.spec.nullHypothesis.label,
          protocol.nullHypothesis.label
        ))
      _ <- protocol.validateDesign(initial, program.spec.design)
      _ <- program.spec.evidence match
        case RequestedEvidence.SignificanceOnly => Right(())
        case other => Left(InferenceError.UnsupportedEvidence(other))
      _ <- validateMonteCarlo(program.spec.monteCarlo, config)
      roots <- protocol.roots(initial)
      limit = Math.min(config.maxSteps.value, roots.length)
      _ <-
        if limit > 0 then Right(())
        else Left(InferenceError.InvalidCount("available ladder roots", limit))
      _ <- requireAxisStatistics(roots, program.spec.units)
      result <- execute(
        initial,
        protocol,
        roots,
        limit,
        Some(program.spec.design),
        config.copy(
          seed = program.spec.seed,
          unitPolicy = program.spec.units
        ),
        program.validity
      )
    yield result

  private def requireAxisStatistics(
      roots: Vector[Double],
      policy: UnitPolicy
  ): Either[InferenceError, Unit] =
    LatentUnitFormation
      .form(roots, Vector.fill(roots.length)(false), policy)
      .flatMap: units =>
        if units.exists(_.unit.identifiability == Identifiability.UnorientedSubspace) then
          Left(InferenceError.UnsupportedUnitPolicy(policy))
        else Right(())

  private def validateMonteCarlo(
      policy: MonteCarloPolicy,
      config: LadderExecutionConfig
  ): Either[InferenceError, Unit] =
    policy match
      case MonteCarloPolicy.Fixed(_) =>
        Left(InferenceError.ExecutionPolicyConflict(
          ExecutionPolicyMismatch.FixedMonteCarloUnsupported
        ))
      case MonteCarloPolicy.Sequential(maxDraws, alpha, batchSize, boundary) =>
        if maxDraws != config.perRung then
          Left(InferenceError.ExecutionPolicyConflict(
            ExecutionPolicyMismatch.PerRungDraws
          ))
        else if alpha != config.alpha then
          Left(InferenceError.ExecutionPolicyConflict(
            ExecutionPolicyMismatch.Alpha
          ))
        else if batchSize != config.batchSize then
          Left(InferenceError.ExecutionPolicyConflict(
            ExecutionPolicyMismatch.BatchSize
          ))
        else if boundary != config.boundary then
          Left(InferenceError.ExecutionPolicyConflict(
            ExecutionPolicyMismatch.SequentialBoundary
          ))
        else Right(())

  private def execute[S, F, K <: TargetKind, N <: NullKind](
      initial: S,
      protocol: ExactLadderProtocol[S, F, K, N],
      roots: Vector[Double],
      limit: Int,
      design: Option[ResamplingDesign[?]],
      config: LadderExecutionConfig,
      validity: ValidityClaim
  ): Either[InferenceError, ExactLadderRun[S, K]] =
    var state = initial
    var budget = BudgetState.initial(config.totalBudget)
    val steps = Vector.newBuilder[LadderStepResult]
    val batches = Vector.newBuilder[Int]
    var rejected = 0
    var index = 0
    var termination = Option.empty[LadderTermination]

    while index < limit && termination.isEmpty do
      val component = acceptedComponent(index)
      val unit = LatentUnit.Axis(UnitId.unsafe(s"u${index + 1}"), component)
      val observed = protocol.observed(state) match
        case Right(value) => value
        case Left(error)  => return Left(error)

      budget.checkout(config.perRung) match
        case Left(_: InferenceError.BudgetExhausted) =>
          steps += LadderStepResult(
            unit,
            observed,
            Evidence.Unavailable(UnavailableReason.Unsupported("global Monte Carlo budget exhausted")),
            selected = false
          )
          termination = Some(LadderTermination.BudgetExhausted(unit.id))
        case Left(error) => return Left(error)
        case Right(grant) =>
          val assignments = ReplicatePlan.forStep(component, config.perRung, grant.allocated) match
            case Right(value) => value
            case Left(error)  => return Left(error)
          var accumulator = MonteCarlo.sequentialAccumulator(
            observed,
            protocol.target.alternative,
            grant.allocated,
            config.alpha,
            config.batchSize,
            config.boundary,
            ReplicateProvenance.Deterministic(
              config.seed,
              protocol.randomizationAlgorithm,
              Seed.derivationAlgorithm
            )
          ) match
            case Right(value) => value
            case Left(error)  => return Left(error)
          while !accumulator.complete do
            val thisBatch = Math.min(
              config.batchSize.value,
              grant.allocated.value - accumulator.consumed
            )
            val statistics = Vector.newBuilder[ReplicateStatistic]
            var draw = accumulator.consumed
            val end = draw + thisBatch
            while draw < end do
              val assignment = assignments(draw)
              val statistic =
                for
                  action <- design match
                    case Some(value) =>
                      protocol.actionForDesign(
                        state,
                        value,
                        assignment.replicate,
                        config.seed
                      )
                    case None =>
                      protocol.action(
                        state,
                        assignment.replicate,
                        config.seed
                      )
                  value <- protocol.nullStatistic(
                    state,
                    component,
                    assignment.replicate,
                    action
                  )
                  statistic <- ReplicateStatistic.from(
                    assignment.replicate,
                    value
                  )
                yield statistic
              statistic match
                case Left(error) =>
                  return Left(
                    InferenceError.ReplicateFailure(
                      assignment.replicate,
                      error
                    )
                  )
                case Right(value) => statistics += value
              draw += 1
            accumulator.offer(statistics.result()) match
              case Right(value) => accumulator = value
              case Left(error)  => return Left(error)
          val receipt = accumulator.result match
            case Right(value) => value
            case Left(error)  => return Left(error)

          budget.record(grant, receipt.consumed.value) match
            case Right(value) => budget = value
            case Left(error)  => return Left(error)
          batches ++= receipt.batchSchedule
          val selected = receipt.pValue.value <= config.alpha.value
          steps += LadderStepResult(unit, observed, Evidence.Computed(receipt), selected)
          if selected then
            rejected += 1
            if index + 1 < limit then
              protocol.remove(state) match
                case Right(value) => state = value
                case Left(error)  => return Left(error)
          else termination = Some(LadderTermination.FirstNonSelection(unit.id))
      index += 1

    val ladder = LadderResult(
      steps.result(),
      rejected,
      budget.used,
      batches.result(),
      budget,
      termination.getOrElse(LadderTermination.Completed)
    )
    val selected = Vector.tabulate(roots.length)(index => index < rejected)
    LatentUnitFormation.form(roots, selected, config.unitPolicy).map { units =>
      ExactLadderRun(protocol.target, roots, units, ladder, state, validity)
    }

  private def acceptedComponent(value: Int): ComponentIx =
    ComponentIx(value) match
      case Right(component) => component
      case Left(error)      => throw IllegalStateException(error.message)

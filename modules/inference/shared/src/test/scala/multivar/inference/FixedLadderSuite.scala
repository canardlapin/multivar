package multivar.inference

import resample4s.kernel.Seed

class FixedLadderSuite extends munit.FunSuite:
  private def accepted[A](value: Either[InferenceError, A]): A = value.fold(error => fail(error.message), identity)
  private val draws = accepted(MonteCarloDraws(4))
  private def config(total: Int = 12): FixedLadderExecutionConfig =
    FixedLadderExecutionConfig(draws, accepted(MonteCarloDraws(total)), accepted(Alpha(.25)),
      accepted(LadderSteps(3)), Seed.fromLong(44))

  private class Protocol(
      failAt: Option[Int] = None,
      nonfinite: Boolean = false,
      validityClaim: ValidityClaim = ValidityClaim.Heuristic
  )
      extends ExactLadderProtocol[Int, PlscFitFamily, TargetKind.CovarianceRoots, NullKind.PairedIndependence]:
    type Action = Int
    val fit: FitDescriptor[PlscFitFamily] = FitDescriptor.Plsc
    val target: TargetSpec[TargetKind.CovarianceRoots] = TargetSpec.CovarianceRoots
    val nullHypothesis: NullSpec[NullKind.PairedIndependence] = NullSpec.BreakX
    def validity: ValidityClaim = validityClaim
    var calls = Vector.empty[Int]
    var seeds = Vector.empty[Seed]
    var removals = 0
    var designedCalls = 0
    def validateDesign(state: Int, design: ResamplingDesign[?]): Either[InferenceError, Unit] = Right(())
    def roots(initial: Int): Either[InferenceError, Vector[Double]] = Right(Vector(10,8,6))
    def observed(state: Int): Either[InferenceError, Double] = Right(if nonfinite then Double.NaN else 10.0 - 2 * state)
    def action(state: Int, replicate: ReplicateId, seed: Seed): Either[InferenceError, Int] =
      calls :+= replicate.value
      seeds :+= seed
      Right(replicate.value)
    def actionForDesign(state: Int, design: ResamplingDesign[?], replicate: ReplicateId, seed: Seed): Either[InferenceError, Int] =
      designedCalls += 1
      action(state, replicate, seed)
    def nullStatistic(state: Int, step: ComponentIx, replicate: ReplicateId, action: Int): Either[InferenceError, Double] =
      if failAt.contains(action) then Left(InferenceError.UnsupportedProblem("injected draw failure"))
      else Right(if state == 0 then 0.0 else 20.0)
    def remove(state: Int): Either[InferenceError, Int] =
      removals += 1
      Right(state + 1)

  private def receipts(run: ExactLadderRun[Int, TargetKind.CovarianceRoots]): Vector[MonteCarloReceipt] =
    run.ladder.steps.map(_.evidence).collect { case Evidence.Computed(value) => value }

  test("fixed rungs consume full allocation and stop after first nonrejection"):
    val protocol = new Protocol()
    val run = accepted(InferenceExecutor.runLadder(0, protocol, config()))
    assertEquals(protocol.calls, (0 until 8).toVector)
    assertEquals(protocol.removals, 1)
    assertEquals(run.finalState, 1)
    assertEquals(run.validity, ValidityClaim.Heuristic)
    assertEquals(run.ladder.rejectedThrough, 1)
    assertEquals(run.ladder.totalDraws, 8)
    assertEquals(run.ladder.termination, LadderTermination.FirstNonSelection(UnitId.unsafe("u2")))
    val evidence = receipts(run)
    assertEqualsDouble(evidence(0).pValue.value, .2, 1e-12)
    assertEqualsDouble(evidence(1).pValue.value, 1.0, 1e-12)
    evidence.foreach { receipt =>
      assertEquals(receipt.allocated.value, 4)
      assertEquals(receipt.consumed.value, 4)
      assertEquals(receipt.boundary, None)
      assertEquals(receipt.stopReason, MonteCarloStopReason.Exhausted)
      assertEquals(receipt.batchSchedule, Vector(4))
    }

  test("insufficient global remainder does not silently shorten a fixed rung"):
    val protocol = new Protocol()
    val run = accepted(InferenceExecutor.runLadder(0, protocol, config(6)))
    assertEquals(protocol.calls, Vector(0,1,2,3))
    assertEquals(run.ladder.totalDraws, 4)
    assertEquals(run.ladder.budget.remaining, 2)
    assertEquals(run.ladder.termination, LadderTermination.BudgetExhausted(UnitId.unsafe("u2")))
    assertEquals(receipts(run).size, 1)

  test("a budget smaller than the first rung generates no draws"):
    val protocol = new Protocol()
    val run = accepted(InferenceExecutor.runLadder(0, protocol, config(3)))
    assertEquals(protocol.calls, Vector.empty)
    assertEquals(run.ladder.totalDraws, 0)
    assertEquals(run.ladder.termination, LadderTermination.BudgetExhausted(UnitId.unsafe("u1")))

  test("draw failures and nonfinite observations fail without partial evidence"):
    val failing = new Protocol(failAt = Some(2))
    InferenceExecutor.runLadder(0, failing, config()) match
      case Left(InferenceError.ReplicateFailure(id, _)) => assertEquals(id.value, 2)
      case other => fail(s"expected replicate failure, got $other")
    assertEquals(failing.calls, Vector(0,1,2))
    val nonfinite = new Protocol(nonfinite = true)
    assert(InferenceExecutor.runLadder(0, nonfinite, config()).isLeft)
    assertEquals(nonfinite.calls, Vector.empty)

  test("compiled fixed policy uses design action and program seed; mismatched mode refuses"):
    val program = InferenceCompiler.compile(InferenceSpec.of(
      FitDescriptor.Plsc, TargetSpec.CovarianceRoots, NullSpec.BreakX,
      ResamplingDesign.exchangeableRows(accepted(RowCount(4))), UnitPolicy.SingleAxes,
      MonteCarloPolicy.Fixed(draws), RequestedEvidence.SignificanceOnly, Seed.fromLong(81)))
    val protocol = new Protocol()
    val run = accepted(InferenceExecutor.runProgram(0, program, protocol, config()))
    assertEquals(run.ladder.totalDraws, 8)
    assertEquals(protocol.designedCalls, 8)
    assert(protocol.seeds.forall(_ == Seed.fromLong(81)))
    val sequential = LadderExecutionConfig(draws, draws, accepted(Alpha(.25)),
      accepted(BatchSize(2)), accepted(LadderSteps(1)), Seed.fromLong(44))
    assertEquals(InferenceExecutor.runProgram(0, program, new Protocol(), sequential),
      Left(InferenceError.ExecutionPolicyConflict(ExecutionPolicyMismatch.MonteCarloMode)))
    assertEquals(InferenceExecutor.runProgram(0, program, new Protocol(), config().copy(perRung = accepted(MonteCarloDraws(5)))),
      Left(InferenceError.ExecutionPolicyConflict(ExecutionPolicyMismatch.PerRungDraws)))

  test("compiled and protocol validity claims compose explicitly before draws"):
    val fixedProgram = InferenceCompiler.compile(InferenceSpec.of(
      FitDescriptor.Plsc, TargetSpec.CovarianceRoots, NullSpec.BreakX,
      ResamplingDesign.exchangeableRows(accepted(RowCount(4))), UnitPolicy.SingleAxes,
      MonteCarloPolicy.Fixed(draws), RequestedEvidence.SignificanceOnly, Seed.fromLong(81)))
    val fixed = new Protocol(validityClaim = ValidityClaim.Heuristic)
    assertEquals(accepted(InferenceExecutor.runProgram(0, fixedProgram, fixed, config())).validity,
      ValidityClaim.Heuristic)

    val batch = accepted(BatchSize(2))
    val sequentialProgram = InferenceCompiler.compile(InferenceSpec.of(
      FitDescriptor.Plsc, TargetSpec.CovarianceRoots, NullSpec.BreakX,
      ResamplingDesign.exchangeableRows(accepted(RowCount(4))), UnitPolicy.SingleAxes,
      MonteCarloPolicy.Sequential(draws, accepted(Alpha(.25)), batch, SequentialBoundary.FromAlpha),
      RequestedEvidence.SignificanceOnly, Seed.fromLong(81)))
    val sequentialConfig = LadderExecutionConfig(draws, accepted(MonteCarloDraws(12)), accepted(Alpha(.25)),
      batch, accepted(LadderSteps(3)), Seed.fromLong(44))
    val sequential = new Protocol(validityClaim = ValidityClaim.Heuristic)
    assertEquals(accepted(InferenceExecutor.runProgram(0, sequentialProgram, sequential, sequentialConfig)).validity,
      ValidityClaim.Heuristic)

    val invalidFixed = new InferenceProgram(fixedProgram.spec, fixedProgram.summary, ValidityClaim.Conditional)
    val refusedFixed = new Protocol(validityClaim = ValidityClaim.Heuristic)
    InferenceExecutor.runProgram(0, invalidFixed, refusedFixed, config()) match
      case Left(InferenceError.InvalidValidity(_)) => ()
      case other => fail(s"expected validity refusal, got $other")
    assertEquals(refusedFixed.calls, Vector.empty)

    val invalidSequential = new InferenceProgram(sequentialProgram.spec, sequentialProgram.summary, ValidityClaim.Conditional)
    val refusedSequential = new Protocol(validityClaim = ValidityClaim.Heuristic)
    InferenceExecutor.runProgram(0, invalidSequential, refusedSequential, sequentialConfig) match
      case Left(InferenceError.InvalidValidity(_)) => ()
      case other => fail(s"expected validity refusal, got $other")
    assertEquals(refusedSequential.calls, Vector.empty)

  test("existing adaptive execution still stops within the nonrejecting rung"):
    val protocol = new Protocol()
    val sequential = LadderExecutionConfig(draws, accepted(MonteCarloDraws(12)), accepted(Alpha(.25)),
      accepted(BatchSize(1)), accepted(LadderSteps(3)), Seed.fromLong(44))
    val run = accepted(InferenceExecutor.runLadder(0, protocol, sequential))
    assertEquals(run.ladder.totalDraws, 6)
    assertEquals(protocol.calls, (0 until 6).toVector)
    assertEquals(receipts(run)(1).stopReason, MonteCarloStopReason.NonRejectionEarly)
    assertEquals(receipts(run)(1).consumed.value, 2)

  test("step limit completes after a selected fixed rung without unnecessary removal"):
    val protocol = new Protocol()
    val run = accepted(InferenceExecutor.runLadder(0, protocol, config().copy(maxSteps = accepted(LadderSteps(1)))))
    assertEquals(run.ladder.termination, LadderTermination.Completed)
    assertEquals(run.ladder.totalDraws, 4)
    assertEquals(protocol.removals, 0)
    assertEquals(run.finalState, 0)

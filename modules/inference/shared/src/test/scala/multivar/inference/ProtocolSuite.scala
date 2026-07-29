package multivar.inference

import gale.linalg.DMat
import resample4s.kernel.{Permutation, Seed}

class ProtocolSuite extends munit.FunSuite:

  import InferenceRReferenceFixtures as R

  private def accepted[A](value: Either[InferenceError, A]): A =
    value.fold(error => fail(error.message), identity)

  private def matrix(value: R.MatrixData): DMat =
    InferenceNumerics.matrixFromRows(Vector.tabulate(value.rows) { row =>
      Vector.tabulate(value.cols)(col => value(row, col))
    })

  private def negate(value: DMat): DMat =
    InferenceNumerics.matrixFromRows(Vector.tabulate(value.rows) { row =>
      Vector.tabulate(value.cols)(col => -value(row, col))
    })

  private def reverseColumns(value: DMat): DMat =
    InferenceNumerics.matrixFromRows(Vector.tabulate(value.rows) { row =>
      Vector.tabulate(value.cols)(col => value(row, value.cols - 1 - col))
    })

  private val config = LadderExecutionConfig(
    accepted(MonteCarloDraws(39)),
    accepted(MonteCarloDraws(117)),
    accepted(Alpha(0.1)),
    accepted(BatchSize(5)),
    accepted(LadderSteps(3)),
    Seed.fromLong(1729L),
    UnitPolicy.GroupNearTies(accepted(RelativeGap(0.01)))
  )

  test("PCA protocol reproduces observed roots and deflated rung statistics") {
    val fixture = R.ladders.find(_.family == "pca").get
    val protocol = PcaVarianceProtocol()
    val initial = accepted(PcaVarianceState.from(matrix(fixture.x)))
    val roots = accepted(protocol.roots(initial))
    val first = accepted(protocol.observed(initial))
    val secondState = accepted(protocol.remove(initial))
    val second = accepted(protocol.observed(secondState))

    fixture.roots.zip(roots).foreach { case (expected, actual) =>
      assertEqualsDouble(actual, expected, 1e-8)
    }
    assertEqualsDouble(first, fixture.steps(0).observed, 1e-12)
    assertEqualsDouble(second, fixture.steps(1).observed, 1e-12)

    val identity =
      Permutation.from(IArray.from(0 until secondState.residual.rows)) match
        case Right(value) => value
        case Left(error)  => fail(error.message)
    val identityNull = accepted(
      protocol.nullStatistic(
        secondState,
        accepted(ComponentIx(1)),
        accepted(ReplicateId(0)),
        Vector.fill(secondState.residual.cols)(identity)
      )
    )
    assertEqualsDouble(identityNull, second, 1e-12)
  }

  test("PLSC protocol reproduces observed cross-covariance roots and deflation") {
    val fixture = R.ladders.find(_.family == "plsc").get
    val protocol = PlscCovarianceProtocol()
    val initial = accepted(PlscCovarianceState.from(matrix(fixture.x), matrix(fixture.y.get)))
    val roots = accepted(protocol.roots(initial))
    val first = accepted(protocol.observed(initial))
    val secondState = accepted(protocol.remove(initial))
    val second = accepted(protocol.observed(secondState))

    fixture.roots.zip(roots).foreach { case (expected, actual) =>
      assertEqualsDouble(actual, expected, 1e-7)
    }
    assertEqualsDouble(first, fixture.steps(0).observed, 1e-8)
    assertEqualsDouble(second, fixture.steps(1).observed, 1e-8)
  }

  test("PCA and PLSC targets are invariant to sign and orthogonal coordinate changes") {
    val pcaFixture = R.ladders.find(_.family == "pca").get
    val pca = PcaVarianceProtocol()
    val pcaOriginal = accepted(pca.observed(accepted(PcaVarianceState.from(matrix(pcaFixture.x)))))
    val pcaNegated = accepted(pca.observed(accepted(PcaVarianceState.from(negate(matrix(pcaFixture.x))))))
    val pcaReversed = accepted(pca.observed(accepted(PcaVarianceState.from(reverseColumns(matrix(pcaFixture.x))))))
    assertEqualsDouble(pcaNegated, pcaOriginal, 1e-12)
    assertEqualsDouble(pcaReversed, pcaOriginal, 1e-12)

    val plscFixture = R.ladders.find(_.family == "plsc").get
    val plsc = PlscCovarianceProtocol()
    val x = matrix(plscFixture.x)
    val y = matrix(plscFixture.y.get)
    val original = accepted(plsc.observed(accepted(PlscCovarianceState.from(x, y))))
    val signChanged = accepted(plsc.observed(accepted(PlscCovarianceState.from(x, negate(y)))))
    val coordinatesChanged = accepted(plsc.observed(accepted(
      PlscCovarianceState.from(reverseColumns(x), reverseColumns(y))
    )))
    assertEqualsDouble(signChanged, original, 1e-8)
    assertEqualsDouble(coordinatesChanged, original, 1e-8)
  }

  test("a declared repeated-root subspace is invariant to rotation within that subspace") {
    val original = InferenceNumerics.matrixFromRows(Vector(
      Vector(1.0, 0.0, 0.0),
      Vector(-1.0, 0.0, 0.0),
      Vector(0.0, 1.0, 0.0),
      Vector(0.0, -1.0, 0.0),
      Vector(0.0, 0.0, 0.5),
      Vector(0.0, 0.0, -0.5)
    ))
    val c = Math.cos(Math.PI / 5.0)
    val s = Math.sin(Math.PI / 5.0)
    val rotation = InferenceNumerics.matrixFromRows(Vector(
      Vector(c, -s, 0.0),
      Vector(s, c, 0.0),
      Vector(0.0, 0.0, 1.0)
    ))
    val rotated = InferenceNumerics.multiply(original, rotation)
    val protocol = PcaVarianceProtocol()
    val originalRoots = accepted(protocol.roots(accepted(PcaVarianceState.from(original))))
    val rotatedRoots = accepted(protocol.roots(accepted(PcaVarianceState.from(rotated))))
    originalRoots.zip(rotatedRoots).foreach { case (left, right) =>
      assertEqualsDouble(left, right, 1e-12)
    }

    val plane = accepted(ComponentSet.from(Vector(accepted(ComponentIx(0)), accepted(ComponentIx(1)))))
    val axis = ComponentSet.one(accepted(ComponentIx(2)))
    val policy = UnitPolicy.Declared(accepted(DeclaredUnitGroups.from(Vector(plane, axis))))
    val originalUnits = accepted(LatentUnitFormation.form(originalRoots, Vector(true, true, false), policy))
    val rotatedUnits = accepted(LatentUnitFormation.form(rotatedRoots, Vector(true, true, false), policy))
    assertEquals(originalUnits.map(_.unit), rotatedUnits.map(_.unit))
    assertEquals(originalUnits.map(_.selected), rotatedUnits.map(_.selected))
    originalUnits.zip(rotatedUnits).foreach { case (left, right) =>
      left.roots.zip(right.roots).foreach { case (leftRoot, rightRoot) =>
        assertEqualsDouble(leftRoot, rightRoot, 1e-12)
      }
    }
    assertEquals(originalUnits.head.unit.identifiability, Identifiability.UnorientedSubspace)
  }

  test("PCA and PLSC execute through one deterministic ladder interpreter") {
    val pcaFixture = R.ladders.find(_.family == "pca").get
    val pcaState = accepted(PcaVarianceState.from(matrix(pcaFixture.x)))
    InferenceExecutor.runLadder(pcaState, PcaVarianceProtocol(), config) match
      case Left(InferenceError.UnsupportedUnitPolicy(UnitPolicy.GroupNearTies(_))) => ()
      case other => fail(s"expected grouped-unit refusal, got $other")

    val axisConfig = config.copy(unitPolicy = UnitPolicy.SingleAxes)
    val pcaFirst = accepted(InferenceExecutor.runLadder(pcaState, PcaVarianceProtocol(), axisConfig))
    val pcaSecond = accepted(InferenceExecutor.runLadder(pcaState, PcaVarianceProtocol(), axisConfig))
    assertEquals(pcaFirst.ladder, pcaSecond.ladder)
    assertEquals(pcaFirst.ladder.rejectedThrough, 2)
    assert(pcaFirst.units.forall(_.unit.identifiability == Identifiability.OrientableAxis))
    assert(pcaFirst.units.take(2).forall(_.selected))
    assert(pcaFirst.units.drop(2).forall(unit => !unit.selected))

    val plscFixture = R.ladders.find(_.family == "plsc").get
    val plscState = accepted(PlscCovarianceState.from(
      matrix(plscFixture.x),
      matrix(plscFixture.y.get)
    ))
    val plscConfig = axisConfig.copy(seed = Seed.fromLong(plscFixture.seed))
    val plscFirst = accepted(InferenceExecutor.runLadder(plscState, PlscCovarianceProtocol(), plscConfig))
    val plscSecond = accepted(InferenceExecutor.runLadder(plscState, PlscCovarianceProtocol(), plscConfig))
    assertEquals(plscFirst.ladder, plscSecond.ladder)
    assertEquals(plscFirst.ladder.rejectedThrough, 2)
  }

  test("compiled programs make their design, seed, and sequential policy authoritative") {
    val rows = accepted(RowCount(6))
    val partition = accepted(RowPartition.from(
      rows,
      Vector(Vector(0, 1, 2), Vector(3, 4, 5))
    ))
    val draws = accepted(MonteCarloDraws(5))
    val alpha = accepted(Alpha(0.5))
    val batch = accepted(BatchSize(5))
    val programSeed: Seed = Seed.fromLong(44L)
    val program = InferenceCompiler.compile(InferenceSpec.of(
      FitDescriptor.Pca,
      TargetSpec.VarianceRoots,
      NullSpec.PermuteRows,
      ResamplingDesign.withinBlocks(partition),
      UnitPolicy.SingleAxes,
      MonteCarloPolicy.Sequential(
        draws,
        alpha,
        batch,
        SequentialBoundary.FromAlpha
      ),
      RequestedEvidence.SignificanceOnly,
      programSeed
    ))
    var protocolTarget: TargetSpec[TargetKind.VarianceRoots] =
      TargetSpec.VarianceRoots
    val protocol = new ExactLadderProtocol[
      Unit,
      PcaFitFamily,
      TargetKind.VarianceRoots,
      NullKind.RowPermutation
    ]:
      type Action = resample4s.kernel.Permutation

      override val fit: FitDescriptor[PcaFitFamily] = FitDescriptor.Pca
      override def target: TargetSpec[TargetKind.VarianceRoots] =
        protocolTarget
      override val nullHypothesis: NullSpec[NullKind.RowPermutation] =
        NullSpec.PermuteRows
      override val validity: ValidityClaim = ValidityClaim.Exact
      override def roots(initial: Unit): Either[InferenceError, Vector[Double]] =
        Right(Vector(1.0))
      override def observed(state: Unit): Either[InferenceError, Double] =
        Right(1.0)
      override def action(
          state: Unit,
          replicate: ReplicateId,
          seed: Seed
      ): Either[InferenceError, Action] =
        Left(InferenceError.UnsupportedProblem("legacy unrestricted action was used"))
      override def actionForDesign(
          state: Unit,
          design: ResamplingDesign[?],
          replicate: ReplicateId,
          seed: Seed
      ): Either[InferenceError, Action] =
        PermutationAction.forDesign(design).flatMap(_.draw(seed, replicate))
      override def validateDesign(
          state: Unit,
          design: ResamplingDesign[?]
      ): Either[InferenceError, Unit] =
        PermutationAction.forDesign(design).map(_ => ())
      override def nullStatistic(
          state: Unit,
          step: ComponentIx,
          replicate: ReplicateId,
          action: Action
      ): Either[InferenceError, Double] =
        val rows = action.toIArray
        val preservesBlocks =
          rows.indices.forall(index => index / 3 == rows(index) / 3)
        if preservesBlocks then Right(0.0)
        else Left(InferenceError.UnsupportedProblem("action crossed a compiled block"))
      override def remove(state: Unit): Either[InferenceError, Unit] =
        Right(())

    val run = accepted(InferenceExecutor.runProgram(
      (),
      program,
      protocol,
      LadderExecutionConfig(
        draws,
        draws,
        alpha,
        batch,
        accepted(LadderSteps(1)),
        Seed.fromLong(999L)
      )
    ))

    val receipt = run.ladder.steps.head.evidence match
      case Evidence.Computed(value) => value
      case other                    => fail(s"expected computed evidence, obtained $other")
    receipt.provenance match
      case ReplicateProvenance.Deterministic(seed, _, _) =>
        assertEquals(seed.value, programSeed.value)
      case other => fail(s"expected deterministic provenance, obtained $other")
    assertEquals(run.validity, program.validity)

    protocolTarget = new TargetSpec[TargetKind.VarianceRoots]:
      override val label: TargetLabel = accepted(TargetLabel("custom-variance-roots"))
      override val alternative: Alternative = Alternative.Greater
      override val invariance: TargetInvariance =
        TargetInvariance.OrthogonalSubspaceRotation

    InferenceExecutor.runProgram(
      (),
      program,
      protocol,
      LadderExecutionConfig(
        draws,
        draws,
        alpha,
        batch,
        accepted(LadderSteps(1)),
        Seed.fromLong(999L)
      )
    ) match
      case Left(InferenceError.TargetProtocolMismatch(compiled, actual)) =>
        assertEquals(compiled, TargetSpec.VarianceRoots.label)
        assertEquals(actual, protocolTarget.label)
      case other =>
        fail(s"expected target-protocol mismatch, obtained $other")
  }

  test("the significance ladder refuses unsupported evidence requests before execution") {
    val fixture = R.ladders.find(_.family == "pca").get
    val state = accepted(PcaVarianceState.from(matrix(fixture.x)))
    val rows = accepted(RowCount(fixture.x.rows))
    val draws = accepted(MonteCarloDraws(5))
    val alpha = accepted(Alpha(0.5))
    val batch = accepted(BatchSize(5))
    val config = LadderExecutionConfig(
      draws,
      draws,
      alpha,
      batch,
      accepted(LadderSteps(1)),
      Seed.fromLong(44L)
    )

    Vector(
      RequestedEvidence.StabilityOnly,
      RequestedEvidence.SignificanceAndStability
    ).foreach { requested =>
      val program = InferenceCompiler.compile(InferenceSpec.of(
        FitDescriptor.Pca,
        TargetSpec.VarianceRoots,
        NullSpec.PermuteRows,
        ResamplingDesign.exchangeableRows(rows),
        UnitPolicy.SingleAxes,
        MonteCarloPolicy.Sequential(
          draws,
          alpha,
          batch,
          SequentialBoundary.FromAlpha
        ),
        requested,
        Seed.fromLong(44L)
      ))

      InferenceExecutor.runProgram(
        state,
        program,
        PcaVarianceProtocol(),
        config
      ) match
        case Left(InferenceError.UnsupportedEvidence(actual)) =>
          assertEquals(actual, requested)
        case other =>
          fail(s"expected unsupported-evidence refusal, obtained $other")
    }
  }

  test("compiled conditioned designs fail before the first replicate without prepared resources") {
    val fixture = R.ladders.find(_.family == "pca").get
    val state = accepted(PcaVarianceState.from(matrix(fixture.x)))
    val rows = accepted(RowCount(fixture.x.rows))
    val draws = accepted(MonteCarloDraws(5))
    val alpha = accepted(Alpha(0.5))
    val batch = accepted(BatchSize(5))
    val reference = accepted(ConditioningRef("nuisance-v1"))
    val program = InferenceCompiler.compile(InferenceSpec.of(
      FitDescriptor.Pca,
      TargetSpec.VarianceRoots,
      NullSpec.PermuteRows,
      ResamplingDesign.nuisanceAdjusted(
        rows,
        reference,
        WhiteningRequirement.Required
      ),
      UnitPolicy.SingleAxes,
      MonteCarloPolicy.Sequential(
        draws,
        alpha,
        batch,
        SequentialBoundary.FromAlpha
      ),
      RequestedEvidence.SignificanceOnly,
      Seed.fromLong(44L)
    ))
    val result = InferenceExecutor.runProgram(
      state,
      program,
      PcaVarianceProtocol(),
      LadderExecutionConfig(
        draws,
        draws,
        alpha,
        batch,
        accepted(LadderSteps(1)),
        Seed.fromLong(44L)
      )
    )

    result match
      case Left(InferenceError.UnsupportedProblem(message)) =>
        assert(message.contains("projection and whitening resources"))
      case other =>
        fail(s"expected conditioned-design preflight refusal, obtained $other")
  }

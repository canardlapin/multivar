package multivar
package family.paired

import multivar.analysis.*
import multivar.core.*

import gale.linalg.Matrix

class TaskDecompositionSuite extends munit.FunSuite:

  private val tol = 1e-8

  private val x = Matrix(6, 3)(
    1.0, 0.0, 2.0,
    2.0, 1.0, 1.0,
    3.0, 1.0, 0.0,
    4.0, 2.0, 1.0,
    5.0, 3.0, 2.0,
    6.0, 4.0, 3.0
  )

  private val heldOut = Matrix(2, 3)(
    2.5, 1.5, 1.0,
    3.5, 2.0, 2.5
  )

  test("public decompose rejects Whole and points callers to Gpca"):
    val result = TaskComponents.decompose(x, TaskResolution.Whole, components = 1)
    assert(result.isLeft)
    result.left.foreach:
      case MultivarError.UnsupportedEstimator(detail) =>
        assert(detail.contains("scientific resolution"))
        assert(detail.contains("Gpca.fit"))
      case other =>
        fail(s"expected UnsupportedEstimator, got $other")

  test("resolved decompose reports observation/feature axes and feature-specified transfer"):
    val featureGroups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents
        .decompose(x, TaskResolution.Target.weights(featureGroups), components = 1)
        .toOption
        .get
    assertEquals(fit.sourceAxis, "observations")
    assertEquals(fit.targetAxis, "features")
    assertEquals(fit.focus.rank, 1)
    assertEquals(fit.remainder.rank, 0)
    assertEquals(fit.resolutionTransfer, TaskResolutionTransfer.FeatureSpecified)

  test("available gain splits across focus and remainder for Target.weights"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents
        .decompose(
          x,
          RankBudget.Split(focus = 1, remainder = 1),
          TaskResolution.Target.weights(groups)
        )
        .toOption
        .get
    assertEqualsDouble(
      fit.value.availableGain,
      fit.focus.value.availableGain + fit.remainder.value.availableGain,
      tol
    )
    assertEqualsDouble(
      fit.value.fittedLoss + fit.value.retainedGain,
      fit.value.baselineLoss,
      tol
    )

  test("new-data reconstruct equals baseline plus focus and remainder contributions"):
    val groups = Matrix(3, 2)(
      1.0, 0.0,
      0.0, 1.0,
      0.0, 0.0
    )
    val fit =
      TaskComponents
        .decompose(
          x,
          RankBudget.Split(1, 1),
          TaskResolution.Target.structure(groups),
          TaskComponents.DecompositionOptions()
        )
        .toOption
        .get
    val reconstructed = fit.reconstruct(heldOut).toOption.get
    val base = fit.baseline(heldOut.rows).toOption.get
    val focus = fit.focusContribution(heldOut).toOption.get
    val rem = fit.remainderContribution(heldOut).toOption.get
    assertMatrixClose(reconstructed, TaskMath.add(base, TaskMath.add(focus, rem)), tol)

  test("training transform reproduces stored canonical scores"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents
        .decompose(
          x,
          RankBudget.Split(1, 1),
          TaskResolution.Target.weights(groups),
          TaskComponents.DecompositionOptions()
        )
        .toOption
        .get
    val scores = fit.transform(x).toOption.get
    val expected = hstack(fit.focus.sourceScores, fit.remainder.sourceScores)
    assertMatrixClose(scores, expected, tol)

  test("training reconstruct matches reconstructTraining and contribution assembly"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents.decompose(x, TaskResolution.Target.weights(groups), 1).toOption.get
    val fromNewApi = fit.reconstruct(x).toOption.get
    val stored = fit.reconstructTraining.toOption.get
    assertMatrixClose(fromNewApi, stored, tol)
    val assembled =
      TaskMath.add(
        fit.baseline(x.rows).toOption.get,
        TaskMath.add(
          fit.focusContribution(x).toOption.get,
          fit.remainderContribution(x).toOption.get
        )
      )
    assertMatrixClose(fromNewApi, assembled, tol)

  test("full-rank branch transfer: training reconstruct matches inductive reconstruct"):
    val groups = Matrix(3, 2)(
      1.0, 0.0,
      0.0, 1.0,
      0.0, 0.0
    )
    val fit =
      TaskComponents
        .decompose(
          x,
          RankBudget.Split(1, 1),
          TaskResolution.Target.weights(groups),
          TaskComponents.DecompositionOptions()
        )
        .toOption
        .get
    assertMatrixClose(fit.reconstruct(x).toOption.get, fit.reconstructTraining.toOption.get, tol)
    assertEquals(fit.focus.rank + fit.remainder.rank, fit.componentIndex.length)
    assert(fit.componentIndex.nonEmpty)

  test("ZF_b recovers Z_b for retained focus and remainder blocks"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val resolution = TaskResolution.Target.weights(groups)
    val options = TaskComponents.DecompositionOptions(
      centering = GpcaCentering.Ordinary
    )
    // Build the same Z/Z_b path the engine uses and check F = Z^+ Z_b.
    val prepared = PreprocessSpec.Center.fit(MatrixView.dense(x)).toOption.get
    val working = prepared.transform(MatrixView.dense(x)).toOption.get.toDense().toOption.get
    val task =
      CanonicalTask
        .fromTable(
          working,
          MetricSpec.identity(working.rows).toOption.get,
          MetricSpec.identity(working.cols).toOption.get,
          DenseSolvers.symmetricEigen,
          1e-12,
          StoragePolicy.AllowDense
        )
        .toOption
        .get
    val core =
      ResolvedTaskProblem(task, resolution).fit(RankBudget.Split(1, 1), DenseSolvers.svd, 1e-12).toOption.get
    val zFocus = task.whitenedReconstruction(core.focusSvd)
    val zRem = task.whitenedReconstruction(core.remainderSvd)
    val fFocus =
      FeatureTransfer.applyPseudoInverse(task.matrix, zFocus, DenseSolvers.svd, 1e-12).toOption.get
    val fRem =
      FeatureTransfer.applyPseudoInverse(task.matrix, zRem, DenseSolvers.svd, 1e-12).toOption.get
    assertMatrixClose(GaleNumerics.multiply(task.matrix, fFocus), zFocus, tol)
    assertMatrixClose(GaleNumerics.multiply(task.matrix, fRem), zRem, tol)
    // Approximate idempotence / complementarity on the feature support of Z.
    assertMatrixClose(GaleNumerics.multiply(fFocus, fFocus), fFocus, 1e-6)
    assertMatrixClose(GaleNumerics.multiply(fRem, fRem), fRem, 1e-6)
    val cross = GaleNumerics.multiply(fFocus, fRem)
    assertEqualsDouble(TaskMath.frobeniusNorm2(cross), 0.0, 1e-6)

  test("transfer idempotence: focus and remainder analysis are complementary on training Z"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents
        .decompose(
          x,
          RankBudget.Split(1, 1),
          TaskResolution.Target.weights(groups),
          TaskComponents.DecompositionOptions()
        )
        .toOption
        .get
    val focusScores = fit.focusContribution(x).toOption.get
    val remScores = fit.remainderContribution(x).toOption.get
    // Branch contributions should not reintroduce each other's retained mass:
    // their working sum equals the retained reconstruction contribution.
    val sum = TaskMath.add(focusScores, remScores)
    val retained =
      TaskMath.add(
        fit.focusContributionTraining.toOption.get,
        fit.remainderContributionTraining.toOption.get
      )
    assertMatrixClose(sum, retained, tol)

  test("Target.weights new-data scores are equivariant under nonsingular reparameterization"):
    val h = Matrix(3, 1)(1.0, 2.0, 0.0)
    val ha = Matrix(3, 1)(2.0, 4.0, 0.0)
    val left =
      TaskComponents.decompose(x, TaskResolution.Target.weights(h), 1).toOption.get
    val right =
      TaskComponents.decompose(x, TaskResolution.Target.weights(ha), 1).toOption.get
    assertMatrixCloseAbsColumns(left.transform(heldOut).toOption.get, right.transform(heldOut).toOption.get, tol)
    assertMatrixClose(left.reconstruct(heldOut).toOption.get, right.reconstruct(heldOut).toOption.get, tol)

  test("Source resolution marks TrainingSourceInduced and still projects new rows"):
    val observationGroups = Matrix(6, 1)(1.0, 1.0, 0.0, 0.0, 0.0, 0.0)
    val fit =
      TaskComponents
        .decompose(x, TaskResolution.Source.covariance(observationGroups), 1)
        .toOption
        .get
    assertEquals(fit.resolutionTransfer, TaskResolutionTransfer.TrainingSourceInduced)
    val scores = fit.transform(heldOut).toOption.get
    assertEquals(scores.rows, heldOut.rows)
    assertEquals(scores.cols, 1)
    val hat = fit.reconstruct(heldOut).toOption.get
    assertEquals(hat.rows, heldOut.rows)
    assertEquals(hat.cols, heldOut.cols)

  test("branch transform matches top-level column slices for Focus/Split"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents
        .decompose(
          x,
          RankBudget.Split(1, 1),
          TaskResolution.Target.weights(groups),
          TaskComponents.DecompositionOptions()
        )
        .toOption
        .get
    val top = fit.transform(heldOut).toOption.get
    val focus = fit.focus.transform(heldOut).toOption.get
    val rem = fit.remainder.transform(heldOut).toOption.get
    assertMatrixClose(focus, GaleNumerics.selectColumns(top, Vector(0)), tol)
    if rem.cols > 0 then
      assertMatrixClose(rem, GaleNumerics.selectColumns(top, Vector(1)), tol)
    assertMatrixClose(
      fit.focus.reconstructContribution(heldOut).toOption.get,
      fit.focusContribution(heldOut).toOption.get,
      tol
    )

  test("validationLoss is held-out Frobenius evaluation distinct from TaskValue"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents.decompose(x, TaskResolution.Target.weights(groups), 1).toOption.get
    val report =
      fit
        .validationLoss(heldOut, DecompositionValidationLoss.ReconstructionFrobenius)
        .toOption
        .get
    val resid = fit.residual(heldOut).toOption.get
    assertEqualsDouble(report.residualFrobenius, TaskMath.frobeniusNorm(resid), tol)
    assertEqualsDouble(report.fittedLoss, TaskMath.frobeniusNorm2(resid), tol)
    assert(report.baselineLoss >= report.fittedLoss - tol)
    // Training gains remain a different vocabulary.
    assert(math.abs(report.fittedLoss - fit.value.fittedLoss) > tol || fit.value.fittedLoss == 0.0)

  test("residual is observation minus reconstruction"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents.decompose(x, TaskResolution.Target.weights(groups), 1).toOption.get
    val resid = fit.residual(heldOut).toOption.get
    val expected = MatrixOps.subtract(heldOut, fit.reconstruct(heldOut).toOption.get)
    assertMatrixClose(resid, expected, tol)

  test("MatrixView diagnostics honor the caller's materialization policy"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents.decompose(x, TaskResolution.Target.weights(groups), 1).toOption.get
    val sparse =
      SparseMatrixView
        .fromRows(
          Vector(
            Vector(2.5, 1.5, 1.0),
            Vector(3.5, 2.0, 2.5)
          )
        )
        .toOption
        .get

    assert(fit.residual(sparse, StoragePolicy.PreserveSparse).isLeft)
    assert(fit.resolve(sparse, StoragePolicy.PreserveSparse).isLeft)
    assert(fit.unsupportedFraction(sparse, StoragePolicy.PreserveSparse).isLeft)
    assert(
      fit
        .validationLoss(
          sparse,
          DecompositionValidationLoss.ReconstructionFrobenius,
          StoragePolicy.PreserveSparse
        )
        .isLeft
    )

  test("resolve parts reassemble observations and isolate unsupported mass"):
    val groups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val fit =
      TaskComponents.decompose(x, TaskResolution.Target.weights(groups), 1).toOption.get
    val resolved = fit.resolve(heldOut).toOption.get
    val assembled =
      TaskMath.add(
        fit.baseline(heldOut.rows).toOption.get,
        TaskMath.add(
          resolved.focus,
          TaskMath.add(resolved.supportedRemainder, resolved.unsupported)
        )
      )
    assertMatrixClose(assembled, heldOut, tol)
    // Training rows lie in the task support: unsupported fraction ≈ 0.
    assertEqualsDouble(fit.unsupportedFraction(x).toOption.get, 0.0, 1e-6)
    // Held-out rows may carry unsupported mass; fraction is in [0, 1].
    val frac = fit.unsupportedFraction(heldOut).toOption.get
    assert(frac >= -tol && frac <= 1.0 + tol)

  test("resolve focus/remainder match full scientific transfers, not only retained rank"):
    val groups = Matrix(3, 2)(
      1.0, 0.0,
      0.0, 1.0,
      0.0, 0.0
    )
    val fit =
      TaskComponents
        .decompose(
          x,
          RankBudget.Focus(1),
          TaskResolution.Target.weights(groups),
          TaskComponents.DecompositionOptions()
        )
        .toOption
        .get
    val resolved = fit.resolve(x).toOption.get
    // Full scientific remainder can be nonzero even when RankBudget.Focus retains none of it.
    assert(TaskMath.frobeniusNorm2(resolved.supportedRemainder) > tol)
    assertEquals(fit.remainder.rank, 0)

  test("Target design shape errors name the feature axis"):
    val observationDesign = Matrix(6, 1)(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    val result =
      TaskComponents.decompose(x, TaskResolution.Target.weights(observationDesign), components = 1)
    assert(result.isLeft)
    result.left.foreach:
      case MultivarError.MatrixShapeMismatch(detail) =>
        assert(detail.contains("Target resolution"))
        assert(detail.contains("3 features"))
        assert(detail.contains("6 rows"))
      case other =>
        fail(s"expected MatrixShapeMismatch, got $other")

  test("Source design shape errors name the observation axis"):
    val featureDesign = Matrix(3, 1)(1.0, 0.0, 0.0)
    val result =
      TaskComponents.decompose(x, TaskResolution.Source.covariance(featureDesign), components = 1)
    assert(result.isLeft)
    result.left.foreach:
      case MultivarError.MatrixShapeMismatch(detail) =>
        assert(detail.contains("Source.covariance"))
        assert(detail.contains("6 observations"))
        assert(detail.contains("3 rows"))
      case other =>
        fail(s"expected MatrixShapeMismatch, got $other")

  test("identity-metric Takane A–D resolutions partition the same available gain"):
    val featureGroups = Matrix(3, 1)(1.0, 0.0, 0.0)
    val observationGroups = Matrix(6, 1)(1.0, 1.0, 0.0, 0.0, 0.0, 0.0)
    val weight =
      TaskComponents.decompose(x, TaskResolution.Target.weights(featureGroups), 1).toOption.get
    val structure =
      TaskComponents
        .decompose(x, TaskResolution.Target.structure(featureGroups), 1)
        .toOption
        .get
    val covariance =
      TaskComponents
        .decompose(x, TaskResolution.Source.covariance(observationGroups), 1)
        .toOption
        .get
    val regression =
      TaskComponents
        .decompose(x, TaskResolution.Source.regression(observationGroups), 1)
        .toOption
        .get
    assertEqualsDouble(weight.value.availableGain, structure.value.availableGain, tol)
    assertEqualsDouble(weight.value.availableGain, covariance.value.availableGain, tol)
    assertEqualsDouble(weight.value.availableGain, regression.value.availableGain, tol)
    assertEquals(weight.resolutionTransfer, TaskResolutionTransfer.FeatureSpecified)
    assertEquals(covariance.resolutionTransfer, TaskResolutionTransfer.TrainingSourceInduced)
    assert(regression.focus.value.availableGain <= covariance.focus.value.availableGain + tol)

  test("Target.weights design equivariance under nonsingular reparameterization"):
    val h = Matrix(3, 1)(1.0, 2.0, 0.0)
    val ha = Matrix(3, 1)(2.0, 4.0, 0.0)
    val left =
      TaskComponents.decompose(x, TaskResolution.Target.weights(h), 1).toOption.get
    val right =
      TaskComponents.decompose(x, TaskResolution.Target.weights(ha), 1).toOption.get
    assertEqualsDouble(left.focus.value.availableGain, right.focus.value.availableGain, tol)
    assertEqualsDouble(left.remainder.value.availableGain, right.remainder.value.availableGain, tol)

  test("internal Whole Focus recovers PCA singular values and transform under identity metrics"):
    val pca = Pca.fit(x, components = 2).toOption.get
    val gmd =
      TaskDecompositionEngine
        .decomposeAllowingWhole(
          MatrixView.dense(x),
          RankBudget.Focus(2),
          TaskComponents.DecompositionOptions()
        )
        .toOption
        .get
    assertEquals(gmd.focus.singularValues.length, pca.singularValues.length)
    var i = 0
    while i < gmd.focus.singularValues.length do
      assertEqualsDouble(gmd.focus.singularValues(i), pca.singularValues(i), tol)
      i += 1
    assertEquals(gmd.remainder.rank, 0)
    assertMatrixCloseAbsColumns(gmd.transform(x).toOption.get, pca.transform(x).toOption.get, tol)

  private def hstack(left: gale.linalg.DMat, right: gale.linalg.DMat): gale.linalg.DMat =
    if right.cols == 0 then left
    else if left.cols == 0 then right
    else
      val out = Matrix.newBuilder(left.rows, left.cols + right.cols)
      var row = 0
      while row < left.rows do
        var col = 0
        while col < left.cols do
          out(row, col) = left(row, col)
          col += 1
        col = 0
        while col < right.cols do
          out(row, left.cols + col) = right(row, col)
          col += 1
        row += 1
      out.result()

  private def assertMatrixClose(left: gale.linalg.DMat, right: gale.linalg.DMat, tolerance: Double): Unit =
    assertEquals(left.rows, right.rows)
    assertEquals(left.cols, right.cols)
    var row = 0
    while row < left.rows do
      var col = 0
      while col < left.cols do
        assertEqualsDouble(left(row, col), right(row, col), tolerance)
        col += 1
      row += 1

  /** Allow column sign flips (SVD gauge). */
  private def assertMatrixCloseAbsColumns(
      left: gale.linalg.DMat,
      right: gale.linalg.DMat,
      tolerance: Double
  ): Unit =
    assertEquals(left.rows, right.rows)
    assertEquals(left.cols, right.cols)
    var col = 0
    while col < left.cols do
      var flip = 1.0
      var row = 0
      while row < left.rows && math.abs(left(row, col)) <= tolerance && math.abs(right(row, col)) <= tolerance do
        row += 1
      if row < left.rows && left(row, col) * right(row, col) < 0.0 then flip = -1.0
      row = 0
      while row < left.rows do
        assertEqualsDouble(left(row, col), flip * right(row, col), tolerance)
        row += 1
      col += 1

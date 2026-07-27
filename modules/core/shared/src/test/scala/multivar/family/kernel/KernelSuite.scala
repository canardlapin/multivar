package multivar
package family.kernel

import multivar.core.*


import gale.linalg.DMat

class KernelSuite extends munit.FunSuite:

  private def ref(id: String, role: SpaceRole, dimension: Int): SpaceRef =
    SpaceRef.of(id, role, dimension).toOption.get

  private def typedInput[Rows <: SemanticSpace, Features <: SemanticSpace](
      values: DMat,
      rows: SpaceEvidence[Rows],
      features: SpaceEvidence[Features],
      id: String
  ): KernelInput[Rows, Features] =
    KernelInput.from(
      MatrixView.dense(values),
      rows,
      features,
      ValueIdentity.source(ValueId.unsafe(id))
    ).toOption.get

  private def assertMatrixClose(actual: DMat, expected: DMat, tol: Double): Unit =
    assertEquals(actual.rows, expected.rows)
    assertEquals(actual.cols, expected.cols)
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assertEqualsDouble(actual(row, col), expected(row, col), tol)
        col += 1
      row += 1

  private def assertMatrixClose(actual: DMat, expected: Vector[Vector[Double]], tol: Double): Unit =
    assertEquals(actual.rows, expected.length)
    assertEquals(actual.cols, expected.headOption.map(_.length).getOrElse(0))
    var row = 0
    while row < actual.rows do
      var col = 0
      while col < actual.cols do
        assertEqualsDouble(actual(row, col), expected(row)(col), tol)
        col += 1
      row += 1

  private def tcross(matrix: DMat): DMat =
    GaleNumerics.multiply(matrix, matrix.transpose)

  private def diag(values: gale.linalg.DVec): DMat =
    MatrixOps.diagonal(values)

  private def explicitNystromKernel(input: DMat, landmarks: Vector[Int]): DMat =
    val landmarkData = RowGeometryOps.selectRows(input, landmarks)
    val c = GaleNumerics.multiply(input, landmarkData.transpose)
    val w = RowGeometryOps.selectRows(c, landmarks)
    val eigen = DenseSolvers.symmetricEigen.decompose(w).toOption.get
    val keep = landmarks.length
    val u = MatrixOps.takeColumns(eigen.vectors, keep)
    val lambda = MatrixOps.takeVector(eigen.values, keep)
    val inv = new Array[Double](keep)
    var i = 0
    while i < keep do
      inv(i) = 1.0 / lambda(i)
      i += 1
    val middle = GaleNumerics.multiply(GaleNumerics.multiply(u, MatrixOps.diagonal(GaleNumerics.vectorFromArray(inv))), u.transpose)
    GaleNumerics.multiply(GaleNumerics.multiply(c, middle), c.transpose)

  test("standard Nyström all-landmark linear fit matches exact kernel eigensystem") {
    val input = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(-1.0, 0.0),
          Vector(0.0, 1.0),
          Vector(1.0, 0.0),
          Vector(0.0, -1.0)
        )
      )
    )

    val fit = Nystrom.fit(
      input,
      ComponentCount(2).toOption.get,
      landmarks = Vector(0, 1, 2, 3),
      preproc = PreprocessSpec.Pass
    ).toOption.get
    val kernel = input.toDense().toOption.get
    val k = GaleNumerics.multiply(kernel, kernel.transpose)
    val residual = RowGeometryOps.subtract(
      GaleNumerics.multiply(k, NystromFit.eigenOf(fit).eigenvectors),
      GaleNumerics.multiply(NystromFit.eigenOf(fit).eigenvectors, diag(NystromFit.eigenOf(fit).eigenvalues))
    )

    assertEquals(NystromFit.methodOf(fit), NystromMethod.Standard)
    assertEquals(NystromFit.eigenOf(fit).components, 2)
    assertEqualsDouble(NystromFit.eigenOf(fit).eigenvalues(0), 2.0, 1e-9)
    assertEqualsDouble(NystromFit.eigenOf(fit).eigenvalues(1), 2.0, 1e-9)
    assertMatrixClose(residual, DMat.zeros(4, 2), 1e-8)
    assertMatrixClose(fit.transform(input).toOption.get, NystromFit.eigenOf(fit).scores, 1e-9)
  }

  test("dense Nyström conveniences preserve the checked fit and validate raw parameters") {
    val input = GaleNumerics.matrixFromRows(
      Vector(
        Vector(-1.0, 0.0),
        Vector(0.0, 1.0),
        Vector(1.0, 0.0),
        Vector(0.0, -1.0)
      )
    )
    val landmarks = Vector(0, 1, 2, 3)
    val convenient = Nystrom.fit(input, components = 2, landmarks).toOption.get
    val canonical =
      Nystrom.fit(MatrixView.dense(input), ComponentCount.unsafe(2), landmarks).toOption.get

    assertMatrixClose(convenient.scores, canonical.scores, 0.0)
    assertMatrixClose(convenient.transform(input).toOption.get, convenient.scores, 0.0)
    assertEquals(convenient.eigenvalues.copyData.toVector, canonical.eigenvalues.copyData.toVector)
    assert(Nystrom.fit(input, components = 0, landmarks).isLeft)
    assert(Nystrom.fitRbf(input, components = 2, landmarks, gamma = -0.1).isLeft)
  }

  test("partial-landmark standard Nyström reconstructs the explicit kernel approximation") {
    val x = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 0.0, 1.0),
        Vector(0.0, 1.0, 1.0),
        Vector(1.0, 1.0, 0.0),
        Vector(2.0, 0.0, 1.0),
        Vector(0.0, 2.0, 1.0)
      )
    )
    val landmarks = Vector(0, 2, 4)
    val fit = Nystrom.fit(
      MatrixView.dense(x),
      ComponentCount(3).toOption.get,
      landmarks = landmarks,
      preproc = PreprocessSpec.Pass
    ).toOption.get
    val explicit = explicitNystromKernel(x, landmarks)

    assertMatrixClose(tcross(NystromFit.eigenOf(fit).scores), explicit, 1e-8)
    assertMatrixClose(fit.transform(MatrixView.dense(x)).toOption.get, NystromFit.eigenOf(fit).scores, 1e-8)
  }

  test("all-landmark standard Nyström matches exact RBF kernel eigenvalues") {
    val x = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(-1.0, 0.0),
          Vector(0.0, 1.0),
          Vector(1.0, 0.0),
          Vector(0.0, -1.0)
        )
      )
    )
    val kernel = RbfKernel(gamma = 0.3).toOption.get
    val fit = Nystrom.fit(
      x,
      ComponentCount(3).toOption.get,
      landmarks = Vector(0, 1, 2, 3),
      kernel = kernel,
      preproc = PreprocessSpec.Pass
    ).toOption.get
    val k = kernel.compute(x, x).toOption.get
    val exact = DenseSolvers.symmetricEigen.decompose(k).toOption.get
    val residual = RowGeometryOps.subtract(
      GaleNumerics.multiply(k, NystromFit.eigenOf(fit).eigenvectors),
      GaleNumerics.multiply(NystromFit.eigenOf(fit).eigenvectors, diag(NystromFit.eigenOf(fit).eigenvalues))
    )

    assertEqualsDouble(NystromFit.eigenOf(fit).eigenvalues(0), exact.values(0), 1e-9)
    assertEqualsDouble(NystromFit.eigenOf(fit).eigenvalues(1), exact.values(1), 1e-9)
    assertEqualsDouble(NystromFit.eigenOf(fit).eigenvalues(2), exact.values(2), 1e-9)
    assertMatrixClose(residual, DMat.zeros(4, 3), 1e-8)
  }

  test("double Nyström with full intermediate rank reconstructs the standard approximation") {
    val x = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 0.0, 1.0),
        Vector(0.0, 1.0, 1.0),
        Vector(1.0, 1.0, 0.0),
        Vector(2.0, 0.0, 1.0),
        Vector(0.0, 2.0, 1.0)
      )
    )
    val landmarks = Vector(0, 2, 4)
    val fit = Nystrom.fit(
      MatrixView.dense(x),
      ComponentCount(3).toOption.get,
      landmarks = landmarks,
      preproc = PreprocessSpec.Pass,
      method = NystromMethod.DoubleNystrom(ComponentCount(3).toOption.get)
    ).toOption.get
    val explicit = explicitNystromKernel(x, landmarks)

    assertEquals(NystromFit.methodOf(fit).label, "double")
    assertMatrixClose(tcross(NystromFit.eigenOf(fit).scores), explicit, 1e-8)
    assertMatrixClose(fit.transform(MatrixView.dense(x)).toOption.get, NystromFit.eigenOf(fit).scores, 1e-8)
  }

  test("double Nyström with truncated intermediate rank matches an explicit truncated reference") {
    val x = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 2.0, 0.5),
        Vector(-1.0, 1.0, 1.5),
        Vector(2.0, -1.0, 1.0),
        Vector(0.5, 0.5, -2.0),
        Vector(-1.5, 2.0, 0.5),
        Vector(1.0, 1.0, 1.0)
      )
    )
    val landmarks = Vector(0, 1, 3, 5)
    val intermediateRank = 2
    val components = 2
    val fit = Nystrom.fit(
      MatrixView.dense(x),
      ComponentCount(components).toOption.get,
      landmarks = landmarks,
      preproc = PreprocessSpec.Pass,
      method = NystromMethod.DoubleNystrom(ComponentCount(intermediateRank).toOption.get)
    ).toOption.get

    // Independent reference: explicit truncated first stage, then the second-stage
    // eigensystem of W'W, using dense linear algebra only.
    val landmarkData = RowGeometryOps.selectRows(x, landmarks)
    val c = GaleNumerics.multiply(x, landmarkData.transpose)
    val kMm = RowGeometryOps.selectRows(c, landmarks)
    val first = DenseSolvers.symmetricEigen.decompose(kMm).toOption.get
    val vSL = MatrixOps.takeColumns(first.vectors, intermediateRank)
    val lambdaL = MatrixOps.takeVector(first.values, intermediateRank)
    val invSqrtL = new Array[Double](intermediateRank)
    var i = 0
    while i < intermediateRank do
      invSqrtL(i) = 1.0 / Math.sqrt(lambdaL(i))
      i += 1
    val firstWeights = GaleNumerics.multiply(vSL, MatrixOps.diagonal(GaleNumerics.vectorFromArray(invSqrtL)))
    val w = GaleNumerics.multiply(c, firstWeights)
    val second = DenseSolvers.symmetricEigen.decompose(GaleNumerics.crossProduct(w)).toOption.get
    val lambdaK = MatrixOps.takeVector(second.values, components)
    val vK = MatrixOps.takeColumns(second.vectors, components)
    val sqrtInvK = new Array[Double](components)
    i = 0
    while i < components do
      sqrtInvK(i) = 1.0 / Math.sqrt(lambdaK(i))
      i += 1
    val eigenWeights = GaleNumerics.multiply(
      firstWeights,
      GaleNumerics.multiply(vK, MatrixOps.diagonal(GaleNumerics.vectorFromArray(sqrtInvK)))
    )
    val sqrtK = new Array[Double](components)
    i = 0
    while i < components do
      sqrtK(i) = Math.sqrt(lambdaK(i))
      i += 1
    val expectedScores = GaleNumerics.multiply(
      GaleNumerics.multiply(c, eigenWeights),
      MatrixOps.diagonal(GaleNumerics.vectorFromArray(sqrtK))
    )

    NystromFit.stateOf(fit) match
      case state: DoubleNystromState =>
        assertEquals(state.firstStageEigenvectors.cols, intermediateRank, "first stage must be genuinely truncated")
      case other =>
        fail(s"expected a double Nyström state, got $other")
    assertEquals(NystromFit.eigenOf(fit).components, components)
    i = 0
    while i < components do
      assertEqualsDouble(NystromFit.eigenOf(fit).eigenvalues(i), lambdaK(i), 1e-9)
      i += 1
    assertMatrixClose(tcross(NystromFit.eigenOf(fit).scores), tcross(expectedScores), 1e-8)
    assertMatrixClose(fit.transform(MatrixView.dense(x)).toOption.get, NystromFit.eigenOf(fit).scores, 1e-8)
  }

  test("Nyström rejects component and intermediate-rank requests beyond the landmark count") {
    val x = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0, 0.0), Vector(0.0, 1.0), Vector(1.0, 1.0))))

    Nystrom.fit(x, ComponentCount(3).toOption.get, landmarks = Vector(0, 1)) match
      case Left(MultivarError.InvalidComponentRequest(requested, limit)) =>
        assertEquals(requested, 3)
        assertEquals(limit, 2)
      case other =>
        fail(s"expected component request rejection, got $other")

    Nystrom.fit(
      x,
      ComponentCount(2).toOption.get,
      landmarks = Vector(0, 1),
      method = NystromMethod.DoubleNystrom(ComponentCount(3).toOption.get)
    ) match
      case Left(MultivarError.InvalidComponentRequest(requested, limit)) =>
        assertEquals(requested, 3)
        assertEquals(limit, 2)
      case other =>
        fail(s"expected intermediate rank rejection, got $other")
  }

  test("Nyström rejects out-of-bounds landmark indices") {
    val x = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0, 0.0), Vector(0.0, 1.0), Vector(1.0, 1.0))))

    Nystrom.fit(x, ComponentCount(1).toOption.get, landmarks = Vector(0, 7)) match
      case Left(MultivarError.IndexOutOfBounds(IndexAxis.Row, 7, 3)) => ()
      case other => fail(s"expected landmark bounds rejection, got $other")

    Nystrom.fit(x, ComponentCount(1).toOption.get, landmarks = Vector(-1)) match
      case Left(MultivarError.IndexOutOfBounds(IndexAxis.Row, -1, 3)) => ()
      case other => fail(s"expected negative landmark rejection, got $other")
  }

  test("linear kernel rejects mismatched feature counts") {
    val left = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0, 2.0))))
    val right = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0, 2.0, 3.0))))

    LinearKernel().compute(left, right) match
      case Left(MultivarError.MatrixShapeMismatch(detail)) =>
        assert(detail.contains("feature counts"), detail)
      case other =>
        fail(s"expected linear kernel feature mismatch, got $other")
  }

  test("RBF kernels require a finite positive gamma") {
    RbfKernel(0.0) match
      case Left(MultivarError.InvalidKernelParameter("gamma", 0.0, _)) =>
      case other => fail(s"expected typed RBF gamma rejection, got $other")
    assert(RbfKernel(-1.0).isLeft)
    assert(RbfKernel(Double.NaN).isLeft)
    assert(RbfKernel(Double.PositiveInfinity).isLeft)
  }

  test("material asymmetry in a user kernel is rejected by default and accepted when unchecked") {
    val asymmetric = new Kernel:
      override def spec: KernelSpec =
        KernelSpec("asymmetric-linear")

      override def compute(left: MatrixView, right: MatrixView): Either[MultivarError, DMat] =
        LinearKernel().compute(left, right).map { out =>
          if out.rows == out.cols && out.cols > 1 then
            val data = out.copyData
            data(1) += 1e-8
            GaleNumerics.matrixFromRowMajor(out.rows, out.cols, data)
          else out
        }

    val x = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(1.0, 0.0),
          Vector(0.0, 1.0),
          Vector(1.0, 1.0),
          Vector(2.0, 1.0)
        )
      )
    )

    Nystrom.fit(x, ComponentCount(2).toOption.get, landmarks = Vector(0, 2, 3), kernel = asymmetric) match
      case Left(MultivarError.NonSymmetricMatrix(_, _, _, _)) => ()
      case other => fail(s"expected default symmetry policy to reject material asymmetry, got $other")

    val unchecked = Nystrom.fit(
      x,
      ComponentCount(2).toOption.get,
      landmarks = Vector(0, 2, 3),
      kernel = asymmetric,
      symmetry = KernelSymmetryPolicy.SymmetrizeUnchecked
    )
    assert(unchecked.isRight, s"expected SymmetrizeUnchecked to accept the kernel, got $unchecked")
    assert(NystromFit.eigenOf(unchecked.toOption.get).eigenvalues(0).isFinite)

    val roundoff = new Kernel:
      override def spec: KernelSpec = KernelSpec("roundoff-asymmetric-linear")
      override def compute(left: MatrixView, right: MatrixView): Either[MultivarError, DMat] =
        LinearKernel().compute(left, right).map { out =>
          if out.rows == out.cols && out.cols > 1 then
            val data = out.copyData
            data(1) += 1e-14
            GaleNumerics.matrixFromRowMajor(out.rows, out.cols, data)
          else out
        }
    assert(
      Nystrom.fit(x, ComponentCount(2).toOption.get, landmarks = Vector(0, 2, 3), kernel = roundoff).isRight,
      "numerical noise within SymmetrizeWithin should still fit"
    )
  }

  test("Require keeps the unsymmetrized kernel while SymmetrizeWithin averages within tolerance") {
    val skewed = GaleNumerics.matrixFromRows(
      Vector(
        Vector(2.0, 1.0 + 1e-11),
        Vector(1.0, 3.0)
      )
    )
    val required = KernelSymmetryPolicy.Require(1e-10).prepare(skewed).toOption.get
    val symmetrized = KernelSymmetryPolicy.SymmetrizeWithin(1e-10).prepare(skewed).toOption.get
    assertEqualsDouble(required(0, 1), skewed(0, 1), 0.0)
    assertEqualsDouble(required(1, 0), skewed(1, 0), 0.0)
    assertEqualsDouble(symmetrized(0, 1), 0.5 * (skewed(0, 1) + skewed(1, 0)), 1e-18)
    assertEqualsDouble(symmetrized(1, 0), symmetrized(0, 1), 0.0)

    val large = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1e6, 1e6 + 1e-7),
        Vector(1e6, 2e6)
      )
    )
    // Relative tolerance: Require accepts magnitude-scaled skew that absolute SymmetrizeWithin rejects.
    assert(KernelSymmetryPolicy.Require(1e-12).prepare(large).isRight)
    KernelSymmetryPolicy.SymmetrizeWithin(1e-8).prepare(large) match
      case Left(MultivarError.NonSymmetricMatrix(0, 1, _, _)) => ()
      case other => fail(s"expected absolute-tolerance rejection, got $other")
  }

  test("centering diagnostics are derived from the fitted preprocessor") {
    val x = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(1.0, 0.0),
          Vector(0.0, 1.0),
          Vector(1.0, 1.0),
          Vector(2.0, 3.0)
        )
      )
    )

    val raw = Nystrom.fit(x, ComponentCount(1).toOption.get, landmarks = Vector(0, 2), preproc = PreprocessSpec.Pass).toOption.get
    assertEquals(NystromFit.centeringOf(raw), KernelCentering.Uncentered)
    assertEquals(raw.diagnostics.centering, KernelCentering.Uncentered)

    val centered = Nystrom.fit(x, ComponentCount(1).toOption.get, landmarks = Vector(0, 2), preproc = PreprocessSpec.Center).toOption.get
    assertEquals(NystromFit.centeringOf(centered), KernelCentering.InputPreprocessed)
    assertEquals(centered.diagnostics.centering, KernelCentering.InputPreprocessed)
  }

  test("out-of-sample projection follows the stored standard Nyström weights") {
    val x = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(1.0, 0.0),
          Vector(0.0, 1.0),
          Vector(1.0, 1.0),
          Vector(2.0, 1.0)
        )
      )
    )
    val newData = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0, 2.0), Vector(3.0, 0.0))))
    val fit = Nystrom.fit(x, ComponentCount(2).toOption.get, landmarks = Vector(0, 2)).toOption.get
    val kNew = LinearKernel().compute(newData, MatrixView.dense(NystromFit.landmarkDataOf(fit))).toOption.get
    val expected = GaleNumerics.multiply(kNew, NystromFit.stateOf(fit).scoreWeights)

    assertMatrixClose(fit.transform(newData).toOption.get, expected, 1e-10)
  }

  test("duplicate landmarks are rejected; canonicalize is explicit; rank-deficient kernels degrade") {
    val x = MatrixView.dense(GaleNumerics.matrixFromRows(Vector.fill(5)(Vector(1.0, 1.0))))
    Nystrom.fit(
      x,
      ComponentCount(3).toOption.get,
      landmarks = Vector(4, 0, 2, 0, 4),
      preproc = PreprocessSpec.Pass
    ) match
      case Left(MultivarError.DuplicateIndex(IndexAxis.Row, 0)) => ()
      case other => fail(s"expected duplicate landmark rejection, got $other")

    val canonical = LandmarkSet.canonicalize(Vector(4, 0, 2, 0, 4), rows = 5).toOption.get
    assertEquals(canonical.indices, Vector(0, 2, 4))

    val fit = Nystrom.fit(
      x,
      ComponentCount(3).toOption.get,
      landmarks = canonical.indices,
      preproc = PreprocessSpec.Pass
    ).toOption.get
    assertEquals(NystromFit.eigenOf(fit).components, 1)
    assert(NystromFit.eigenOf(fit).standardDeviations(0).isFinite)
  }

  private final class CountingKernel extends Kernel:
    var invocations = 0

    override def spec: KernelSpec =
      KernelSpec("counting-linear")

    override def compute(left: MatrixView, right: MatrixView): Either[MultivarError, DMat] =
      invocations += 1
      LinearKernel().compute(left, right)

  test("standard Nyström computes the landmark and all-landmark kernels exactly once each") {
    val x = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(1.0, 0.0, 1.0),
          Vector(0.0, 1.0, 1.0),
          Vector(1.0, 1.0, 0.0),
          Vector(2.0, 0.0, 1.0),
          Vector(0.0, 2.0, 1.0)
        )
      )
    )
    val kernel = new CountingKernel

    val fit = Nystrom.fit(x, ComponentCount(2).toOption.get, landmarks = Vector(0, 2, 4), kernel = kernel)

    assert(fit.isRight, s"expected counting-kernel fit to succeed, got $fit")
    assertEquals(kernel.invocations, 2)
  }

  test("double Nyström computes the n x m kernel once, not once per stage") {
    val x = MatrixView.dense(
      GaleNumerics.matrixFromRows(
        Vector(
          Vector(1.0, 0.0, 1.0),
          Vector(0.0, 1.0, 1.0),
          Vector(1.0, 1.0, 0.0),
          Vector(2.0, 0.0, 1.0),
          Vector(0.0, 2.0, 1.0)
        )
      )
    )
    val kernel = new CountingKernel

    val fit = Nystrom.fit(
      x,
      ComponentCount(2).toOption.get,
      landmarks = Vector(0, 2, 4),
      kernel = kernel,
      method = NystromMethod.DoubleNystrom(ComponentCount(2).toOption.get)
    )

    assert(fit.isRight, s"expected counting-kernel double fit to succeed, got $fit")
    assertEquals(kernel.invocations, 2)
  }

  test("non-finite kernel outputs fail before eigendecomposition") {
    val badKernel = new Kernel:
      override def spec: KernelSpec =
        KernelSpec("bad")

      override def compute(left: MatrixView, right: MatrixView): Either[MultivarError, DMat] =
        Right(GaleNumerics.matrixFromRows(Vector.fill(left.rows)(Vector.fill(right.rows)(Double.NaN))))

    val x = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0), Vector(2.0), Vector(3.0))))
    val result = Nystrom.fit(x, ComponentCount(1).toOption.get, landmarks = Vector(0, 1), kernel = badKernel)

    assert(result.swap.toOption.exists {
      case MultivarError.NonFiniteValue(_, _, _) => true
      case _                                     => false
    })
  }

  test("new-sample projection rejects wrong feature counts") {
    val x = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0, 0.0), Vector(0.0, 1.0), Vector(1.0, 1.0))))
    val fit = Nystrom.fit(x, ComponentCount(2).toOption.get, landmarks = Vector(0, 1)).toOption.get
    val bad = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0))))

    assert(fit.transform(bad).swap.toOption.exists(_.message.contains("expected 2 columns")))
  }

  test("typed Nyström artifacts preserve kernel roles, evidence, and low-rank storage") {
    val rows = ref("kernel.typed.training", SpaceRole.Samples, 4)
    val features = ref("kernel.typed.features", SpaceRole.Observed, 2)
    val values = GaleNumerics.matrixFromRows(
      Vector(
        Vector(1.0, 0.0),
        Vector(0.0, 1.0),
        Vector(1.0, 1.0),
        Vector(2.0, 1.0)
      )
    )
    val input = typedInput(values, rows.evidence, features.evidence, "kernel.typed.input")
    val fit = Nystrom.fitTyped(
      input,
      ComponentCount.unsafe(2),
      landmarks = Vector(0, 2, 3)
    ).toOption.get
    val operators = NystromFit.operatorOf(fit)

    assertEquals(operators.trainingRows.descriptor, rows.descriptor)
    assertEquals(operators.featureSpace.descriptor, features.descriptor)
    assertEquals(operators.landmarkKernel.role.value, OperatorRole.Kernel)
    assertEquals(operators.landmarkKernel.certificate.status, EvidenceStatus.Certified)
    assertEquals(operators.extensionKernel.role.value, OperatorRole.Kernel)
    assertEquals(operators.extensionKernel.certificate.status, EvidenceStatus.Unchecked)
    assertEquals(operators.approximateKernel.role.value, OperatorRole.Kernel)
    assertEquals(operators.approximateKernel.certificate.status, EvidenceStatus.Certified)
    assertEquals(operators.approximateKernel.representation, OperatorRepresentation.LowRank)
    assertEquals(operators.trainingScores.role.value, OperatorRole.Score)
    assertMatrixClose(operators.trainingScores.toDense.toOption.get, NystromFit.eigenOf(fit).scores, 1e-10)
    assertMatrixClose(operators.approximateKernel.toDense.toOption.get, tcross(NystromFit.eigenOf(fit).scores), 1e-10)
  }

  test("typed out-of-sample transforms enforce feature identity and retain row-space provenance") {
    val trainingRows = ref("kernel.transform.training", SpaceRole.Samples, 4)
    val features = ref("kernel.transform.features", SpaceRole.Observed, 2)
    val training = typedInput(
      GaleNumerics.matrixFromRows(
        Vector(Vector(1.0, 0.0), Vector(0.0, 1.0), Vector(1.0, 1.0), Vector(2.0, 1.0))
      ),
      trainingRows.evidence,
      features.evidence,
      "kernel.transform.input"
    )
    val fit = Nystrom.fitTyped(training, ComponentCount.unsafe(2), Vector(0, 2)).toOption.get
    val newRows = ref("kernel.transform.new-rows", SpaceRole.Samples, 2)
    val newValues = GaleNumerics.matrixFromRows(Vector(Vector(1.0, 2.0), Vector(3.0, 0.0)))
    val newInput = typedInput(
      newValues,
      newRows.evidence,
      features.evidence,
      "kernel.transform.new-input"
    )
    val transformed = NystromFit
      .operatorOf(fit)
      .transform(newInput, NystromFit.preprocessorOf(fit), NystromFit.kernelOf(fit))
      .toOption
      .get

    assertEquals(transformed.scores.codomain.descriptor.space, newRows.descriptor)
    assertEquals(
      transformed.scores.domain.descriptor.space,
      NystromFit.operatorOf(fit).componentSpace.descriptor
    )
    assertMatrixClose(
      transformed.values,
      fit.transform(MatrixView.dense(newValues)).toOption.get,
      1e-10
    )
    assert(transformed.extensionKernel.provenance.events.exists {
      case SemanticProvenanceEvent.Derived("kernel-extension", _) => true
      case _                                                       => false
    })

    val foreignFeatures = ref("kernel.transform.foreign-features", SpaceRole.Observed, 2)
    val foreignInput = typedInput(
      newValues,
      newRows.evidence,
      foreignFeatures.evidence,
      "kernel.transform.foreign-input"
    )
    assert(
      NystromFit
        .operatorOf(fit)
        .transform(foreignInput, NystromFit.preprocessorOf(fit), NystromFit.kernelOf(fit))
        .swap
        .toOption
        .exists(_.message.contains("does not match fitted space"))
    )
  }

  test("indefinite landmark kernels fail the PSD boundary and invalid tolerances fail early") {
    val indefinite = new Kernel:
      override def spec: KernelSpec = KernelSpec("indefinite")

      override def compute(left: MatrixView, right: MatrixView): Either[MultivarError, DMat] =
        if left.rows == 2 && right.rows == 2 then
          Right(GaleNumerics.matrixFromRows(Vector(Vector(1.0, 2.0), Vector(2.0, 1.0))))
        else LinearKernel().compute(left, right)

    val input = MatrixView.dense(GaleNumerics.matrixFromRows(Vector(Vector(1.0), Vector(2.0))))
    val rejected = Nystrom.fit(
      input,
      ComponentCount.unsafe(1),
      landmarks = Vector(0, 1),
      kernel = indefinite
    )
    assert(rejected.swap.toOption.exists {
      case MultivarError.InvalidKernelFit(detail) => detail.contains("not certified PSD")
      case _                                      => false
    })

    assert(Nystrom.fit(input, ComponentCount.unsafe(1), Vector(0), tolerance = -1.0).swap.toOption.exists {
      case MultivarError.InvalidTolerance(_, _) => true
      case _                                     => false
    })
  }

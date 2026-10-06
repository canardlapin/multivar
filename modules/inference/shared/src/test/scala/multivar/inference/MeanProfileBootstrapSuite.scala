package multivar.inference

import gale.backend.Backend.given
import gale.linalg.{DMat, DVec, Matrix}
import multivar.core.SpaceId
import resample4s.kernel.Seed

class MeanProfileBootstrapSuite extends munit.FunSuite:
  private def ok[A](value: Either[InferenceError, A]): A = value.fold(e => fail(e.message), identity)
  private def mat(rows: Int, cols: Int, values: String): DMat =
    val xs = if values.isEmpty then Vector.empty else values.split(",").toVector.map(_.toDouble)
    InferenceNumerics.matrixFromRows(Vector.tabulate(rows)(r => Vector.tabulate(cols)(c => xs(r * cols + c))))
  private def close(a: DMat, b: DMat, tolerance: Double = 3e-10): Unit =
    assertEquals((a.rows, a.cols), (b.rows, b.cols))
    for r <- 0 until a.rows; c <- 0 until a.cols do assertEqualsDouble(a(r, c), b(r, c), tolerance)
  private def close(a: DVec, b: DVec, tolerance: Double): Unit =
    assertEquals(a.length, b.length)
    for i <- 0 until a.length do assertEqualsDouble(a(i), b(i), tolerance)
  private def summary(state: ReducedLoadingMoments): ReducedMomentSummary = state.result match
    case Evidence.Computed(value) => value
    case other => fail(s"expected covariance, got $other")
  private val key = StabilityKey(ok(UnitId("axis")), SpaceId.unsafe("profile-coordinates"), StabilityChannel.Loadings)

  test("416 independent raw-cell SVD bootstrap fits preserve beta/FIR maps, moments and context") {
    var checked = 0
    MeanProfileBootstrapOracle.fixtures.foreach { f =>
      val profiles = mat(f.units * f.d, f.p, f.profiles)
      val raw = mat(f.units * f.cells, f.p, f.raw)
      val h = mat(f.cells, f.d, f.h)
      val core = ok(MeanProfileBootstrap.prepare(profiles, f.d, ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced)))
      val direct = ok(MeanProfileBootstrap.prepare(profiles, f.d, ProfileBootstrapOptions(mode = ProfileBootstrapMode.Direct)))
      val observation = ok(core.observed())
      val refMaps = ok(observation.loadingBlock(0, f.p))
      val projection = ok(core.project(raw))
      val voxels = f.p / f.lags
      val lagProjections = Vector.tabulate(f.lags)(lag => ok(core.project(raw.slice(0, raw.rows, lag * voxels, (lag + 1) * voxels), lag * voxels)))
      var moments = Vector.fill(f.d)(ok(ReducedLoadingMoments.empty(key, core.coordinateDimension)))
      var scaledMoments = moments
      val maps = Vector.newBuilder[DMat]
      val scaledMaps = Vector.newBuilder[DMat]
      f.draws.foreach { draw =>
        val counts = ok(ProfileBootstrapCounts.fromIndices(f.units, draw.indices))
        val result = ok(core.refit(counts))
        val reference = ok(direct.refit(counts))
        val expectedS = draw.roots.split(",").map(_.toDouble).toVector
        assertEquals(result.rank, expectedS.length, f.name)
        result.singularValues.toSeq.zip(expectedS).foreach((a, b) => assertEqualsDouble(a, b, 3e-10))
        close(result.singularValues, reference.singularValues, 3e-10)
        val expectedMaps = mat(f.p, result.rank, draw.loadings)
        val fittedMaps = ok(result.loadingBlock(0, f.p))
        close(ok(ComponentAlignment.alignSigns(expectedMaps, fittedMaps)), expectedMaps)
        close(ok(ComponentAlignment.alignSigns(expectedMaps, ok(reference.loadingBlock(0, f.p)))), expectedMaps)
        val aligned = ok(ComponentAlignment.alignSigns(observation.loadingCoordinates, result.loadingCoordinates))
        val alignedExpected = ok(ComponentAlignment.alignSigns(refMaps, expectedMaps))
        close(ok(core.liftBlock(aligned, 0, f.p)), alignedExpected)
        val expectedScaled = Matrix.newBuilder(f.p, result.rank)
        for r <- 0 until f.p; c <- 0 until result.rank do expectedScaled(r, c) = alignedExpected(r, c) * expectedS(c)
        maps += alignedExpected
        scaledMaps += expectedScaled.result()
        moments = moments.indices.map(i => ok(moments(i).add(aligned.col(i)))).toVector
        scaledMoments = scaledMoments.indices.map { i =>
          val scaled = InferenceNumerics.vectorFromSeq(aligned.col(i).toSeq.map(_ * result.singularValues(i)).toVector)
          ok(scaledMoments(i).add(scaled))
        }.toVector
        val scores = ok(projection.scores(result))
        close(scores, raw * fittedMaps)
        val means = Matrix.newBuilder(f.cells, f.p)
        val scoreMeans = Matrix.newBuilder(f.cells, result.rank)
        for cell <- 0 until f.cells do
          for p <- 0 until f.p do means(cell, p) = draw.indices.map(i => raw(i * f.cells + cell, p)).sum / f.units
          for c <- 0 until result.rank do scoreMeans(cell, c) = draw.indices.map(i => scores(i * f.cells + cell, c)).sum / f.units
        val full = scoreMeans.result()
        val rawMeans = means.result()
        close(full, rawMeans * fittedMaps)
        val coordinates = ProfileBootstrapMath.scaled(h.t * rawMeans, 1.0 / Math.sqrt(f.cells))
        val componentScores = coordinates * fittedMaps
        val leftScaled = Matrix.newBuilder(f.d, result.rank)
        for r <- 0 until f.d; c <- 0 until result.rank do leftScaled(r, c) = result.leftDirections(r, c) * result.singularValues(c)
        close(componentScores, leftScaled.result())
        val effect = ProfileBootstrapMath.scaled(h * componentScores, Math.sqrt(f.cells))
        val context = full - effect
        close(context + effect, full)
        val lagScores = lagProjections.map(p => ok(p.scores(result)))
        val sum = lagScores.reduce(_ + _)
        close(sum, scores)
        lagScores.zipWithIndex.foreach { (values, lag) =>
          close(values, raw.slice(0, raw.rows, lag * voxels, (lag + 1) * voxels) * fittedMaps.slice(lag * voxels, (lag + 1) * voxels, 0, result.rank))
        }
        checked += 1
      }
      val retained = maps.result()
      val retainedScaled = scaledMaps.result()
      for (states, values) <- Vector(moments -> retained, scaledMoments -> retainedScaled); c <- states.indices do
        val recovered = ok(summary(states(c)).projectRows(ok(core.basisBlock(0, f.p))))
        for p <- 0 until f.p do
          val mean = values.map(_(p, c)).sum / values.size
          val variance = values.map(v => Math.pow(v(p, c) - mean, 2)).sum / (values.size - 1)
          assertEqualsDouble(recovered.mean(p), mean, 3e-10)
          assertEqualsDouble(recovered.variance(p), variance, 3e-10)
    }
    assertEquals(checked, 416)
  }

  test("original group rank does not truncate bootstrap variation; cancellation is refined") {
    val core = ok(MeanProfileBootstrap.prepare(Matrix(2, 3)(1, 1, 0, 1, -1, 0), 1,
      ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced)))
    val fit = ok(core.refit(ok(ProfileBootstrapCounts.fromIndices(2, Vector(0, 0)))))
    val v = ok(fit.loadingBlock(0, 3))
    assertEqualsDouble(Math.abs(v(1, 0)), 1.0 / Math.sqrt(2), 1e-12)
    val cancellation = ok(MeanProfileBootstrap.prepare(Matrix(2, 3)(1e8, 1, 0, -1e8, 1, 0), 1,
      ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced)))
    val mean = ok(cancellation.observed())
    assertEquals(mean.path, ProfileBootstrapPath.QrDirectRefinement)
    assertEqualsDouble(mean.singularValues(0), 1.0, 1e-10)
    val map = ok(mean.loadingBlock(0, 3))
    assertEqualsDouble(Math.abs(map(1, 0)), 1.0, 1e-10)
    assertEqualsDouble(map(0, 0), 0.0, 1e-10)
    // Re-aggregating already normalized input loses this half-unit difference.
    // Refinement must retain the original supplied entries.
    for mode <- Vector(ProfileBootstrapMode.Direct, ProfileBootstrapMode.Reduced) do
      val shifted = ok(MeanProfileBootstrap.prepare(Matrix(2, 3)(1e8 + 1, 1, 0, -1e8, 1, 0), 1,
        ProfileBootstrapOptions(mode = mode)))
      val refined = ok(shifted.observed())
      assertEqualsDouble(refined.singularValues(0), Math.sqrt(1.25), 1e-10)
      val expected = Matrix(3, 1)(.5 / Math.sqrt(1.25), 1.0 / Math.sqrt(1.25), 0)
      close(ok(refined.loadingBlock(0, 3)), expected, 1e-10)
  }

  test("complete copies precede division so fractional weights do not corrupt cancellation") {
    val z = Matrix(3, 3)(1e8 + 1, 1, 0, -1e8, 1, 0, 0, 1, 0)
    val expected = Matrix(1, 3)(1.0 / 3.0, 1, 0)
    for mode <- Vector(ProfileBootstrapMode.Direct, ProfileBootstrapMode.Reduced) do
      val core = ok(MeanProfileBootstrap.prepare(z, 1, ProfileBootstrapOptions(mode = mode)))
      val fit = ok(core.observed())
      close(fit.leftDirections * ok(fit.loadingBlock(0, 3, scaled = true)).t, expected, 1e-12)
      val repeated = Matrix(3, 3)(1e8 + 1, 1, 0, -2e8, 1, 0, 0, 1, 0)
      val duplicateCore = ok(MeanProfileBootstrap.prepare(repeated, 1, ProfileBootstrapOptions(mode = mode)))
      val duplicate = ok(duplicateCore.refit(ok(ProfileBootstrapCounts.fromIndices(3, Vector(0, 0, 1)))))
      close(duplicate.leftDirections * ok(duplicate.loadingBlock(0, 3, scaled = true)).t, Matrix(1, 3)(2.0 / 3.0, 1, 0), 1e-12)
      val tiny = ok(MeanProfileBootstrap.prepare(Matrix(2, 1)(Double.MinPositiveValue, Double.MinPositiveValue), 1,
        ProfileBootstrapOptions(mode = mode, absoluteZeroTolerance = 0.0)))
      assertEquals(ok(tiny.observed()).singularValues(0), Double.MinPositiveValue)
      val unrepresentable = ok(MeanProfileBootstrap.prepare(Matrix(3, 1)(Double.MinPositiveValue, 0, 0), 1,
        ProfileBootstrapOptions(mode = mode, absoluteZeroTolerance = 0.0)))
      assert(unrepresentable.observed().isLeft)
    assert(MeanProfileBootstrap.prepare(Matrix(2, 2)(1e300, 1e-300, -1e300, 1e-300), 1,
      ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced)).isLeft)
  }

  test("duplicate units, valid zero and tied subspaces keep their rank semantics") {
    val z = Matrix(4, 3)(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
    val core = ok(MeanProfileBootstrap.prepare(z, 2, ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced)))
    val full = ok(core.observed())
    val one = ok(core.refit(ok(ProfileBootstrapCounts.fromIndices(2, Vector(0, 0)))))
    assertEquals((full.rank, one.rank), (2, 1))
    val tied = ok(full.loadingBlock(0, 3))
    val rotation = Matrix(2, 2)(.6, -.8, .8, .6)
    val angles = ok(PrincipalAngles.between(tied, tied * rotation))
    assert(angles.values.forall(_ <= 3e-8))
    val zero = ok(MeanProfileBootstrap.prepare(Matrix(2, 3)(1, 2, 3, -1, -2, -3), 1,
      ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced)))
    assertEquals(ok(zero.observed()).rank, 0)
  }

  test("rank boundaries use stable full SVD and refuse numerically indeterminate cutoffs") {
    for scale <- Vector(1e-8, 1.0, 1e8); ratio <- Vector(5e-7, 2e-6, 1e-6) do
      val z = ProfileBootstrapMath.scaled(Matrix(4, 3)(1, 0, 0, 0, ratio, 0, 1, 0, 0, 0, ratio, 0), scale)
      val core = ok(MeanProfileBootstrap.prepare(z, 2, ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced)))
      if ratio == 1e-6 then assert(core.observed().isLeft)
      else assertEquals(ok(core.observed()).rank, if ratio < 1e-6 then 1 else 2)
  }

  test("full covariance preserves rotated projections; chunked merges match retained moments") {
    val values = Vector.tabulate(37)(i => InferenceNumerics.vectorFromSeq(Vector(i * .2, i * .3 + (i % 3) * .1)))
    val empty = ok(ReducedLoadingMoments.empty(key, 2))
    val all = ok(empty.addBatch(values.iterator))
    val parts = values.grouped(7).map(xs => ok(empty.addBatch(xs.iterator))).toVector
    val combined = parts.foldLeft(empty)((a, b) => ok(a.combine(b)))
    close(summary(all).covariance, summary(combined).covariance, 2e-13)
    val q = Matrix(3, 2)(.6, -.8, .8, .6, .3, -.2)
    val recovered = ok(summary(all).projectRows(q))
    val projected = values.map(v => q * InferenceNumerics.matrixFromRows(v.toSeq.map(x => Vector(x))))
    for i <- 0 until q.rows do
      val mean = projected.map(_(i, 0)).sum / values.size
      val variance = projected.map(x => Math.pow(x(i, 0) - mean, 2)).sum / (values.size - 1)
      assertEqualsDouble(recovered.mean(i), mean, 1e-12)
      assertEqualsDouble(recovered.variance(i), variance, 1e-12)
    val wrong = .36 * summary(all).covariance(0, 0) + .64 * summary(all).covariance(1, 1)
    assert(Math.abs(wrong - recovered.variance(0)) > 1.0)
    assertEquals(empty.count, 0L)
    assert(empty.result.isInstanceOf[Evidence.Unavailable[?]])
    assert(ok(empty.add(values.head)).result.isInstanceOf[Evidence.Unavailable[?]])
  }

  test("counts, preparation, projection, moments and cancellation validate their boundaries") {
    assert(ProfileBootstrapCounts.fromIndices(2, Vector(0)).isLeft)
    assert(ProfileBootstrapCounts.fromIndices(2, Vector(0, 2)).isLeft)
    val action = BootstrapAction.rows(ok(RowCount(7)))
    val a = ok(ProfileBootstrapCounts.fromDraw(ok(action.draw(Seed.fromLong(3), ok(ReplicateId(5))))))
    val b = ok(ProfileBootstrapCounts.fromDraw(ok(action.draw(Seed.fromLong(3), ok(ReplicateId(5))))))
    assertEquals(a.multiplicities, b.multiplicities)
    assertEquals(a.multiplicities.sum, 7)
    val z = Matrix(2, 3)(1, 2, 3, 2, 3, 5)
    assert(MeanProfileBootstrap.prepare(z, 3).isLeft)
    assert(MeanProfileBootstrap.prepare(z, 1, ProfileBootstrapOptions(maxWorkingBytes = 1)).isLeft)
    assert(MeanProfileBootstrap.prepare(z, 1, cancelled = () => true).isLeft)
    var polls = 0
    assert(MeanProfileBootstrap.prepare(z, 1, ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced), () => { polls += 1; polls > 3 }).isLeft)
    assert(MeanProfileBootstrap.prepare(Matrix(2, 1)(1, Double.NaN), 1).isLeft)
    val core = ok(MeanProfileBootstrap.prepare(z, 1))
    assert(core.observed(() => true).isLeft)
    assert(core.project(z, featureFrom = 1).isLeft)
    assert(core.project(z, cancelled = () => true).isLeft)
    assert(core.basisBlock(0, 4).isLeft)
    val another = ok(MeanProfileBootstrap.prepare(z, 1))
    assert(ok(core.project(z)).scores(ok(another.observed())).isLeft)
    val overflowCore = ok(MeanProfileBootstrap.prepare(Matrix(2, 3)(1, 1, 1, 1, 1, 1), 1,
      ProfileBootstrapOptions(mode = ProfileBootstrapMode.Direct)))
    val huge = Matrix(1, 3)(Double.MaxValue, Double.MaxValue, Double.MaxValue)
    assert(ok(overflowCore.project(huge)).scores(ok(overflowCore.observed())).isLeft)
    val empty = ok(ReducedLoadingMoments.empty(key, 2))
    assert(ReducedLoadingMoments.empty(key, Int.MaxValue).isLeft)
    assert(empty.add(InferenceNumerics.vectorFromSeq(Vector(1.0))).isLeft)
    assert(empty.add(InferenceNumerics.vectorFromSeq(Vector(1.0, Double.PositiveInfinity))).isLeft)
    assert(empty.addBatch(Iterator(InferenceNumerics.vectorFromSeq(Vector(1.0, 2.0))), () => true).isLeft)
    assert(empty.combine(ok(ReducedLoadingMoments.empty(key, 3))).isLeft)
  }

  test("public count means preserve complete readouts, cancellation and admission") {
    val input = Matrix(6, 2)(1e16, 1, 10, 20, -1e16, 1, 30, 40, 2, 4, 50, 60)
    val counts = ok(ProfileBootstrapCounts.identity(3))
    close(ok(counts.meanProfiles(input, 2)), Matrix(2, 2)(2.0/3, 2, 30, 40))
    val duplicates = ok(ProfileBootstrapCounts.fromIndices(3, Vector(2, 2, 1)))
    close(ok(duplicates.meanProfiles(input, 2)), Matrix(2, 2)((-1e16+4)/3, 3, 130.0/3, 160.0/3))
    assert(counts.meanProfiles(input, 1).isLeft)
    assert(counts.meanProfiles(input, 2, maxWorkingBytes = 1).isLeft)
    assert(counts.meanProfiles(input, 2, cancelled = () => true).isLeft)
    assert(counts.meanProfiles(Matrix(6, 1)(1, 2, 3, 4, 5, Double.NaN), 2).isLeft)
    var checks = 0
    assert(counts.meanProfiles(input, 2, cancelled = () => { checks += 1; checks > input.rows }).isLeft)
  }

package multivar.inference

import gale.backend.Backend.given
import gale.linalg.{DMat, Matrix}

class ProfileDeletionSuite extends munit.FunSuite:
  private def ok[A](value: Either[InferenceError, A]): A = value.fold(e => fail(e.message), identity)
  private def close(a: DMat, b: DMat, tolerance: Double = 3e-10): Unit =
    assertEquals((a.rows, a.cols), (b.rows, b.cols))
    for r <- 0 until a.rows; c <- 0 until a.cols do
      assertEqualsDouble(a(r, c), b(r, c), tolerance)
  private def reconstructed(fit: MeanProfileFit, p: Int): DMat =
    fit.leftDirections * ok(fit.loadingBlock(0, p, scaled = true)).t

  test("delete one copy preserves domain, source multiset and N-of-N validation") {
    val source = ok(ProfileBootstrapCounts.fromIndices(4, Vector(0, 0, 2, 3)))
    val deletion = ok(source.deleteOne(0))
    assertEquals(deletion.multiplicities, Vector(1, 0, 1, 1))
    assertEquals(deletion.source, source)
    assertEquals((deletion.units, deletion.sampleSize, deletion.removedUnit), (4, 3, 0))
    assertEquals(source.multiplicities, Vector(2, 0, 1, 1))
    assert(source.deleteOne(1).isLeft)
    assert(source.deleteOne(-1).isLeft)
    assert(source.deleteOne(4).isLeft)
    assert(ProfileBootstrapCounts.fromIndices(4, Vector(0, 2, 3)).isLeft)
    val tiny = ok(ProfileBootstrapCounts.identity(2)).deleteOne(0)
    assertEquals(ok(tiny).sampleSize, 1)
  }

  test("all N=3 multisets and deletions match literal profiles, means, maps and readouts") {
    val profiles = Matrix.newBuilder(6, 11)
    for r <- 0 until 6; c <- 0 until 11 do
      profiles(r, c) = Math.sin((r + 1) * (c + 1) * .17) + r * .2 - c * .03
    val z = profiles.result()
    var checked = 0
    for mode <- Vector(ProfileBootstrapMode.Direct, ProfileBootstrapMode.Reduced) do
      val core = ok(MeanProfileBootstrap.prepare(z, 2, ProfileBootstrapOptions(mode = mode)))
      val projected = ok(core.project(z))
      for a <- 0 until 3; b <- 0 until 3; c <- 0 until 3 do
        val copies = Vector(a, b, c)
        val source = ok(ProfileBootstrapCounts.fromIndices(3, copies))
        for removed <- copies.indices do
          val deletion = ok(source.deleteOne(copies(removed)))
          val retained = copies.patch(removed, Vector.empty, 1)
          val literal = z.gatherRows(retained.flatMap(i => Vector(i * 2, i * 2 + 1)))
          val direct = ok(MeanProfileBootstrap.prepare(literal, 2,
            ProfileBootstrapOptions(mode = ProfileBootstrapMode.Direct)))
          val expected = ok(direct.observed())
          val fit = ok(core.refitDeletion(deletion))
          assertEquals(fit.rank, expected.rank)
          close(reconstructed(fit, 11), reconstructed(expected, 11))
          close(ok(deletion.meanProfiles(z, 2)), reconstructed(expected, 11))
          val maps = ok(fit.loadingBlock(0, 11))
          close(ok(ComponentAlignment.alignSigns(ok(expected.loadingBlock(0, 11)), maps)),
            ok(expected.loadingBlock(0, 11)))
          close(ok(projected.scores(fit)), z * maps)
          close(ok(projected.scores(fit, scaled = true)), z * ok(fit.loadingBlock(0, 11, scaled = true)))
          checked += 1
    assertEquals(checked, 162)
  }

  test("deletion fits preserve participant and feature identity under permutation and scale") {
    val z = Matrix.newBuilder(8, 13)
    for r <- 0 until 8; c <- 0 until 13 do z(r, c) = Math.sin((r + 1) * (c + 2) * .13) + r * .1
    val profiles = z.result()
    val copies = Vector(0, 0, 2, 3)
    val units = Vector(2, 0, 3, 1)
    val features = (0 until 13).reverse.toVector
    for mode <- Vector(ProfileBootstrapMode.Direct, ProfileBootstrapMode.Reduced) do
      val options = ProfileBootstrapOptions(mode = mode)
      val baseline = ok(MeanProfileBootstrap.prepare(profiles, 2, options))
      val baseFit = ok(baseline.refitDeletion(ok(ok(ProfileBootstrapCounts.fromIndices(4, copies)).deleteOne(0))))
      val baseMaps = ok(baseFit.loadingBlock(0, 13))
      for scale <- Vector(1e-3, 1.0, 1e3) do
        val changed = Matrix.newBuilder(8, 13)
        for u <- units.indices; d <- 0 until 2; c <- features.indices do
          changed(u * 2 + d, c) = scale * profiles(units(u) * 2 + d, features(c))
        val core = ok(MeanProfileBootstrap.prepare(changed.result(), 2, options))
        val counts = ok(ProfileBootstrapCounts.fromIndices(4, copies.map(units.indexOf)))
        val fit = ok(core.refitDeletion(ok(counts.deleteOne(units.indexOf(0)))))
        val expectedMaps = Matrix.newBuilder(13, baseMaps.cols)
        for r <- features.indices; c <- 0 until baseMaps.cols do expectedMaps(r, c) = baseMaps(features(r), c)
        val expected = expectedMaps.result()
        close(ok(ComponentAlignment.alignSigns(expected, ok(fit.loadingBlock(0, 13)))), expected)
        for axis <- 0 until fit.rank do
          assertEqualsDouble(fit.singularValues(axis) / scale, baseFit.singularValues(axis), 3e-10)
    }

  test("deletion refinement preserves cancellation, tiny values and zero/tied ranks") {
    for mode <- Vector(ProfileBootstrapMode.Direct, ProfileBootstrapMode.Reduced) do
      val z = Matrix(3, 3)(1e8 + 1, 1, 0, -1e8, 1, 0, 4, 5, 6)
      val core = ok(MeanProfileBootstrap.prepare(z, 1, ProfileBootstrapOptions(mode = mode)))
      val deletion = ok(ok(ProfileBootstrapCounts.identity(3)).deleteOne(2))
      val fit = ok(core.refitDeletion(deletion))
      close(reconstructed(fit, 3), Matrix(1, 3)(.5, 1, 0), 2e-10)
      if mode == ProfileBootstrapMode.Reduced then assertEquals(fit.path, ProfileBootstrapPath.QrDirectRefinement)
      val tiny = ok(MeanProfileBootstrap.prepare(Matrix(2, 1)(Double.MinPositiveValue, 0), 1,
        ProfileBootstrapOptions(mode = mode, absoluteZeroTolerance = 0)))
      val one = ok(ok(ProfileBootstrapCounts.identity(2)).deleteOne(1))
      assertEquals(ok(tiny.refitDeletion(one)).singularValues(0), Double.MinPositiveValue)
      val tied = Matrix(6, 3)(1,0,0, 0,1,0, 1,0,0, 0,1,0, 0,0,0, 0,0,0)
      val tiedCore = ok(MeanProfileBootstrap.prepare(tied, 2, ProfileBootstrapOptions(mode = mode)))
      val tiedFit = ok(tiedCore.refitDeletion(deletion))
      assertEquals(tiedFit.rank, 2)
      assertEqualsDouble(tiedFit.singularValues(0), tiedFit.singularValues(1), 2e-12)
      close(reconstructed(tiedFit, 3), Matrix(2, 3)(1,0,0, 0,1,0))
      val zeroCore = ok(MeanProfileBootstrap.prepare(Matrix(3, 2)(1,2, -1,-2, 2,1), 1,
        ProfileBootstrapOptions(mode = mode)))
      assertEquals(ok(zeroCore.refitDeletion(deletion)).rank, 0)
  }

  test("deletions retain owner, source-domain, cancellation, budget and finite-data checks") {
    val z = Matrix(3, 2)(1,2, 3,4, 5,6)
    val core = ok(MeanProfileBootstrap.prepare(z, 1))
    val other = ok(MeanProfileBootstrap.prepare(z, 1))
    val counts = ok(ok(ProfileBootstrapCounts.identity(3)).deleteOne(0))
    val fit = ok(core.refitDeletion(counts))
    assert(ok(other.project(z)).scores(fit).isLeft)
    assert(core.refitDeletion(ok(ok(ProfileBootstrapCounts.identity(4)).deleteOne(0))).isLeft)
    assert(core.refitDeletion(counts, () => true).isLeft)
    assert(counts.meanProfiles(z, 1, cancelled = () => true).isLeft)
    assert(counts.meanProfiles(z, 1, maxWorkingBytes = 1).isLeft)
    assert(counts.meanProfiles(z, 2).isLeft)
    assert(counts.meanProfiles(Matrix(3, 1)(1, Double.NaN, 3), 1).isLeft)
  }

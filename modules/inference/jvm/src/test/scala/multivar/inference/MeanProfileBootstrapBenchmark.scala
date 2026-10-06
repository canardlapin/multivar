package multivar.inference

import gale.linalg.{DMat, Matrix}
import multivar.core.SpaceId
import java.io.{BufferedOutputStream, DataOutputStream, FileOutputStream}
import java.nio.file.{Files, Path}
import java.util.Random

/** Opt-in native timing: identical fixed profiles/counts, mean + variance maps
  * and every singular value. No high-dimensional draws are retained.
  */
object MeanProfileBootstrapBenchmark:
  private def ok[A](value: Either[InferenceError, A]): A = value.fold(e => throw new IllegalStateException(e.message), identity)
  def main(args: Array[String]): Unit =
    require(args.length == 6, "route units coordinates features draws output-prefix")
    val route = args(0)
    require(route == "direct" || route == "reduced")
    val n = args(1).toInt
    val d = args(2).toInt
    val p = args(3).toInt
    val b = args(4).toInt
    require(n >= 2 && d >= 1 && p >= n * d && b >= 2)
    val prefix = args(5)
    val random = new Random(20260908L)
    val builder = Matrix.newBuilder(n * d, p)
    for row <- 0 until n * d; feature <- 0 until p do
      val coordinate = row % d
      builder(row, feature) = Math.sin((feature + 1) * (.013 * (coordinate + 1))) / (coordinate + 1) + .1 * random.nextGaussian()
    val input = builder.result()
    val draws = Vector.fill(b)(ok(ProfileBootstrapCounts.fromIndices(n, Vector.fill(n)(random.nextInt(n)))))
    val started = System.nanoTime()
    val core = ok(MeanProfileBootstrap.prepare(input, d, ProfileBootstrapOptions(
      mode = if route == "direct" then ProfileBootstrapMode.Direct else ProfileBootstrapMode.Reduced,
      expectedReplicates = b, maxWorkingBytes = 2L * 1024 * 1024 * 1024)))
    val reference = ok(core.observed())
    require(reference.rank == d)
    val prepared = System.nanoTime()
    val key = StabilityKey(ok(UnitId("axis")), SpaceId.unsafe("benchmark-core"), StabilityChannel.Loadings)
    var states = if route == "reduced" then Vector.fill(d)(ok(ReducedLoadingMoments.empty(key, core.coordinateDimension))) else Vector.empty
    val mean = new Array[Double](p * d)
    val m2 = new Array[Double](p * d)
    val roots = new Array[Double](b * d)
    var refinementCount = 0
    draws.zipWithIndex.foreach { (counts, index) =>
      val fit = ok(core.refit(counts))
      require(fit.rank == d)
      if fit.path == ProfileBootstrapPath.QrDirectRefinement then refinementCount += 1
      val aligned = ok(ComponentAlignment.alignSigns(reference.loadingCoordinates, fit.loadingCoordinates))
      for component <- 0 until d do roots(index * d + component) = fit.singularValues(component)
      if route == "reduced" then
        states = states.indices.map(c => ok(states(c).add(aligned.col(c)))).toVector
      else
        var feature = 0
        while feature < p do
          var c = 0
          while c < d do
            val position = feature * d + c
            val delta = aligned(feature, c) - mean(position)
            mean(position) += delta / (index + 1)
            m2(position) += delta * (aligned(feature, c) - mean(position))
            c += 1
          feature += 1
    }
    val loopEnded = System.nanoTime()
    if route == "reduced" then
      val summaries = states.map(_.result match
        case Evidence.Computed(value) => value
        case other => throw new IllegalStateException(other.toString))
      var from = 0
      while from < p do
        val until = Math.min(p, from + 256)
        val basis = ok(core.basisBlock(from, until))
        summaries.zipWithIndex.foreach { (summary, c) =>
          val block = ok(summary.projectRows(basis))
          for row <- from until until do
            mean(row * d + c) = block.mean(row - from)
            m2(row * d + c) = block.variance(row - from)
        }
        from = until
    else for i <- m2.indices do m2(i) /= (b - 1)
    val ended = System.nanoTime()
    require(mean.forall(_.isFinite) && m2.forall(x => x.isFinite && x >= 0))
    val output = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(prefix + ".bin")))
    try
      roots.foreach(output.writeDouble)
      mean.foreach(output.writeDouble)
      m2.foreach(output.writeDouble)
    finally output.close()
    def seconds(a: Long, z: Long): Double = (z - a).toDouble / 1e9
    val text = s"""{"route":"$route","units":$n,"coordinates":$d,"features":$p,"draws":$b,"seed":20260908,"coordinateDimension":${core.coordinateDimension},"prepareSeconds":${seconds(started, prepared)},"replicateSeconds":${seconds(prepared, loopEnded)},"liftSeconds":${seconds(loopEnded, ended)},"totalSeconds":${seconds(started, ended)},"estimatedWorkingBytes":${core.preparation.estimatedWorkingBytes},"refinements":$refinementCount,"binaryFormat":"big-endian doubles: roots B*d, mean P*d, variance P*d"}"""
    Files.writeString(Path.of(prefix + ".json"), text + "\n")
    println(s"PASS $route N=$n d=$d P=$p B=$b total=${seconds(started, ended)}")

package multivar.smoke

import gale.linalg.Matrix
import multivar.core.SpaceId
import multivar.inference.*

object MeanProfileBootstrapSmoke:
  def main(args: Array[String]): Unit =
    require(args.isEmpty)
    val core = MeanProfileBootstrap.prepare(Matrix(2, 3)(1, 1, 0, 1, -1, 0), 1,
      ProfileBootstrapOptions(mode = ProfileBootstrapMode.Reduced)).toOption.get
    val mean = ProfileBootstrapCounts.identity(2).toOption.get
      .meanProfiles(Matrix(2, 2)(2, 4, 4, 8), 1).toOption.get
    require(mean(0, 0) == 3 && mean(0, 1) == 6)
    val first = core.refit(ProfileBootstrapCounts.fromIndices(2, Vector(0, 0)).toOption.get).toOption.get
    val second = core.refit(ProfileBootstrapCounts.fromIndices(2, Vector(1, 1)).toOption.get).toOption.get
    val key = StabilityKey(UnitId("axis").toOption.get, SpaceId.unsafe("reduced"), StabilityChannel.Loadings)
    val moments = ReducedLoadingMoments.empty(key, core.coordinateDimension).toOption.get
      .addBatch(Iterator(first.loadingCoordinates.col(0), second.loadingCoordinates.col(0))).toOption.get
    val summary = moments.result match
      case Evidence.Computed(value) => value
      case other => throw new IllegalStateException(other.toString)
    val maps = summary.projectRows(core.basisBlock(0, 3).toOption.get).toOption.get
    require(Math.abs(first.singularValues(0) - Math.sqrt(2)) < 1e-10)
    require(Math.abs(maps.variance(1) - 1.0) < 1e-10)
    val deletion = ProfileBootstrapCounts.identity(2).toOption.get.deleteOne(1).toOption.get
    val deletedFit = core.refitDeletion(deletion).toOption.get
    require(deletion.sampleSize == 1 && Math.abs(deletedFit.singularValues(0) - Math.sqrt(2)) < 1e-10)
    require(Math.abs(JackknifeStandardError.compute(Vector(1.0, 3.0), 2).toOption.get - 1.0) < 1e-14)
    val ids = Vector.tabulate(5)(i => ReplicateId(i).toOption.get)
    val observed = StudentizedValue(10.0, 2.0).toOption.get
    val draws = ids.zip(Vector(7.0, 9.0, 10.0, 11.0, 13.0)).iterator.map { (id, value) =>
      StudentizedDraw(id, StudentizedValue(value, 1.0))
    }
    val interval = StudentizedInterval.compute(observed, ids, draws, StudentizedIntervalOptions(.95)).toOption.get
    require(Math.abs(interval.lower - 4.4) < 1e-12 && Math.abs(interval.upper - 15.6) < 1e-12)
    println("PASS: published deletion refits, jackknife SE and studentized intervals")
    println("PASS: published exact profile bootstrap and reduced map moments")

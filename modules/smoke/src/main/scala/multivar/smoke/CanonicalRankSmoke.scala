package multivar.smoke

import gale.linalg.DMat
import multivar.inference.{Alpha, CanonicalRankSpectrum, GaussianCanonicalRank, GaussianCanonicalRankResult, InferenceError, MonteCarloDraws}
import resample4s.kernel.Seed

object CanonicalRankSmoke:
  def reference(): Either[InferenceError, GaussianCanonicalRankResult] =
    val x = DMat.dense(4,1,Vector(1.0,2.0,3.0,4.0))
    val y = DMat.dense(4,1,Vector(3.0,1.0,4.0,2.0))
    for
      spectrum <- CanonicalRankSpectrum.from(x,y)
      draws <- MonteCarloDraws(19)
      alpha <- Alpha(.05)
      result <- GaussianCanonicalRank.run(spectrum,Seed.fromLong(18237L),draws,alpha)
    yield result

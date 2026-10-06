package multivar.inference

/** Multiplicity over a declared family, preserving its input order. */
object Multiplicity:
  def adjust(raw: Vector[PValue], method: MultiplicityMethod): Either[InferenceError,Vector[PValue]] =
    if raw.isEmpty then Left(InferenceError.InvalidCount("multiplicity family",0))
    else
      adjustArray(raw.map(_.value).toArray,method).foldLeft[Either[InferenceError,Vector[PValue]]](Right(Vector.empty)) {
        (values,p) => for previous <- values; next <- PValue(p) yield previous :+ next
      }

  private[inference] def adjustArray(raw: Array[Double], method: MultiplicityMethod): Array[Double] =
    method match
      case MultiplicityMethod.Bonferroni =>
        raw.map(value => Math.min(1.0, value * raw.length))
      case MultiplicityMethod.Holm =>
        val order = raw.indices.sortBy(raw).toVector
        val out = new Array[Double](raw.length)
        var running = 0.0
        var rank = 0
        while rank < order.length do
          val index = order(rank)
          running = Math.max(running, (raw.length - rank).toDouble * raw(index))
          out(index) = Math.min(1.0, running)
          rank += 1
        out
      case MultiplicityMethod.BenjaminiHochberg =>
        val order = raw.indices.sortBy(raw).toVector
        val out = new Array[Double](raw.length)
        var running = 1.0
        var rank = order.length - 1
        while rank >= 0 do
          val index = order(rank)
          val candidate = raw(index) * raw.length.toDouble / (rank + 1.0)
          running = Math.min(running, candidate)
          out(index) = Math.min(1.0, running)
          rank -= 1
        out

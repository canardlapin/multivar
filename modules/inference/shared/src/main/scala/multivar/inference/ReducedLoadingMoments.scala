package multivar.inference

import gale.linalg.{DMat, DVec, Matrix}

/** Full coordinate covariance, retained so a later linear lift preserves all
  * off-diagonal terms. Immutable states can be combined across worker chunks.
  */
final class ReducedLoadingMoments private (
    val key: StabilityKey,
    val dimension: Int,
    val count: Long,
    private val means: Array[Double],
    private val sums: Array[Double]
):
  def add(value: DVec): Either[InferenceError, ReducedLoadingMoments] = addBatch(Iterator.single(value))

  /** Clone the small accumulator once per batch, not once per replicate. */
  def addBatch(
      values: Iterator[DVec],
      cancelled: () => Boolean = () => false
  ): Either[InferenceError, ReducedLoadingMoments] =
    val mu = means.clone()
    val m2 = sums.clone()
    val delta = new Array[Double](dimension)
    var n = count
    while values.hasNext do
      if cancelled() then return ProfileBootstrapMath.cancellation
      val value = values.next()
      if value.length != dimension then
        return Left(InferenceError.RowCountMismatch("reduced moment vector", dimension, value.length))
      if n == Long.MaxValue then return Left(InferenceError.UnsupportedProblem("bootstrap moment count overflow"))
      n += 1L
      var i = 0
      while i < dimension do
        if !value(i).isFinite then return Left(InferenceError.NonFiniteStatistic("reduced moment value", value(i)))
        delta(i) = value(i) - mu(i)
        mu(i) += delta(i) / n.toDouble
        i += 1
      // Symmetric outer product with (n-1)/n is the Welford cross-moment update.
      val weight = (n - 1L).toDouble / n.toDouble
      i = 0
      while i < dimension do
        var j = 0
        while j <= i do
          val updated = m2(i * dimension + j) + delta(i) * (delta(j) * weight)
          if !updated.isFinite then return Left(InferenceError.NonFiniteStatistic("reduced covariance accumulator", updated))
          m2(i * dimension + j) = updated
          m2(j * dimension + i) = updated
          j += 1
        i += 1
    if mu.exists(x => !x.isFinite) then Left(InferenceError.NumericalFailure("reduced moments", "mean exceeds finite range"))
    else Right(new ReducedLoadingMoments(key, dimension, n, mu, m2))

  def combine(other: ReducedLoadingMoments): Either[InferenceError, ReducedLoadingMoments] =
    if key != other.key || dimension != other.dimension then
      Left(InferenceError.InvalidReplicatePlan("cannot combine different stability keys or coordinate dimensions"))
    else if count == 0L then Right(other)
    else if other.count == 0L then Right(this)
    else if Long.MaxValue - count < other.count then Left(InferenceError.UnsupportedProblem("bootstrap moment count overflow"))
    else
      val n = count + other.count
      val weight = other.count.toDouble / n.toDouble
      val crossWeight = count.toDouble * weight
      val delta = Array.tabulate(dimension)(i => other.means(i) - means(i))
      val mu = Array.tabulate(dimension)(i => means(i) + weight * delta(i))
      val m2 = new Array[Double](sums.length)
      var i = 0
      while i < dimension do
        var j = 0
        while j < dimension do
          m2(i * dimension + j) = sums(i * dimension + j) + other.sums(i * dimension + j) + delta(i) * (delta(j) * crossWeight)
          j += 1
        i += 1
      if mu.exists(x => !x.isFinite) || m2.exists(x => !x.isFinite) then
        Left(InferenceError.NumericalFailure("reduced moment merge", "moments exceed finite range"))
      else Right(new ReducedLoadingMoments(key, dimension, n, mu, m2))

  def result: Evidence[ReducedMomentSummary] =
    if count < 2L then Evidence.Unavailable(UnavailableReason.Unsupported("sample covariance requires at least two valid replicates"))
    else
      val covariance = Matrix.newBuilder(dimension, dimension)
      for i <- 0 until dimension; j <- 0 until dimension do covariance(i, j) = sums(i * dimension + j) / (count - 1L).toDouble
      Evidence.Computed(new ReducedMomentSummary(key, count, InferenceNumerics.vectorFromArray(means), covariance.result()))

object ReducedLoadingMoments:
  def empty(
      key: StabilityKey,
      dimension: Int,
      maxWorkingBytes: Long = 256L * 1024 * 1024
  ): Either[InferenceError, ReducedLoadingMoments] =
    if dimension <= 0 then Left(InferenceError.InvalidCount("reduced moment dimension", dimension))
    else ProfileBootstrapMath.allocation(8 * (BigInt(3) * dimension * dimension + BigInt(5) * dimension), maxWorkingBytes)
      .map(_ => new ReducedLoadingMoments(key, dimension, 0L, new Array[Double](dimension), new Array[Double](dimension * dimension)))

final case class ProjectedMomentSummary(
    key: StabilityKey,
    replicates: Long,
    mean: DVec,
    variance: DVec,
    standardDeviation: DVec
)

final class ReducedMomentSummary private[inference] (
    val key: StabilityKey,
    val replicates: Long,
    val mean: DVec,
    val covariance: DMat
):
  /** Exact moments for rows of a linear lift, with O(r) row scratch. Call with
    * feature blocks to bound output storage. This does not compute quantiles.
    */
  def projectRows(
      lift: DMat,
      cancelled: () => Boolean = () => false
  ): Either[InferenceError, ProjectedMomentSummary] =
    val d = mean.length
    if lift.cols != d then return Left(InferenceError.RowCountMismatch("moment lift coordinates", d, lift.cols))
    val mu = new Array[Double](lift.rows)
    val variance = new Array[Double](lift.rows)
    val sd = new Array[Double](lift.rows)
    var row = 0
    while row < lift.rows do
      if cancelled() then return ProfileBootstrapMath.cancellation
      var v = 0.0
      var compensation = 0.0
      var absolute = 0.0
      var i = 0
      while i < d do
        if !lift(row, i).isFinite then return Left(InferenceError.NonFiniteStatistic("moment lift", lift(row, i)))
        mu(row) += lift(row, i) * mean(i)
        var j = 0
        while j < d do
          val term = lift(row, i) * covariance(i, j) * lift(row, j)
          val next = v + term
          compensation += (if Math.abs(v) >= Math.abs(term) then (v - next) + term else (term - next) + v)
          absolute += Math.abs(term)
          v = next
          j += 1
        i += 1
      v += compensation
      val roundoff = 64.0 * Math.ulp(1.0) * absolute
      if !v.isFinite || !mu(row).isFinite || v < -roundoff then
        return Left(InferenceError.NumericalFailure("moment projection", s"invalid variance $v or mean ${mu(row)}"))
      variance(row) = Math.max(0.0, v)
      sd(row) = Math.sqrt(variance(row))
      row += 1
    Right(ProjectedMomentSummary(key, replicates, InferenceNumerics.vectorFromArray(mu),
      InferenceNumerics.vectorFromArray(variance), InferenceNumerics.vectorFromArray(sd)))

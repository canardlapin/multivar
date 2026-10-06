package multivar
package core

import gale.backend.Backend.given
import gale.linalg.{DMat, Matrix}

/** Stable identities and declared level order, independent of rendered labels. */
final case class FactorialFactor(id: String, levels: Vector[String])

/** A complete categorical reference grid with uniform cell mass.
  * No observations, response weights, centering, or inference are inferred here.
  */
final class FactorialDesign private (
    val factors: Vector[FactorialFactor],
    val cells: Vector[Vector[String]],
    private val positions: Vector[Vector[Int]]
):
  val cellMass: Double = 1.0 / cells.size

  /** Compile a union of distinct marginal terms; the empty subset is the intercept.
    * Overlapping requests remain separate questions, not an additive partition.
    * The budget bounds the K by d basis before allocation, without feature identities.
    */
  def effect(terms: Vector[Vector[String]], maxBasisElements: Long = 1000000L): Either[MultivarError, FactorialEffect] =
    def fail(message: String) = Left(MultivarError.InvalidRowGeometry(message))
    val ids = factors.map(_.id)
    if terms.isEmpty || terms.exists(t => t.distinct.size != t.size || !t.forall(ids.contains)) then
      return fail("Effect requires nonempty distinct valid factor subsets (an empty subset denotes the intercept)")
    val canonical = terms.map(t => ids.filter(t.contains))
    if canonical.distinct.size != canonical.size then return fail("Repeated marginal term in effect union")
    val dimensions = canonical.map(t => factors.map(f => if t.contains(f.id) then f.levels.size - 1 else 1))
    val widths = dimensions.map(_.foldLeft(1L)(_ * _))
    val width = widths.sum
    if maxBasisElements <= 0 || width * cells.size > maxBasisElements || width > Int.MaxValue then
      return fail("Factorial effect exceeds the declared basis-element budget")
    val builder = Matrix.newBuilder(cells.size,width.toInt)
    var offset = 0
    for (term, dims) <- canonical.zip(dimensions) do
      val count = dims.product
      for column <- 0 until count do
        var remainder = column
        val coordinates = dims.indices.reverse.map { j =>
          val value = remainder % dims(j)
          remainder /= dims(j)
          j -> value
        }.toMap
        for row <- cells.indices do
          var value = 1.0
          for j <- factors.indices do
            val levels = factors(j).levels.size
            if !term.contains(factors(j).id) then value /= math.sqrt(levels.toDouble)
            else
              val contrast = coordinates(j) + 1
              val level = positions(row)(j)
              val denominator = math.sqrt(contrast.toDouble * (contrast + 1))
              value *= (if level < contrast then 1.0 / denominator
                else if level == contrast then -contrast.toDouble / denominator else 0.0)
          builder(row,offset + column) = value
      offset += count
    Right(new FactorialEffect(canonical,builder.result(),cellMass))

object FactorialDesign:
  def compile(factors: Vector[FactorialFactor], cells: Vector[Vector[String]], maxCells: Int = 4096): Either[MultivarError, FactorialDesign] =
    def fail(message: String) = Left(MultivarError.InvalidRowGeometry(message))
    if factors.isEmpty || factors.map(_.id).distinct.size != factors.size ||
        factors.exists(f => f.id.trim.isEmpty || f.levels.size < 2 || f.levels.distinct.size != f.levels.size || f.levels.exists(_.trim.isEmpty)) then
      return fail("Factors require unique nonempty IDs and at least two unique nonempty levels")
    // Bound the product incrementally: no overflowing factorial grid or phantom cells.
    var size = 1L
    var factorIndex = 0
    while factorIndex < factors.size do
      size *= factors(factorIndex).levels.size
      if maxCells <= 0 || size > maxCells then return fail("Factorial grid exceeds the declared cell budget")
      factorIndex += 1
    if cells.size != size || cells.distinct.size != cells.size || cells.exists(c =>
        c.size != factors.size || c.zip(factors).exists((level,factor) => !factor.levels.contains(level))) then
      return fail("Cells must identify the complete Cartesian grid exactly once")
    val positions = cells.map(c => c.zip(factors).map((level,factor) => factor.levels.indexOf(level)))
    Right(new FactorialDesign(factors,cells,positions))

/** Compact normalized effect coordinates: R = sqrt(cellMass) Q' Y.
  * Multiplication delegates to Gale. No feature by feature identity or covariance
  * is constructed; callers can pass R directly to an existing small-side SVD.
  */
final class FactorialEffect private[core] (
    val terms: Vector[Vector[String]],
    val basis: DMat,
    val cellMass: Double
):
  val degreesOfFreedom: Int = basis.cols
  private def validate(response: DMat): Either[MultivarError,Unit] =
    var row = 0
    while row < response.rows do
      var column = 0
      while column < response.cols do
        if !response(row,column).isFinite then
          return Left(MultivarError.NonFiniteValue("factorial response",row*response.cols+column,response(row,column)))
        column += 1
      row += 1
    Right(())

  def coordinates(response: DMat): Either[MultivarError,DMat] =
    if response.rows != basis.rows || response.cols == 0 then
      Left(MultivarError.MatrixShapeMismatch("Response rows must match the compiled reference grid"))
    else validate(response).map { _ =>
      val result = basis.t * response
      val builder = Matrix.newBuilder(result.rows,result.cols)
      for r <- 0 until result.rows; c <- 0 until result.cols do builder(r,c) = result(r,c) * math.sqrt(cellMass)
      builder.result()
    }
  def contribution(response: DMat): Either[MultivarError,DMat] =
    if response.rows != basis.rows || response.cols == 0 then
      Left(MultivarError.MatrixShapeMismatch("Response rows must match the compiled reference grid"))
    else validate(response).map(_ => basis * (basis.t * response))

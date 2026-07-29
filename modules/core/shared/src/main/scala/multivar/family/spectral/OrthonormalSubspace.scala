package multivar
package family.spectral

import multivar.core.*

import gale.linalg.DMat

/** Orthonormal column basis for a subspace, applied without materialising \(EE^\top\).
  *
  * Prefer this over Gram-eigendecomposition projectors when the generators are
  * narrow scientific designs: range extraction uses SVD, and projection is
  * \(Z \mapsto E(E^\top Z)\).
  */
private[multivar] final class OrthonormalSubspace private (
    val basis: DMat
):
  def rank: Int =
    basis.cols

  def dim: Int =
    basis.rows

  def project(input: DMat): Either[MultivarError, DMat] =
    if input.rows != dim then
      Left(
        MultivarError.MatrixShapeMismatch(
          s"orthonormal subspace of dimension $dim cannot project a ${input.rows}x${input.cols} matrix"
        )
      )
    else if rank == 0 then Right(DMat.zeros(input.rows, input.cols))
    else
      Right(
        GaleNumerics.multiply(
          basis,
          GaleNumerics.multiply(basis.transpose, input)
        )
      )

  def residual(input: DMat): Either[MultivarError, DMat] =
    project(input).map(focused => MatrixOps.subtract(input, focused))

  def split(input: DMat): Either[MultivarError, OrthonormalSubspace.Split] =
    for
      focus <- project(input)
      remainder = MatrixOps.subtract(input, focus)
    yield OrthonormalSubspace.Split(focus, remainder)

private[multivar] object OrthonormalSubspace:
  final case class Split(focus: DMat, remainder: DMat)

  def empty(dim: Int): OrthonormalSubspace =
    require(dim >= 0, "subspace dimension must be non-negative")
    new OrthonormalSubspace(DMat.zeros(dim, 0))

  /** Orthonormal basis for the column space of `generators`, truncated by singular-value support. */
  def span(
      generators: DMat,
      solver: SvdSolver,
      tolerance: Double
  ): Either[MultivarError, OrthonormalSubspace] =
    if generators.cols == 0 || generators.rows == 0 then Right(empty(generators.rows))
    else
      val limit = math.min(generators.rows, generators.cols)
      for
        _ <- RowGeometryOps.requireTolerance("orthonormal subspace tolerance", tolerance)
        request <- ComponentCount(limit)
        svd <- solver.decompose(MatrixView.dense(generators), request)
        kept = supportedRank(svd.singularValues, tolerance)
      yield
        if kept == 0 then empty(generators.rows)
        else new OrthonormalSubspace(MatrixOps.takeColumns(svd.u, kept))

  private def supportedRank(values: gale.linalg.DVec, tolerance: Double): Int =
    var largest = 0.0
    var i = 0
    while i < values.length do
      val absolute = math.abs(values(i))
      if absolute > largest then largest = absolute
      i += 1
    val cutoff = tolerance * math.max(1.0, largest)
    var kept = 0
    i = 0
    while i < values.length do
      if math.abs(values(i)) > cutoff then kept += 1
      i += 1
    kept

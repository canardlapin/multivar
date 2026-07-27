package multivar
package core

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.MutableDVec

/** Internal adapter for semantic operators whose natural implementation acts
  * on a block of parameter columns.
  *
  * The public numerical contract remains Gale's `DoubleLinearOperator`; this
  * adapter only preserves block-aware downstream implementations.
  */
private[multivar] trait MatrixActionOperator extends DoubleLinearOperator:
  protected def forwardMatrix(input: DMat): Either[LinAlgError, DMat]
  protected def transposeMatrix(input: DMat): Either[LinAlgError, DMat]

  final override def applyTo(input: DMat): Either[LinAlgError, DMat] =
    if input.rows != cols then
      Left(LinAlgError.InvalidArgument(s"operator expected $cols input rows, got ${input.rows}"))
    else forwardMatrix(input)

  final override def transposeApplyTo(input: DMat): Either[LinAlgError, DMat] =
    if input.rows != rows then
      Left(LinAlgError.InvalidArgument(s"adjoint expected $rows input rows, got ${input.rows}"))
    else transposeMatrix(input)

  final override def applyTo(input: DVec, output: MutableDVec): Unit =
    applyColumn(input, output, cols, rows, forwardMatrix)

  final override def transposeApplyTo(input: DVec, output: MutableDVec): Unit =
    applyColumn(input, output, rows, cols, transposeMatrix)

  private def applyColumn(
      input: DVec,
      output: MutableDVec,
      expectedInput: Int,
      expectedOutput: Int,
      action: DMat => Either[LinAlgError, DMat]
  ): Unit =
    if input.length != expectedInput then
      throw LinAlgError.VectorLengthMismatch(expectedInput, input.length)
    if output.length != expectedOutput then
      throw LinAlgError.VectorLengthMismatch(expectedOutput, output.length)
    val column = DMat.tabulate(input.length, 1): (row, _) =>
      input(row)
    val result = action(column).fold(error => throw error, identity)
    if result.rows != expectedOutput || result.cols != 1 then
      throw LinAlgError.InvalidArgument(
        s"operator returned ${result.rows}x${result.cols}, expected ${expectedOutput}x1"
      )
    var row = 0
    while row < expectedOutput do
      output(row) = result(row, 0)
      row += 1

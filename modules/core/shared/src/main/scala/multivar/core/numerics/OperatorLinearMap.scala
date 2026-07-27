package multivar
package core

import gale.linalg.DoubleLinearOperator
import gale.linalg.DVec
import gale.linalg.MutableDVec

/** Narrow adapter from the semantic operator algebra to the reusable linalg
  * operator capability. Program validation and solver compilation share it.
  */
private[multivar] final class OperatorLinearMap[
    From <: Coordinate,
    To <: Coordinate,
    Role <: OperatorRoleTag,
    Evidence <: OperatorEvidence
](
    operator: Op[From, To, Role, Evidence]
) extends DoubleLinearOperator:
  private val delegate = operator.kernel.linearMap

  val rows: Int = delegate.rows
  val cols: Int = delegate.cols

  def applyTo(input: DVec, output: MutableDVec): Unit =
    delegate.applyTo(input, output)

  override def transposeApplyTo(input: DVec, output: MutableDVec): Unit =
    delegate.transposeApplyTo(input, output)

  override def adjoint: DoubleLinearOperator =
    delegate.adjoint

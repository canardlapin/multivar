package multivar.smoke

import gale.linalg.Matrix
import multivar.inference.*
import resample4s.kernel.Seed

/** Runs against publishedLocal artifacts, without source-project dependencies. */
object SignFlipSmoke:
  private def checked[A](value: Either[InferenceError,A]): A =
    value.fold(e => throw IllegalArgumentException(e.message),identity)
  def main(args: Array[String]): Unit =
    require(args.isEmpty)
    val rows = checked(RowCount(4))
    val action = SignFlipAction.forPartition(checked(RowPartition.from(rows,Vector(Vector(0,2),Vector(1,3)))))
    val draw = checked(action.draw(Seed.fromLong(17L),checked(ReplicateId(3))))
    val data = Matrix(4,2)(1,2,3,4,5,6,7,8)
    val result = checked(draw.applyTo(data))
    require(draw.signs(0) == draw.signs(2) && draw.signs(1) == draw.signs(3))
    for r <- 0 until 4; c <- 0 until 2 do require(result(r,c) == data(r,c)*draw.signs(r))
    println("PASS: published sign-flip action preserves block rows and transforms the matrix")

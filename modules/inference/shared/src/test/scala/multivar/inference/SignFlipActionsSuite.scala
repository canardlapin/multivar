package multivar.inference

import resample4s.kernel.Seed

class SignFlipActionsSuite extends munit.FunSuite:
  private def accepted[A](value: Either[InferenceError,A]): A = value.fold(e => fail(e.message),identity)
  private def pattern(values: Int*): RowSignFlip = accepted(RowSignFlip.from(values.toVector))

  test("explicit signs form the complete row-sign group and preserve every column") {
    val input = InferenceNumerics.matrixFromRows(Vector(Vector(2.0,-3.0),Vector(5.0,7.0),Vector(-11.0,13.0)))
    val orbit = Vector.tabulate(8)(mask => accepted(RowSignFlip.from(Vector.tabulate(3)(r => if (mask & (1 << r)) == 0 then 1 else -1))))
    orbit.foreach { action =>
      val transformed = accepted(action.applyTo(input))
      for r <- 0 until 3; c <- 0 until 2 do
        assertEquals(transformed(r,c),input(r,c)*action.signs(r))
      assertEquals(accepted(action.compose(action)).signs,Vector(1,1,1))
      orbit.foreach { other =>
        val combined = accepted(accepted(action.compose(other)).applyTo(input))
        val sequential = accepted(action.applyTo(accepted(other.applyTo(input))))
        for r <- 0 until 3; c <- 0 until 2 do assertEquals(combined(r,c),sequential(r,c))
      }
    }
  }

  test("Scala callers cannot bypass the sign smart constructor") {
    assert(compileErrors("new multivar.inference.RowSignFlip(multivar.inference.RowCount(1).toOption.get, Vector(0))").nonEmpty)
  }

  test("invalid signs and mismatched shapes are refused") {
    assert(RowSignFlip.from(Vector.empty).isLeft)
    assert(RowSignFlip.from(Vector(1,0)).isLeft)
    assert(RowSignFlip.from(Vector(1,2)).isLeft)
    assert(pattern(1,-1).compose(pattern(1)).isLeft)
    assert(pattern(1,-1).applyTo(InferenceNumerics.matrixFromRows(Vector(Vector(1.0)))).isLeft)
  }

  test("partition signs preserve noncontiguous row groups and replay by replicate") {
    val partition = accepted(RowPartition.from(accepted(RowCount(4)),Vector(Vector(0,2),Vector(1,3))))
    val action = SignFlipAction.forPartition(partition)
    val root = Seed.fromLong(78123L)
    val draws = Vector.tabulate(256) { index =>
      val replicate = accepted(ReplicateId(index))
      val draw = accepted(action.draw(root,replicate))
      assertEquals(draw.signs,accepted(action.draw(root,replicate)).signs)
      assertEquals(draw.signs(0),draw.signs(2))
      assertEquals(draw.signs(1),draw.signs(3))
      draw.signs
    }
    assertEquals(draws.distinct.size,4)
    val reversed = (0 until 256).reverse.map(i => accepted(action.draw(root,accepted(ReplicateId(i)))).signs).reverse.toVector
    assertEquals(draws,reversed)
  }

  test("independent row signs explore the full orbit without enforced balance") {
    val action = accepted(SignFlipAction.independent(accepted(RowCount(4))))
    val draws = Vector.tabulate(256)(i => accepted(action.draw(Seed.fromLong(913L),accepted(ReplicateId(i)))).signs)
    assertEquals(draws.distinct.size,16)
    assert(draws.contains(Vector(1,1,1,1)))
    assert(draws.contains(Vector(-1,-1,-1,-1)))
  }

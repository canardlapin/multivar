package multivar.core

import gale.backend.Backend.given
import gale.linalg.{DMat, Matrix}

class FactorialDesignSuite extends munit.FunSuite:
  private def checked[A](value: Either[MultivarError,A]): A = value.fold(e => fail(e.message),identity)
  private val factors = Vector(FactorialFactor("A",Vector("0","1")),FactorialFactor("B",Vector("low","mid","high")))
  private val cells = for a <- factors(0).levels; b <- factors(1).levels yield Vector(a,b)
  private val data = Matrix(6,2)(10,1,12,4,14,2,20,8,26,3,21,6)
  private def close(a: DMat,b: DMat): Unit =
    assertEquals((a.rows,a.cols),(b.rows,b.cols))
    for r <- 0 until a.rows; c <- 0 until a.cols do assertEqualsDouble(a(r,c),b(r,c),1e-11)
  test("all marginal terms match independent cell averaging and reconstruct") {
    val design = checked(FactorialDesign.compile(factors,cells))
    val terms = Vector(Vector.empty[String],Vector("A"),Vector("B"),Vector("A","B"))
    val contributions = terms.map(t => checked(checked(design.effect(Vector(t))).contribution(data)))
    for r <- cells.indices; p <- 0 until data.cols do
      val grand = cells.indices.map(data(_,p)).sum / 6
      val a = cells.indices.filter(i => cells(i)(0) == cells(r)(0)).map(data(_,p)).sum / 3 - grand
      val b = cells.indices.filter(i => cells(i)(1) == cells(r)(1)).map(data(_,p)).sum / 2 - grand
      val expected = Vector(grand,a,b,data(r,p)-grand-a-b)
      for e <- terms.indices do assertEqualsDouble(contributions(e)(r,p),expected(e),1e-11)
    close(contributions.reduce(_ + _),data)
    val union = checked(design.effect(Vector(Vector("A"),Vector("A","B"))))
    close(checked(union.contribution(data)),contributions(1)+contributions(3))
    close(union.basis.t*union.basis,DMat.eye(3))
  }
  test("cell order and level coding change coordinates but preserve effects") {
    val order = Vector(5,2,0,4,1,3)
    val reordered = Matrix.newBuilder(6,2)
    for r <- order.indices; c <- 0 until 2 do reordered(r,c) = data(order(r),c)
    val first = checked(checked(FactorialDesign.compile(factors,cells)).effect(Vector(Vector("B"))))
    val second = checked(checked(FactorialDesign.compile(factors.map(f => f.copy(levels=f.levels.reverse)),order.map(cells))).effect(Vector(Vector("B"))))
    val expected = checked(first.contribution(data))
    val actual = checked(second.contribution(reordered.result()))
    for r <- order.indices; c <- 0 until 2 do assertEqualsDouble(actual(r,c),expected(order(r),c),1e-11)
  }
  test("invalid grids, duplicate terms, unknown factors and budgets are refused") {
    assert(FactorialDesign.compile(factors,cells.drop(1)).isLeft)
    assert(FactorialDesign.compile(factors,cells.updated(0,cells(1))).isLeft)
    assert(FactorialDesign.compile(factors,cells,maxCells=5).isLeft)
    assert(FactorialDesign.compile(Vector.fill(100)(FactorialFactor("A",Vector("a","b"))),cells).isLeft)
    val d = checked(FactorialDesign.compile(factors,cells))
    assert(d.effect(Vector(Vector("Q"))).isLeft)
    assert(d.effect(Vector(Vector("A","B"),Vector("B","A"))).isLeft)
    assert(d.effect(Vector(Vector("A","B")),maxBasisElements=11).isLeft)
  }
  test("overflowed effect coordinates are refused as typed nonfinite output") {
    val design = checked(FactorialDesign.compile(Vector(FactorialFactor("A", Vector("low", "high"))),
      Vector(Vector("low"), Vector("high"))))
    val intercept = checked(design.effect(Vector(Vector.empty)))
    intercept.coordinates(Matrix(2, 1)(Double.MaxValue, Double.MaxValue)) match
      case Left(MultivarError.NonFiniteValue("factorial coordinates", _, value)) => assert(!value.isFinite)
      case other => fail(s"expected typed nonfinite factorial coordinate error, got $other")
  }

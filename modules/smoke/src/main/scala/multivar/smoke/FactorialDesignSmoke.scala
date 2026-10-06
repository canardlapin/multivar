package multivar.smoke

import multivar.core.{FactorialDesign,FactorialFactor}
import gale.linalg.Matrix

object FactorialDesignSmoke:
  def main(args: Array[String]): Unit =
    require(args.isEmpty)
    val d = FactorialDesign.compile(Vector(FactorialFactor("A",Vector("a","b"))),Vector(Vector("a"),Vector("b"))).toOption.get
    val e = d.effect(Vector(Vector("A"))).toOption.get
    val result = e.contribution(Matrix(2,1)(3,7)).toOption.get
    require(math.abs(result(0,0)+2) < 1e-12 && math.abs(result(1,0)-2) < 1e-12)
    println("PASS: published factorial coordinates preserve the centered mean effect")

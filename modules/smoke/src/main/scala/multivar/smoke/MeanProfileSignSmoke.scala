package multivar.smoke

import gale.linalg.Matrix
import multivar.inference.*

object MeanProfileSignSmoke:
  def main(args: Array[String]): Unit =
    require(args.isEmpty)
    val test = MeanProfileSignTest.prepare(Matrix(3,1)(1.0,2.0,3.0)).toOption.get
    val result = test.run(MeanProfileSampling.Exact).toOption.get
    require(result.pValue.value == .25 && result.draws.size == 8 && result.identityOccurrences == 1)
    val adjusted = Multiplicity.adjust(Vector(result.pValue),MultiplicityMethod.Holm).toOption.get
    require(adjusted.head.value == .25)
    println("PASS: published mean-profile sign test and family adjustment")

package multivar
package analysis

import scala.compiletime.testing.typeCheckErrors

class AnalysisSurfaceSuite extends munit.FunSuite:

  test("the analysis façade resolves the ordinary dense entry points"):
    val errors = typeCheckErrors("""
      import multivar.analysis.*

      val pca: Pca.type = Pca
      val svd: Svd.type = Svd
      val gpca: Gpca.type = Gpca
      val plsc: Plsc.type = Plsc
      val cca: Cca.type = Cca
      val rrr: ReducedRankRegression.type = ReducedRankRegression
      val lda: FisherDiscriminant.type = FisherDiscriminant
      val within: WithinScatterPolicy = WithinScatterPolicy.defaultTraceScaled
      val gpcaCentering: GpcaCentering = GpcaCentering.Auto
      val symmetry: KernelSymmetryPolicy = KernelSymmetryPolicy.SymmetrizeWithin()
      val cpca: Cpca.type = Cpca
      val nystrom: Nystrom.type = Nystrom
      val preproc: PreprocessSpec = PreprocessSpec.Center
      val convention: VarianceConvention = VarianceConvention.Sample
      val err: MultivarError = MultivarError.InvalidDimension("component count", 0)
    """)

    assertEquals(errors, List.empty)

  test("importing only analysis does not bring internal types into scope"):
    val semanticSpace = typeCheckErrors("""
      import multivar.analysis.*
      val space = null.asInstanceOf[SemanticSpace]
    """)
    val op = typeCheckErrors("""
      import multivar.analysis.*
      val operator = null.asInstanceOf[Op[?, ?, ?, ?]]
    """)
    val programFit = typeCheckErrors("""
      import multivar.analysis.*
      val fit = null.asInstanceOf[OperatorProgramFit]
    """)
    val svdResult = typeCheckErrors("""
      import multivar.analysis.*
      val result = null.asInstanceOf[SvdResult]
    """)
    val frame = typeCheckErrors("""
      import multivar.analysis.*
      val transform = null.asInstanceOf[FittedFrameTransform]
    """)
    val certificate = typeCheckErrors("""
      import multivar.analysis.*
      val evidence = null.asInstanceOf[VariationalFrameCertificate]
    """)

    assert(semanticSpace.nonEmpty)
    assert(op.nonEmpty)
    assert(programFit.nonEmpty)
    assert(svdResult.nonEmpty)
    assert(frame.nonEmpty)
    assert(certificate.nonEmpty)

  test("the analysis façade does not resurrect the flat root namespace"):
    val errors = typeCheckErrors("""
      import multivar.analysis.*
      val pca = multivar.Pca
    """)

    assert(errors.nonEmpty)

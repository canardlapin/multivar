package multivar
package workflow

import scala.compiletime.testing.typeCheckErrors

class PackageHierarchySuite extends munit.FunSuite:

  test("public concepts resolve from their semantic owners"):
    val errors = typeCheckErrors("""
      val matrix = null.asInstanceOf[multivar.core.MatrixView]
      val contract = null.asInstanceOf[multivar.contract.MathematicalModelContract]
      val program = null.asInstanceOf[multivar.optimization.OperatorProgram]
      val evidence = null.asInstanceOf[multivar.solver.VariationalFrameCertificate]
      val fitted = null.asInstanceOf[
        multivar.lifecycle.FittedModel[?, ?, ?]
      ]
      val projection = null.asInstanceOf[multivar.capability.FittedFrameTransform]
      val glrm = multivar.family.glrm.GeneralizedLowRankProgram
      val multiblock = multivar.family.multiblock.ExactMultiblockPrograms
      val spec = multivar.workflow.ModelSpec
    """)

    assertEquals(errors, List.empty)

  test("the former flat namespace cannot silently regrow"):
    val errors = typeCheckErrors("""
      val matrix = null.asInstanceOf[multivar.MatrixView]
      val contract = null.asInstanceOf[multivar.MathematicalModelContract]
      val program = null.asInstanceOf[multivar.OperatorProgram]
      val fitted = null.asInstanceOf[multivar.FittedModel[?, ?, ?]]
      val glrm = multivar.GeneralizedLowRankProgram
      val spec = multivar.ModelSpec
    """)

    assert(errors.nonEmpty)

  test("family composition is named at the owning vertical boundary"):
    val obsoleteMultiblockEntry = typeCheckErrors("""
      multivar.family.spectral.ExactSpectralPrograms.multisetQuadratic
    """)
    val pairedTransferInGenericCapability = typeCheckErrors("""
      val transfer = multivar.capability.PairedTransfer
    """)
    val currentOwners = typeCheckErrors("""
      val multiblock = multivar.family.multiblock.ExactMultiblockPrograms
      val transfer = multivar.family.paired.PairedTransfer
    """)

    assert(obsoleteMultiblockEntry.nonEmpty)
    assert(pairedTransferInGenericCapability.nonEmpty)
    assertEquals(currentOwners, List.empty)

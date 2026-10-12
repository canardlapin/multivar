package multivar.optimization

import gale.optim.FirstOrderStoppingStatus
import multivar.contract.OptimizationClaimClass
import multivar.core.CertificateTolerance

class GaleStoppingSuite extends munit.FunSuite:
  private val certificate = SplitResidualCertificate(0.0, 0.0, 0.0, 0.0, 1, CertificateTolerance.strict)

  test("execution stops cannot acquire convergence claims even with zero residuals"):
    val stops = Vector(
      FirstOrderStoppingStatus.NumericalStagnation,
      FirstOrderStoppingStatus.EvaluationLimit,
      FirstOrderStoppingStatus.Cancelled,
      FirstOrderStoppingStatus.LineSearchFailed
    )
    for status <- stops do
      assertEquals(GaleStopping.termination(status), NumericalTermination.NumericalFailure)
      assertEquals(GaleStopping.stationaryClaim(status), OptimizationClaimClass.Unresolved)
      assertEquals(GaleStopping.split(status, certificate), SplitStoppingStatus.NumericalFailure(s"Gale stopped: $status"))

  test("convergence and iteration-limit evidence remain distinct"):
    assertEquals(GaleStopping.termination(FirstOrderStoppingStatus.Converged), NumericalTermination.Converged)
    assertEquals(GaleStopping.stationaryClaim(FirstOrderStoppingStatus.Converged), OptimizationClaimClass.Stationary)
    assertEquals(GaleStopping.split(FirstOrderStoppingStatus.Converged, certificate), SplitStoppingStatus.Converged(certificate))
    assertEquals(GaleStopping.termination(FirstOrderStoppingStatus.IterationLimit), NumericalTermination.IterationLimit)
    assertEquals(GaleStopping.stationaryClaim(FirstOrderStoppingStatus.IterationLimit), OptimizationClaimClass.Unresolved)
    assertEquals(GaleStopping.split(FirstOrderStoppingStatus.IterationLimit, certificate), SplitStoppingStatus.IterationLimit(certificate))

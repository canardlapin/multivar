package multivar.optimization

import gale.optim.FirstOrderStoppingStatus
import multivar.contract.OptimizationClaimClass

/** Preserve non-convergence when adapting Gale's execution stops to semantic evidence. */
private[multivar] object GaleStopping:
  def termination(status: FirstOrderStoppingStatus): NumericalTermination = status match
    case FirstOrderStoppingStatus.Converged => NumericalTermination.Converged
    case FirstOrderStoppingStatus.IterationLimit => NumericalTermination.IterationLimit
    case FirstOrderStoppingStatus.NumericalStagnation | FirstOrderStoppingStatus.EvaluationLimit |
        FirstOrderStoppingStatus.Cancelled | FirstOrderStoppingStatus.LineSearchFailed =>
      NumericalTermination.NumericalFailure

  def split(status: FirstOrderStoppingStatus, certificate: SplitResidualCertificate): SplitStoppingStatus =
    termination(status) match
      case NumericalTermination.Converged => SplitStoppingStatus.Converged(certificate)
      case NumericalTermination.IterationLimit => SplitStoppingStatus.IterationLimit(certificate)
      case NumericalTermination.NumericalFailure | NumericalTermination.Infeasible =>
        SplitStoppingStatus.NumericalFailure(s"Gale stopped: $status")

  def stationaryClaim(status: FirstOrderStoppingStatus): OptimizationClaimClass =
    if termination(status) == NumericalTermination.Converged then OptimizationClaimClass.Stationary
    else OptimizationClaimClass.Unresolved

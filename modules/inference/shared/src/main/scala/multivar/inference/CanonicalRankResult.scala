package multivar.inference

import gale.spectral.SpectralDiagnostics
import resample4s.kernel.Seed

enum CanonicalRankMethod:
  case ScoreOrthogonalPermutationV2
  case GaussianInterlacingWilksV1

  def identity: String = this match
    case ScoreOrthogonalPermutationV2 => "canonical-rank-score-orthogonal-permutation/v2"
    case GaussianInterlacingWilksV1 => "canonical-rank-gaussian-interlacing-wilks/v1"

/** Arithmetic shared by two explicitly different reference laws. Neither
  * supplies downstream study qualification or verifies distributional assumptions. */
sealed trait CanonicalRankResult:
  def method: CanonicalRankMethod
  def correlations: Vector[Double]
  def receipts: Vector[MonteCarloReceipt]
  def adjustedPValues: Vector[PValue]
  def detectableRank: Int
  def completedCompactFits: Long
  def observedCanonicalDiagnostics: SpectralDiagnostics
  def qrRankTolerance: Double

final class FixedCanonicalRankResult private[inference] (
    val correlations: Vector[Double], val receipts: Vector[MonteCarloReceipt],
    val adjustedPValues: Vector[PValue], val detectableRank: Int,
    val refitMethod: CanonicalRankRefitMethod, val completedCompactFits: Long, val identityDraws: Int,
    val sampling: CanonicalRankSampling, val action: PermutationAction,
    val candidateDraws: Int, val duplicateDraws: Int, val sampledCandidateIds: Vector[ReplicateId],
    val observedCanonicalDiagnostics: SpectralDiagnostics, val qrRankTolerance: Double,
    val maximumNullSvdResidual: Double, val maximumNullOrthogonalityError: Double,
    val allNullFitsExtremeCertified: Boolean
) extends CanonicalRankResult:
  val method: CanonicalRankMethod = CanonicalRankMethod.ScoreOrthogonalPermutationV2


/** Conservative finite-Gaussian reference, conditional on iid zero-mean joint
  * Gaussian rows and fixed candidate spaces. Not a fitted-tail permutation law. */
final class GaussianCanonicalRankResult private[inference] (
    val correlations: Vector[Double], val receipts: Vector[MonteCarloReceipt],
    val adjustedPValues: Vector[PValue], val detectableRank: Int,
    val references: Vector[GaussianCanonicalReference], val hypothesisSeeds: Vector[Seed],
    val completedCompactFits: Long, val plannedElements: Long,
    val observedCanonicalDiagnostics: SpectralDiagnostics, val qrRankTolerance: Double,
    val maximumNullSvdResidual: Double, val maximumNullOrthogonalityError: Double
) extends CanonicalRankResult:
  val method: CanonicalRankMethod = CanonicalRankMethod.GaussianInterlacingWilksV1

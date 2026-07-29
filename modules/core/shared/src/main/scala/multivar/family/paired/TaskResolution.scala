package multivar
package family.paired

import gale.linalg.DMat

/** Selects a scientifically meaningful subspace of a least-squares task.
  *
  * Designs are supplied in original working coordinates. Metric standardization
  * and spectral transport are internal.
  */
sealed trait TaskResolution:
  def label: String

object TaskResolution:
  /** Analyze the complete task without an external split. */
  case object Whole extends TaskResolution:
    val label: String = "whole task"

  /** Constraints stated in the target (response / feature) space. */
  object Target:
    /** `design` specifies combinations of target variables (Takane weights). */
    def weights(design: DMat): TaskResolution =
      TargetWeights(design)

    /** `design` specifies target covariance or loading patterns (Takane structure). */
    def structure(
        design: DMat,
        dual: DualSpectrum = DualSpectrum.moorePenrose
    ): TaskResolution =
      TargetStructure(design, dual)

  /** External information stated in the source space (Phase 2). */
  object Source:
    def covariance(design: DMat): TaskResolution =
      SourceCovariance(design)

    def regression(design: DMat): TaskResolution =
      SourceRegression(design)

  private[paired] final case class TargetWeights(design: DMat) extends TaskResolution:
    val label: String = "target weights"

  private[paired] final case class TargetStructure(
      design: DMat,
      dual: DualSpectrum
  ) extends TaskResolution:
    val label: String = "target structure"

  private[paired] final case class SourceCovariance(design: DMat) extends TaskResolution:
    val label: String = "source covariance"

  private[paired] final case class SourceRegression(design: DMat) extends TaskResolution:
    val label: String = "source regression"

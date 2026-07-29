package multivar
package family.paired

/** How a fitted resolution becomes a feature-side map for new observations. */
enum TaskResolutionTransfer:
  /** Design lived on the feature / target axis; transfer is feature-specified. */
  case FeatureSpecified

  /** Design lived on training rows; transfer is the induced feature operator from training rows. */
  case TrainingSourceInduced

  def label: String =
    this match
      case FeatureSpecified      => "feature-specified"
      case TrainingSourceInduced => "training-source-induced"

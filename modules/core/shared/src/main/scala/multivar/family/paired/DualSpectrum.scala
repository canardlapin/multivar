package multivar
package family.paired

/** Spectral policy for the dual task operator used by structure resolution. */
enum DualSpectrum:
  /** Moore–Penrose inversion on the numerically supported task range. */
  case MoorePenrose(rankTolerance: Double)
  /** Invert only the first `components` supported task directions. */
  case Truncated(components: Int, rankTolerance: Double)
  /** Replace `1 / d` by `d / (d² + lambda)`. */
  case Ridge(lambda: Double)

  def label: String =
    this match
      case MoorePenrose(_) => "moore-penrose"
      case Truncated(_, _) => "truncated"
      case Ridge(_)        => "ridge"

object DualSpectrum:
  private val DefaultRankTolerance = 1e-12

  /** Moore–Penrose inversion with Multivar's ordinary numerical tolerance. */
  def moorePenrose: DualSpectrum =
    MoorePenrose(DefaultRankTolerance)

  /** Truncation with Multivar's ordinary numerical tolerance. */
  def truncated(components: Int): DualSpectrum =
    Truncated(components, DefaultRankTolerance)

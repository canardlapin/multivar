package multivar
package advanced

import multivar.capability.{FittedCoefficientTransform, FittedFrameTransform}
import multivar.core.{FittedInvertiblePreprocessor, SvdResult}
import multivar.family.canonical.FisherDiscriminantFit
import multivar.family.cpca.{CpcaBlockResult, CpcaFit, PreparedCpcaOperatorFit}
import multivar.family.kernel.{
  KernelEigenArtifact,
  KernelSpec,
  LandmarkSet,
  NystromFit,
  NystromMethod,
  NystromOperatorFit,
  NystromState
}
import multivar.family.paired.{CcaFit, PlscFit, PlsRegressionFit, ReducedRankRegressionFit}
import multivar.family.spectral.{GpcaFit, PcaFit, SvdFit}
import gale.linalg.DMat

/** Expert escape hatches from ordinary fitted results.
  *
  * Import this package when a typed frame, raw SVD payload, or other lifecycle
  * artifact is required. Ordinary analysis code should stay in
  * [[multivar.analysis]].
  */
extension (fit: PcaFit)
  /** The typed functional frame behind this PCA fit. */
  def typedFrame: FittedFrameTransform =
    PcaFit.coreOf(fit).frameTransform

  /** The underlying dense SVD of the preprocessed training matrix. */
  def svdResult: SvdResult =
    PcaFit.coreOf(fit).result

extension (fit: SvdFit)
  /** The typed functional frame behind this SVD fit. */
  def typedFrame: FittedFrameTransform =
    SvdFit.coreOf(fit).frameTransform

  /** The underlying dense SVD of the preprocessed training matrix. */
  def svdResult: SvdResult =
    SvdFit.coreOf(fit).result

extension (fit: GpcaFit)
  def typedFrame: FittedFrameTransform =
    GpcaFit.frameOf(fit)

  def operatorFit =
    GpcaFit.operatorOf(fit)

extension (fit: PlscFit)
  def sourceFrame: FittedFrameTransform = PlscFit.sourceOf(fit)
  def targetFrame: FittedFrameTransform = PlscFit.targetOf(fit)
  def svdResult: SvdResult = PlscFit.resultOf(fit)

extension (fit: CcaFit)
  def sourceFrame: FittedFrameTransform = CcaFit.sourceOf(fit)
  def targetFrame: FittedFrameTransform = CcaFit.targetOf(fit)
  def svdResult: SvdResult = CcaFit.resultOf(fit)

extension (fit: ReducedRankRegressionFit)
  def sourceFrame: FittedFrameTransform = ReducedRankRegressionFit.sourceOf(fit)
  def targetFrame: FittedFrameTransform = ReducedRankRegressionFit.targetOf(fit)
  def coefficientTransform: FittedCoefficientTransform =
    ReducedRankRegressionFit.coefficientTransformOf(fit)

  /** Unconstrained working-space coefficient map retained for diagnostics. */
  def unconstrainedWorkingCoefficients: DMat =
    ReducedRankRegressionFit.unconstrainedWorkingCoefficients(fit)

extension (fit: PlsRegressionFit)
  def sourceFrame: FittedFrameTransform = PlsRegressionFit.frameOf(fit)
  def coefficientTransform: FittedCoefficientTransform =
    PlsRegressionFit.coefficientTransformOf(fit)

extension (fit: FisherDiscriminantFit)
  def typedFrame: FittedFrameTransform =
    FisherDiscriminantFit.frameOf(fit)

extension (fit: CpcaFit)
  def blockResult: CpcaBlockResult =
    CpcaFit.blockOf(fit)

  def operatorFit: PreparedCpcaOperatorFit =
    CpcaFit.operatorOf(fit)

  def fittedPreprocessor: FittedInvertiblePreprocessor =
    CpcaFit.preprocessorOf(fit)

extension (fit: NystromFit)
  def landmarks: LandmarkSet = NystromFit.landmarksOf(fit)
  def operatorFit: NystromOperatorFit = NystromFit.operatorOf(fit)
  def state: NystromState = NystromFit.stateOf(fit)
  def eigen: KernelEigenArtifact = NystromFit.eigenOf(fit)
  def kernelSpec: KernelSpec = NystromFit.kernelSpecOf(fit)
  def method: NystromMethod = NystromFit.methodOf(fit)

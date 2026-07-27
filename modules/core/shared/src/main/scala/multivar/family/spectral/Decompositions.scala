package multivar
package family.spectral


import multivar.core.*
import multivar.capability.*

import gale.linalg.DMat
import gale.linalg.DVec

final case class SvdFit(result: SvdResult, transform: FittedFrameTransform):
  def scores: DMat = transform.scores

  def loadings: DMat = transform.weights

  def singularValues: DVec = result.singularValues

  def project(input: DMat): Either[MultivarError, DMat] =
    project(MatrixView.dense(input))

  def project(input: MatrixView): Either[MultivarError, DMat] =
    transform.project(input)

object Svd:
  def fit(
      input: DMat,
      components: Int
  ): Either[MultivarError, SvdFit] =
    fit(input, components, PreprocessSpec.Pass)

  def fit(
      input: DMat,
      components: Int,
      preproc: PreprocessSpec
  ): Either[MultivarError, SvdFit] =
    ComponentCount(components).flatMap: checked =>
      fit(MatrixView.dense(input), checked, preproc)

  def fit(
      input: MatrixView,
      components: ComponentCount,
      preproc: PreprocessSpec = PreprocessSpec.Pass,
      solver: SvdSolver = DenseSolvers.svd
  ): Either[MultivarError, SvdFit] =
    fitFrame(input, components, preproc, solver, "svd").map((result, transform) => SvdFit(result, transform))

final case class PcaFit(result: SvdResult, transform: FittedFrameTransform):
  /** Training observations expressed in principal-component coordinates. */
  def scores: DMat = transform.scores

  /** Principal axes, with one component per column. */
  def loadings: DMat = transform.weights

  def singularValues: DVec = result.singularValues

  def project(input: DMat): Either[MultivarError, DMat] =
    project(MatrixView.dense(input))

  def project(input: MatrixView): Either[MultivarError, DMat] =
    transform.project(input)

object Pca:
  /** Fit ordinary dense PCA while validating the requested component count.
    *
    * This is the direct entry point for dense data. It expands exactly to the
    * `MatrixView` and `ComponentCount` API below.
    */
  def fit(
      input: DMat,
      components: Int
  ): Either[MultivarError, PcaFit] =
    fit(input, components, PreprocessSpec.Center, DenseSolvers.svd)

  def fit(
      input: DMat,
      components: Int,
      preproc: PreprocessSpec
  ): Either[MultivarError, PcaFit] =
    fit(input, components, preproc, DenseSolvers.svd)

  def fit(
      input: DMat,
      components: Int,
      preproc: PreprocessSpec,
      solver: SvdSolver
  ): Either[MultivarError, PcaFit] =
    ComponentCount(components).flatMap: checked =>
      fit(MatrixView.dense(input), checked, preproc, solver)

  def fit(
      input: MatrixView,
      components: ComponentCount,
      preproc: PreprocessSpec = PreprocessSpec.Center,
      solver: SvdSolver = DenseSolvers.svd
  ): Either[MultivarError, PcaFit] =
    fitFrame(input, components, preproc, solver, "pca").map((result, transform) => PcaFit(result, transform))

private def fitFrame(
    input: MatrixView,
    components: ComponentCount,
    preproc: PreprocessSpec,
    solver: SvdSolver,
    method: String
): Either[MultivarError, (SvdResult, FittedFrameTransform)] =
  for
    fitted <- preproc.fit(input)
    transformed <- fitted.transform(input)
    svd <- solver.decompose(transformed, components)
    _ <- requireComponents(svd)
    transform <- FittedFrameTransform.fromTraining(
      input,
      svd.v,
      fitted,
      method,
      components,
      Some(svd.singularValues)
    )
  yield (svd, transform)

private def requireComponents(svd: SvdResult): Either[MultivarError, Unit] =
  if svd.singularValues.length == 0 then Left(MultivarError.SolverFailed("no singular values above tolerance"))
  else Right(())

package multivar
package family.spectral

import multivar.capability.FittedFrameTransform
import multivar.core.*

import gale.linalg.DMat
import gale.linalg.DVec

/** Dense, dynamically named entry point for generalized PCA.
  *
  * The semantic-diagram API remains the authority when row or feature identity
  * must be carried explicitly. This adapter supplies stable local identities,
  * validates raw component counts, and delegates to the same `GpcaProblem`.
  */
final case class GpcaFit(
    result: GpcaOperatorFit[? <: SemanticSpace, ? <: SemanticSpace, ? <: SemanticSpace],
    transform: FittedFrameTransform
):
  def scores: DMat = transform.scores

  /** Generalized component weights, one component per column. */
  def weights: DMat = transform.weights

  def eigenvalues: DVec = result.generalizedEigenvalues

  def singularValues: DVec = result.singularValues

  def project(input: DMat): Either[MultivarError, DMat] =
    project(MatrixView.dense(input))

  def project(input: MatrixView): Either[MultivarError, DMat] =
    transform.project(input)

object Gpca:
  def fit(
      input: DMat,
      components: Int,
      rowMetric: MetricSpec,
      featureMetric: MetricSpec
  ): Either[MultivarError, GpcaFit] =
    fit(input, components, rowMetric, featureMetric, PreprocessSpec.Center)

  def fit(
      input: DMat,
      components: Int,
      rowMetric: MetricSpec,
      featureMetric: MetricSpec,
      preproc: PreprocessSpec
  ): Either[MultivarError, GpcaFit] =
    for
      checked <- ComponentCount(components)
      _ <-
        if rowMetric.dim == input.rows then Right(())
        else Left(MultivarError.MetricShapeMismatch(IndexAxis.Row, input.rows, rowMetric.dim))
      _ <-
        if featureMetric.dim == input.cols then Right(())
        else Left(MultivarError.MetricShapeMismatch(IndexAxis.Feature, input.cols, featureMetric.dim))
      fitted <- preproc.fit(MatrixView.dense(input))
      prepared <- fitted.transform(MatrixView.dense(input))
      rows <- MvSpace.of("gpca.rows", SpaceRole.Samples, input.rows)
      features <- MvSpace.of("gpca.features", SpaceRole.Observed, input.cols)
      problem <- DynamicGpcaProblem.from(
        prepared,
        rows,
        features,
        rowMetric,
        featureMetric,
        ValueIdentity.source(ValueId.unsafe("gpca.input")),
        SemanticProvenance.source("dense-gpca")
      )
      operator <- problem.fit(checked)
      weights <- gpcaWeights(operator)
      transform <- FittedFrameTransform.fromTraining(
        MatrixView.dense(input),
        weights,
        fitted,
        "gpca",
        checked,
        Some(operator.singularValues)
      )
    yield GpcaFit(operator, transform)

  def fit(
      input: DMat,
      components: Int,
      rowWeights: DVec,
      featureWeights: DVec
  ): Either[MultivarError, GpcaFit] =
    for
      rowMetric <- MetricSpec.diagonal(rowWeights)
      featureMetric <- MetricSpec.diagonal(featureWeights)
      fit <- fit(input, components, rowMetric, featureMetric)
    yield fit

private def gpcaWeights(
    fit: GpcaOperatorFit[? <: SemanticSpace, ? <: SemanticSpace, ? <: SemanticSpace]
): Either[MultivarError, DMat] =
  fit.functionalFrame.weights.toDense.left.map:
    case SemanticError.MultivarFailure(error)  => error
    case SemanticError.LinearMapFailure(error) => LinalgErrorAdapter.toMultivarError(error)
    case error                                 => MultivarError.SolverFailed(error.message)

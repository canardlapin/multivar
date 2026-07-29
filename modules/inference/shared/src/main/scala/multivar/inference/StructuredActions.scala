package multivar.inference

import multivar.core.{RowProjector, RowWhitening}

import gale.linalg.DMat
import resample4s.PermutationTest
import resample4s.kernel.{Permutation, Seed}

extension (permutation: Permutation)
  private[inference] def applyTo(input: DMat): Either[InferenceError, DMat] =
    if input.rows != permutation.domain then
      Left(InferenceError.RowCountMismatch(
        "row permutation input",
        permutation.domain,
        input.rows
      ))
    else
      val out = new Array[Double](input.rows * input.cols)
      val sourceRows = permutation.toIArray
      var targetRow = 0
      while targetRow < input.rows do
        val sourceRow = sourceRows(targetRow)
        var col = 0
        while col < input.cols do
          out(targetRow * input.cols + col) = input(sourceRow, col)
          col += 1
        targetRow += 1
      Right(InferenceNumerics.matrixFromRows(Vector.tabulate(input.rows) { row =>
        Vector.tabulate(input.cols)(col => out(row * input.cols + col))
      }))

sealed trait PermutationAction:
  def rowCount: RowCount
  def draw(seed: Seed, replicate: ReplicateId): Either[InferenceError, Permutation]

object PermutationAction:
  final case class Unrestricted private[inference] (
      rowCount: RowCount
  ) extends PermutationAction:
    override def draw(
        seed: Seed,
        replicate: ReplicateId
    ): Either[InferenceError, Permutation] =
      ResamplingPlans.permutation(seed, replicate, rowCount.value)

  final case class WithinGroups private[inference] (
      partition: RowPartition
  ) extends PermutationAction:
    override def rowCount: RowCount = partition.rowCount

    override def draw(
        seed: Seed,
        replicate: ReplicateId
    ): Either[InferenceError, Permutation] =
      for
        replicateSeed <- ResamplingPlans.replicateSeed(seed, replicate)
        plan <- adapt(PermutationTest
          .within(partition.blocks, 1)
          .plan(rowCount.value, replicateSeed.value))
      yield plan.compiled.plan.first

  final case class WholeClusters private[inference] (
      partition: ClusterPartition
  ) extends PermutationAction:
    override def rowCount: RowCount = partition.rowCount

    override def draw(
        seed: Seed,
        replicate: ReplicateId
    ): Either[InferenceError, Permutation] =
      for
        replicateSeed <- ResamplingPlans.replicateSeed(seed, replicate)
        plan <- adapt(PermutationTest
          .wholeGroups(partition.groups, 1)
          .plan(rowCount.value, replicateSeed.value))
      yield plan.compiled.plan.first

  def unrestricted(rows: RowCount): PermutationAction =
    Unrestricted(rows)

  def withinBlocks(partition: RowPartition): PermutationAction =
    WithinGroups(partition)

  def withinStrata(partition: StrataPartition): PermutationAction =
    WithinGroups(partition.value)

  def wholeClusters(partition: ClusterPartition): Either[InferenceError, PermutationAction] =
    val sizes = partition.clusters.map(_.length).distinct
    if sizes.length != 1 then
      Left(InferenceError.UnsupportedProblem(
        "whole-cluster permutation requires equal cluster sizes"
      ))
    else Right(WholeClusters(partition))

  def forDesign[K <: DesignKind](
      design: ResamplingDesign[K]
  ): Either[InferenceError, PermutationAction] =
    design.conditioning match
      case Conditioning.Unadjusted =>
        forExchangeability(design)
      case Conditioning.Nuisance(_, _, _) =>
        Left(InferenceError.UnsupportedProblem(
          "conditioned permutation requires prepared projection and whitening resources"
        ))

  private[inference] def forExchangeability(
      design: ResamplingDesign[?]
  ): Either[InferenceError, PermutationAction] =
    (design.units, design.exchangeability) match
      case (SamplingUnits.Rows(rows), Exchangeability.Unrestricted) =>
        Right(unrestricted(rows))
      case (SamplingUnits.Rows(_), Exchangeability.WithinBlocks(partition)) =>
        Right(withinBlocks(partition))
      case (SamplingUnits.Rows(_), Exchangeability.WithinStrata(partition)) =>
        Right(withinStrata(partition))
      case (SamplingUnits.Clusters(partition), Exchangeability.Unrestricted) =>
        wholeClusters(partition)
      case _ =>
        Left(InferenceError.UnsupportedProblem(
          "cluster permutation within row-level blocks or strata has no proved action"
        ))

  private def adapt[A](
      value: Either[resample4s.spi.DesignError, A]
  ): Either[InferenceError, A] =
    value.left.map(InferenceError.Resampling.apply)

private[inference] object ProtocolDesign:
  def permutation(
      design: ResamplingDesign[?],
      expectedRows: Int,
      role: String
  ): Either[InferenceError, PermutationAction] =
    for
      action <- PermutationAction.forDesign(design)
      _ <-
        if action.rowCount.value == expectedRows then Right(())
        else Left(InferenceError.RowCountMismatch(
          role,
          expectedRows,
          action.rowCount.value
        ))
    yield action

final case class ResidualPermutationAction private (
    reference: ConditioningRef,
    rowCount: RowCount,
    nuisance: RowProjector,
    residual: RowProjector,
    whitening: RowWhitening,
    permutation: PermutationAction
):
  def draw(
      input: DMat,
      seed: Seed,
      replicate: ReplicateId
  ): Either[InferenceError, DMat] =
    if input.rows != rowCount.value then
      Left(InferenceError.RowCountMismatch(
        "conditioned permutation input",
        rowCount.value,
        input.rows
      ))
    else
      firstNonFinite(input) match
        case Some((index, value)) =>
          Left(InferenceError.NonFiniteStatistic(
            s"conditioned permutation input entry $index",
            value
          ))
        case None =>
          for
            whitened <- adapt("row whitening", whitening.whiten(input))
            fitted <- adapt("nuisance projection", nuisance.project(whitened))
            residuals <- adapt("residual projection", residual.project(whitened))
            rowOrder <- permutation.draw(seed, replicate)
            permuted <- rowOrder.applyTo(residuals)
            randomized = add(fitted, permuted)
            restored <- adapt("row unwhitening", whitening.unwhiten(randomized))
          yield restored

  private def firstNonFinite(input: DMat): Option[(Int, Double)] =
    val values = input.copyData
    var i = 0
    while i < values.length do
      if !values(i).isFinite then return Some((i, values(i)))
      i += 1
    None

  private def add(left: DMat, right: DMat): DMat =
    val out = new Array[Double](left.rows * left.cols)
    var row = 0
    while row < left.rows do
      var col = 0
      while col < left.cols do
        out(row * left.cols + col) = left(row, col) + right(row, col)
        col += 1
      row += 1
    InferenceNumerics.matrixFromRows(Vector.tabulate(left.rows) { row =>
      Vector.tabulate(left.cols)(col => out(row * left.cols + col))
    })

  private def adapt[A](
      role: String,
      value: Either[multivar.core.MultivarError, A]
  ): Either[InferenceError, A] =
    value.left.map(error => InferenceError.NumericalFailure(role, error.message))

object ResidualPermutationAction:
  def from(
      design: ResamplingDesign[DesignKind.ConditionedRows],
      reference: ConditioningRef,
      nuisance: RowProjector,
      whitening: Option[RowWhitening]
  ): Either[InferenceError, ResidualPermutationAction] =
    design.conditioning match
      case Conditioning.Unadjusted =>
        Left(InferenceError.UnsupportedProblem(
          "residual permutation requires a nuisance-conditioned design"
        ))
      case Conditioning.Nuisance(expectedReference, rows, requirement) =>
        if reference != expectedReference then
          Left(InferenceError.UnsupportedProblem(
            s"conditioning resource '${reference.value}' does not match '${expectedReference.value}'"
          ))
        else if nuisance.rows != rows.value then
          Left(InferenceError.RowCountMismatch(
            "nuisance projector",
            rows.value,
            nuisance.rows
          ))
        else if requirement == WhiteningRequirement.Required && whitening.isEmpty then
          Left(InferenceError.UnsupportedProblem(
            "conditioned design requires an explicit row-whitening capability"
          ))
        else
          for
            permutation <- PermutationAction.forExchangeability(design)
            checkedWhitening <- whitening match
              case Some(value) if value.rows != rows.value =>
                Left(InferenceError.RowCountMismatch(
                  "row whitening",
                  rows.value,
                  value.rows
                ))
              case Some(value) => Right(value)
              case None => adaptIdentity(RowWhitening.identity(rows.value))
          yield ResidualPermutationAction(
            reference,
            rows,
            nuisance,
            nuisance.complement,
            checkedWhitening,
            permutation
          )

  private def adaptIdentity(
      value: Either[multivar.core.MultivarError, RowWhitening]
  ): Either[InferenceError, RowWhitening] =
    value.left.map(error => InferenceError.NumericalFailure("identity row whitening", error.message))

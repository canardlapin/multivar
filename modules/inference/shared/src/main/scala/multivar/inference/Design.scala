package multivar.inference

import resample4s.{Blocks, Groups, Strata}
import resample4s.kernel.{LabelRefinement, Labels}

final case class RowPartition private (
    rowCount: RowCount,
    labels: Labels
):
  lazy val groups: Vector[Vector[RowIx]] =
    val out = Array.fill(labels.cardinality)(Vector.newBuilder[RowIx])
    val codes = labels.toIArray
    var row = 0
    while row < codes.length do
      out(codes(row)) += RowIx.unsafe(row)
      row += 1
    out.iterator.map(_.result()).toVector

  def blocks: Blocks = Blocks.from(labels)

  def groupIndexByRow: Array[Int] =
    val canonical = labels.toIArray
    val out = new Array[Int](canonical.length)
    var index = 0
    while index < canonical.length do
      out(index) = canonical(index)
      index += 1
    out

object RowPartition:
  def fromLabels(
      rows: RowCount,
      labels: Labels
  ): Either[InferenceError, RowPartition] =
    if labels.size != rows.value then
      Left(InferenceError.RowCountMismatch(
        "row partition labels",
        rows.value,
        labels.size
      ))
    else Right(RowPartition(rows, labels))

  def from(rows: RowCount, groups: Iterable[Iterable[Int]]): Either[InferenceError, RowPartition] =
    val raw = groups.iterator.map(_.toVector).toVector
    if raw.isEmpty then Left(InferenceError.InvalidPartition("at least one group is required"))
    else if raw.exists(_.isEmpty) then Left(InferenceError.InvalidPartition("groups must be non-empty"))
    else
      val seen = Array.fill(rows.value)(false)
      val codes = new Array[Int](rows.value)
      var groupIndex = 0
      var error = Option.empty[InferenceError]
      while groupIndex < raw.length && error.isEmpty do
        val group = raw(groupIndex)
        var i = 0
        while i < group.length && error.isEmpty do
          val row = group(i)
          if row < 0 || row >= rows.value then
            error = Some(InferenceError.InvalidPartition(s"row $row is outside [0, ${rows.value})"))
          else if seen(row) then
            error = Some(InferenceError.InvalidPartition(s"row $row appears more than once"))
          else
            seen(row) = true
            codes(row) = groupIndex
          i += 1
        groupIndex += 1
      var row = 0
      while row < rows.value && error.isEmpty do
        if !seen(row) then error = Some(InferenceError.InvalidPartition(s"row $row is not covered"))
        row += 1
      error match
        case Some(value) => Left(value)
        case None =>
          Labels
            .dense(IArray.unsafeFromArray(codes))
            .left
            .map(InferenceError.Resampling.apply)
            .map(RowPartition(rows, _))

final case class ClusterPartition private (value: RowPartition):
  def rowCount: RowCount = value.rowCount
  def clusters: Vector[Vector[RowIx]] = value.groups
  def groups: Groups = Groups.from(value.labels)

object ClusterPartition:
  def fromLabels(
      rows: RowCount,
      groups: Groups
  ): Either[InferenceError, ClusterPartition] =
    RowPartition.fromLabels(rows, groups.labels).map(ClusterPartition(_))

  def from(rows: RowCount, clusters: Iterable[Iterable[Int]]): Either[InferenceError, ClusterPartition] =
    RowPartition.from(rows, clusters).map(ClusterPartition(_))

final case class StrataPartition private (value: RowPartition):
  def rowCount: RowCount = value.rowCount
  def strata: Vector[Vector[RowIx]] = value.groups
  def labels: Strata = Strata.from(value.labels)

object StrataPartition:
  def fromLabels(
      rows: RowCount,
      strata: Strata
  ): Either[InferenceError, StrataPartition] =
    RowPartition.fromLabels(rows, strata.labels).map(StrataPartition(_))

  def from(rows: RowCount, strata: Iterable[Iterable[Int]]): Either[InferenceError, StrataPartition] =
    RowPartition.from(rows, strata).map(StrataPartition(_))

final case class NestedPartitions private (
    inner: RowPartition,
    outer: RowPartition,
    refinement: LabelRefinement
)

object NestedPartitions:
  def from(
      inner: RowPartition,
      outer: RowPartition
  ): Either[InferenceError, NestedPartitions] =
    if inner.rowCount != outer.rowCount then
      Left(InferenceError.RowCountMismatch(
        "nested outer partition",
        inner.rowCount.value,
        outer.rowCount.value
      ))
    else
      inner.labels
        .refines(outer.labels)
        .left
        .map(InferenceError.Resampling.apply)
        .map(NestedPartitions(inner, outer, _))

enum WhiteningRequirement:
  case NotRequired
  case Required

enum SamplingUnits:
  case Rows(count: RowCount)
  case Clusters(partition: ClusterPartition)

  def rowCount: RowCount =
    this match
      case Rows(count)         => count
      case Clusters(partition) => partition.rowCount

enum Exchangeability:
  case Unrestricted
  case WithinBlocks(partition: RowPartition)
  case WithinStrata(partition: StrataPartition)

enum Conditioning:
  case Unadjusted
  case Nuisance(
      reference: ConditioningRef,
      rows: RowCount,
      whitening: WhiteningRequirement
  )

sealed trait DesignKind

object DesignKind:
  sealed trait ExchangeableRows extends DesignKind
  sealed trait WithinBlockRows extends DesignKind
  sealed trait WithinStrataRows extends DesignKind
  sealed trait ClusterRows extends DesignKind
  sealed trait ConditionedRows extends DesignKind

final case class ResamplingDesign[K <: DesignKind] private (
    units: SamplingUnits,
    exchangeability: Exchangeability,
    conditioning: Conditioning
):
  def rowCount: RowCount = units.rowCount

object ResamplingDesign:
  def exchangeableRows(rows: RowCount): ResamplingDesign[DesignKind.ExchangeableRows] =
    ResamplingDesign(SamplingUnits.Rows(rows), Exchangeability.Unrestricted, Conditioning.Unadjusted)

  def withinBlocks(partition: RowPartition): ResamplingDesign[DesignKind.WithinBlockRows] =
    ResamplingDesign(
      SamplingUnits.Rows(partition.rowCount),
      Exchangeability.WithinBlocks(partition),
      Conditioning.Unadjusted
    )

  def withinStrata(partition: StrataPartition): ResamplingDesign[DesignKind.WithinStrataRows] =
    ResamplingDesign(
      SamplingUnits.Rows(partition.rowCount),
      Exchangeability.WithinStrata(partition),
      Conditioning.Unadjusted
    )

  def clustered(partition: ClusterPartition): ResamplingDesign[DesignKind.ClusterRows] =
    ResamplingDesign(SamplingUnits.Clusters(partition), Exchangeability.Unrestricted, Conditioning.Unadjusted)

  def nuisanceAdjusted(
      rows: RowCount,
      reference: ConditioningRef,
      whitening: WhiteningRequirement
  ): ResamplingDesign[DesignKind.ConditionedRows] =
    ResamplingDesign(
      SamplingUnits.Rows(rows),
      Exchangeability.Unrestricted,
      Conditioning.Nuisance(reference, rows, whitening)
    )

  def nuisanceAdjustedWithinBlocks(
      partition: RowPartition,
      reference: ConditioningRef,
      whitening: WhiteningRequirement
  ): ResamplingDesign[DesignKind.ConditionedRows] =
    ResamplingDesign(
      SamplingUnits.Rows(partition.rowCount),
      Exchangeability.WithinBlocks(partition),
      Conditioning.Nuisance(reference, partition.rowCount, whitening)
    )

  def nuisanceAdjustedWithinStrata(
      partition: StrataPartition,
      reference: ConditioningRef,
      whitening: WhiteningRequirement
  ): ResamplingDesign[DesignKind.ConditionedRows] =
    ResamplingDesign(
      SamplingUnits.Rows(partition.rowCount),
      Exchangeability.WithinStrata(partition),
      Conditioning.Nuisance(reference, partition.rowCount, whitening)
    )

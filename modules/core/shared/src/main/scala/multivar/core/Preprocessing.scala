package multivar
package core

import gale.linalg.DMat
import gale.linalg.DVec

/** Denominator convention for the standard deviations used to standardize columns.
  *
  * The convention changes the fitted scale, and therefore the fitted model, so it is
  * part of the request rather than an implementation detail.
  */
enum VarianceConvention:
  /** Divide the centered sum of squares by `n - 1`. Requires at least two rows. */
  case Sample

  /** Divide the centered sum of squares by `n`. */
  case Population

  def label: String =
    this match
      case Sample     => "sample"
      case Population => "population"

  /** Rows required before the convention is defined. */
  def minimumRows: Int =
    this match
      case Sample     => 2
      case Population => 1

  def denominator(count: Int): Double =
    this match
      case Sample     => (count - 1).toDouble
      case Population => count.toDouble

enum PreprocessSpec:
  case Pass
  case Center

  /** Multiply each column by the corresponding weight.
    *
    * Weights must be finite but may be negative or zero. A zero weight collapses its
    * column, which is a valid transform but not an invertible one, so a fitted
    * preprocessor built from one cannot satisfy [[FittedPreprocessor.requireInvertible]].
    */
  case MultiplyColumns(weights: DVec)

  /** Standardize columns to zero mean and unit variance under `convention`.
    *
    * Degenerate (numerically constant) columns — those whose standard deviation
    * is at most `MatrixView.DegenerateScaleEpsilon` relative to the magnitude of the
    * column mean — are centered but not rescaled (their scale weight is 1), matching
    * scikit-learn's `StandardScaler` convention for zero-variance features. The test
    * is relative to the column's own magnitude, so genuine variation in tiny-valued
    * columns is still standardized while rounding noise on huge-valued constant
    * columns is not mistaken for signal.
    */
  case Standardize(convention: VarianceConvention = VarianceConvention.Sample)

  def fit(input: MatrixView): Either[MultivarError, FittedPreprocessor] =
    if input.cols <= 0 then Left(MultivarError.InvalidDimension("preprocessing input columns", input.cols))
    else
      this match
        case Pass =>
          Right(FittedColumnAffine(input.cols, MatrixView.ones(input.cols), MatrixView.zeros(input.cols)))
        case Center =>
          input.columnStats.flatMap(_.means).map { means =>
            FittedColumnAffine(input.cols, MatrixView.ones(input.cols), MatrixView.negate(means))
          }
        case MultiplyColumns(weights) =>
          for
            _ <- MatrixView.requireVectorLength("preprocessing weights", weights, input.cols)
            _ <- MatrixView.requireFinite("preprocessing weights", weights)
          yield FittedColumnAffine(input.cols, weights, MatrixView.zeros(input.cols))
        case Standardize(convention) =>
          for
            stats <- input.columnStats
            means <- stats.means
            sds <- stats.standardDeviations(convention)
            scale <- MatrixView.invert(PreprocessSpec.safeStandardizeScale(means, sds))
          yield FittedColumnAffine(input.cols, scale, MatrixView.multiply(MatrixView.negate(means), scale))

  /** Fit, and prove at the same time that the fitted transform can be undone.
    *
    * Response-side preprocessing in a regression must be invertible, because a
    * prediction is only meaningful in the original response units. Requiring the proof
    * here turns "fitting succeeded and prediction failed later" into a single failure
    * at the point where the offending scale was supplied.
    */
  def fitInvertible(input: MatrixView): Either[MultivarError, FittedInvertiblePreprocessor] =
    fit(input).flatMap(_.requireInvertible)

object PreprocessSpec:
  def multiplyColumns(weights: Seq[Double]): Either[MultivarError, PreprocessSpec] =
    val vector = DVec.fromSeq(weights)
    MatrixView.requireFinite("preprocessing weights", vector).map(_ => PreprocessSpec.MultiplyColumns(vector))

  /** Column scale for standardization: the sample standard deviation, or 1.0 for
    * degenerate columns whose spread is negligible relative to their mean magnitude
    * (see the `Standardize` case documentation).
    */
  private def safeStandardizeScale(means: DVec, sds: DVec): DVec =
    val out = new Array[Double](sds.length)
    var col = 0
    while col < sds.length do
      val sd = sds(col)
      out(col) =
        if sd > MatrixView.DegenerateScaleEpsilon * Math.abs(means(col)) then sd
        else 1.0
      col += 1
    GaleNumerics.vectorFromArray(out)

/** Fitted column preprocessing that can be applied to new observations.
  *
  * Applying preprocessing is always possible; undoing it is not. A zero column scale
  * transforms perfectly well and cannot be reversed. That distinction lives in the
  * type: only [[FittedInvertiblePreprocessor]] can restore original coordinates, and
  * [[requireInvertible]] is the one way to obtain it.
  */
trait FittedPreprocessor:
  def inputCols: Int

  def transform(
      input: MatrixView,
      columns: Option[IndexSet] = None,
      policy: StoragePolicy = StoragePolicy.Operator
  ): Either[MultivarError, MatrixView]

  def restrict(columns: IndexSet): Either[MultivarError, FittedPreprocessor]

  /** Prove that this fitted preprocessing can be undone, or explain why it cannot.
    *
    * Callers that need original coordinates should demand the invertible type in their
    * own signatures and require the proof once, rather than discovering a
    * non-invertible scale on every reconstruction.
    */
  def requireInvertible: Either[MultivarError, FittedInvertiblePreprocessor]

/** Fitted preprocessing whose inverse was computed and checked at construction, so
  * `inverseTransform` fails only on a shape mismatch, never on a zero scale.
  */
trait FittedInvertiblePreprocessor extends FittedPreprocessor:
  def inverseTransform(
      input: MatrixView,
      columns: Option[IndexSet] = None,
      policy: StoragePolicy = StoragePolicy.Operator
  ): Either[MultivarError, MatrixView]

  /** Closed-form dense inverse when the fitted map is column affine; otherwise
    * falls back to [[inverseTransform]] followed by densification.
    */
  private[multivar] def inverseTransformDense(
      working: DMat,
      columns: Option[IndexSet] = None
  ): Either[MultivarError, DMat] =
    inverseTransform(MatrixView.dense(working), columns, StoragePolicy.AllowDense)
      .flatMap(_.toDense(StoragePolicy.AllowDense))

  /** Map a processed contribution back to original coordinates, cancelling any
    * fitted affine shift. Column affines use a single scaled pass; other
    * invertible preprocessors retain the difference-of-inverses definition.
    */
  private[multivar] def inverseContributionDense(processed: DMat): Either[MultivarError, DMat] =
    val zero = DMat.zeros(processed.rows, processed.cols)
    for
      original <- inverseTransformDense(processed)
      originalZero <- inverseTransformDense(zero)
    yield MatrixOps.subtract(original, originalZero)

  override def requireInvertible: Either[MultivarError, FittedInvertiblePreprocessor] =
    Right(this)

/** Fitted centering and scaling in the familiar `(x - center) / scale` form.
  *
  * The stored form is a multiplier and an offset, which composes better but reads
  * less like the statistics it represents. This is the reporting form.
  */
final case class ColumnAffineSummary(center: DVec, scale: DVec):
  require(center.length == scale.length, "affine summary center and scale must have equal length")

object ColumnAffineSummary:
  /** Describe fitted preprocessing as centering and scaling, when it can be.
    *
    * Preprocessing that collapses a column has no finite scale to report, and
    * preprocessing that is not a column affine has no such description at all.
    */
  def of(preprocessor: FittedPreprocessor): Option[ColumnAffineSummary] =
    preprocessor.requireInvertible.toOption.collect { case affine: InvertibleColumnAffine =>
      affine.summary
    }

final case class FittedColumnAffine(
    inputCols: Int,
    scale: DVec,
    shift: DVec
) extends FittedPreprocessor:
  require(inputCols > 0, "fitted preprocessor input columns must be positive")
  require(scale.length == inputCols, "scale length must match input columns")
  require(shift.length == inputCols, "shift length must match input columns")

  override def transform(
      input: MatrixView,
      columns: Option[IndexSet],
      policy: StoragePolicy
  ): Either[MultivarError, MatrixView] =
    parametersFor(input, columns).flatMap { case (selectedScale, selectedShift) =>
      MatrixView.affine(input, selectedScale, selectedShift, policy, "preprocessing transform")
    }

  /** The inverse map, which is itself an affine, when every column scale is invertible.
    *
    * Undoing `x * scale + shift` is applying `x * (1/scale) + (-shift/scale)`, so the
    * inverse needs no separate machinery and column selection behaves identically.
    */
  override def requireInvertible: Either[MultivarError, FittedInvertiblePreprocessor] =
    MatrixView.invert(scale).map { inverseScale =>
      InvertibleColumnAffine(
        this,
        FittedColumnAffine(
          inputCols,
          inverseScale,
          MatrixView.multiply(MatrixView.negate(shift), inverseScale)
        )
      )
    }

  override def restrict(columns: IndexSet): Either[MultivarError, FittedPreprocessor] =
    MatrixView.requireColumnIndexSet(columns, inputCols).map { checked =>
      FittedColumnAffine(
        inputCols = checked.length,
        scale = MatrixView.selectVector(scale, checked),
        shift = MatrixView.selectVector(shift, checked)
      )
    }

  private def parametersFor(
      input: MatrixView,
      columns: Option[IndexSet]
  ): Either[MultivarError, (DVec, DVec)] =
    columns match
      case None =>
        if input.cols != inputCols then
          Left(MultivarError.MatrixShapeMismatch(s"input has ${input.cols} columns but preprocessor expects $inputCols"))
        else Right((scale, shift))
      case Some(indices) =>
        MatrixView.requireColumnIndexSet(indices, inputCols).flatMap { checked =>
          if input.cols != checked.length then
            Left(
              MultivarError.MatrixShapeMismatch(
                s"input has ${input.cols} columns but column selection has ${checked.length}"
              )
            )
          else Right((MatrixView.selectVector(scale, checked), MatrixView.selectVector(shift, checked)))
        }

/** A [[FittedColumnAffine]] paired with the affine that undoes it.
  *
  * Holding the inverse rather than recomputing it is what makes the invertibility a
  * property of the value: it could not have been constructed if any column scale were
  * zero or non-finite.
  */
final class InvertibleColumnAffine private[core] (
    val forward: FittedColumnAffine,
    private val inverse: FittedColumnAffine
) extends FittedInvertiblePreprocessor:
  override def inputCols: Int =
    forward.inputCols

  override def transform(
      input: MatrixView,
      columns: Option[IndexSet],
      policy: StoragePolicy
  ): Either[MultivarError, MatrixView] =
    forward.transform(input, columns, policy)

  override def inverseTransform(
      input: MatrixView,
      columns: Option[IndexSet],
      policy: StoragePolicy
  ): Either[MultivarError, MatrixView] =
    inverse.transform(input, columns, policy)

  override private[multivar] def inverseTransformDense(
      working: DMat,
      columns: Option[IndexSet]
  ): Either[MultivarError, DMat] =
    columns match
      case None =>
        if working.cols != inverse.inputCols then
          Left(
            MultivarError.MatrixShapeMismatch(
              s"input has ${working.cols} columns but preprocessor expects ${inverse.inputCols}"
            )
          )
        else Right(MatrixView.materializeAffine(working, inverse.scale, inverse.shift))
      case Some(indices) =>
        for
          checked <- MatrixView.requireColumnIndexSet(indices, inputCols)
          _ <-
            if working.cols != checked.length then
              Left(
                MultivarError.MatrixShapeMismatch(
                  s"input has ${working.cols} columns but column selection has ${checked.length}"
                )
              )
            else Right(())
        yield
          MatrixView.materializeAffine(
            working,
            MatrixView.selectVector(inverse.scale, checked),
            MatrixView.selectVector(inverse.shift, checked)
          )

  override private[multivar] def inverseContributionDense(processed: DMat): Either[MultivarError, DMat] =
    if processed.cols != inverse.inputCols then
      Left(
        MultivarError.MatrixShapeMismatch(
          s"input has ${processed.cols} columns but preprocessor expects ${inverse.inputCols}"
        )
      )
    else
      Right(
        MatrixView.materializeAffine(processed, inverse.scale, MatrixView.zeros(inverse.scale.length))
      )

  override def restrict(columns: IndexSet): Either[MultivarError, FittedPreprocessor] =
    forward.restrict(columns)

  /** Fitted centering and scaling as `(x - center) / scale`. */
  def summary: ColumnAffineSummary =
    ColumnAffineSummary(
      center = MatrixView.negate(MatrixView.multiply(forward.shift, inverse.scale)),
      scale = inverse.scale
    )

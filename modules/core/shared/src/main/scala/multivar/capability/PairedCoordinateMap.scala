package multivar
package capability

import multivar.core.*

import gale.linalg.DMat
import gale.linalg.DVec

/** Working-space coefficient map decoded into original predictor/response units.
  *
  * With forward preprocessing `x ↦ x ⊙ D + a`, the identity
  * `y = ((x ⊙ D_x + a_x) B_w − a_y) ⊙ D_y^{-1}` rearranges to
  * `y = x B_raw + intercept` with `B_raw = D_x B_w D_y^{-1}` and
  * `intercept = (a_x B_w − a_y) ⊙ D_y^{-1}`.
  */
private[multivar] final case class RawCoefficientMap(
    coefficients: DMat,
    intercept: DVec
)

private[multivar] object PairedCoordinateMap:
  def decode(
      working: DMat,
      predictors: FittedPreprocessor,
      responses: FittedInvertiblePreprocessor
  ): Either[MultivarError, RawCoefficientMap] =
    for
      predictorAffine <- columnAffine(predictors, "predictor")
      responseAffine <- responses match
        case affine: InvertibleColumnAffine => Right(affine)
        case other =>
          Left(
            MultivarError.InvalidMap(
              s"response preprocessor must be a column affine to expose raw coefficients, got ${other.getClass.getName}"
            )
          )
    yield
      val dx = predictorAffine.scale
      val ax = predictorAffine.shift
      val dyInv = responseAffine.summary.scale
      val ay = responseAffine.forward.shift
      val rowScaled = MatrixView.scaleRows(working, dx)
      val raw = MetricOperator.scaleColumnsDense(rowScaled, dyInv)
      val intercept = interceptFrom(ax, working, ay, dyInv)
      RawCoefficientMap(raw, intercept)

  private def columnAffine(
      preprocessor: FittedPreprocessor,
      role: String
  ): Either[MultivarError, FittedColumnAffine] =
    preprocessor match
      case affine: FittedColumnAffine     => Right(affine)
      case affine: InvertibleColumnAffine => Right(affine.forward)
      case other =>
        Left(
          MultivarError.InvalidMap(
            s"$role preprocessor must be a column affine to expose raw coefficients, got ${other.getClass.getName}"
          )
        )

  private def interceptFrom(ax: DVec, working: DMat, ay: DVec, dyInv: DVec): DVec =
    val out = new Array[Double](working.cols)
    var col = 0
    while col < working.cols do
      var sum = 0.0
      var row = 0
      while row < working.rows do
        sum += ax(row) * working(row, col)
        row += 1
      out(col) = (sum - ay(col)) * dyInv(col)
      col += 1
    GaleNumerics.vectorFromArray(out)

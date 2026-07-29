package multivar.inference

import multivar.core.{ComponentCount, MatrixView}
import multivar.family.paired.{ReducedRankRegression, RegressionRegularization}

import gale.linalg.DMat
import resample4s.kernel.{CodomainMismatch, IndexSpace, Permutation, Seed, Selection, Split}
import resample4s.spi.DesignError

final case class HeldOutSplit private (
    rowCount: RowCount,
    split: Split[Selection]
):
  def training: Selection = split.analysis
  def test: Selection = split.assessment

object HeldOutSplit:
  def fromSplit(
      split: Split[Selection]
  ): Either[InferenceError, HeldOutSplit] =
    if split.assessment.domain == 0 then
      Left(InferenceError.InvalidPartition("test split must be non-empty"))
    else if split.analysis.domain + split.assessment.domain != split.analysis.codomain then
      Left(InferenceError.InvalidPartition(
        "held-out split must assign every population row exactly once"
      ))
    else
      RowCount(split.analysis.codomain).map(HeldOutSplit(_, split))

  def from(
      rows: RowCount,
      training: Iterable[Int],
      test: Iterable[Int]
  ): Either[InferenceError, HeldOutSplit] =
    val train = training.iterator.toVector
    val heldOut = test.iterator.toVector
    if train.isEmpty then Left(InferenceError.InvalidPartition("training split must be non-empty"))
    else if heldOut.isEmpty then Left(InferenceError.InvalidPartition("test split must be non-empty"))
    else
      val seen = Array.fill(rows.value)(false)
      def validate(values: Vector[Int], role: String): Option[InferenceError] =
        var i = 0
        while i < values.length do
          val row = values(i)
          if row < 0 || row >= rows.value then
            return Some(InferenceError.InvalidPartition(
              s"$role row $row is outside [0, ${rows.value})"
            ))
          else if seen(row) then
            return Some(InferenceError.InvalidPartition(
              s"row $row occurs more than once across the held-out split"
            ))
          seen(row) = true
          i += 1
        None
      validate(train, "training")
        .orElse(validate(heldOut, "test")) match
        case Some(error) => Left(error)
        case None =>
          var row = 0
          while row < rows.value do
            if !seen(row) then
              return Left(InferenceError.InvalidPartition(
                s"row $row is not assigned to training or test"
              ))
            row += 1
          for
            space <- indexSpace(rows.value)
            trainingSelection <- adapt(Selection.from(
              IArray.from(train.sorted),
              space
            ))
            testSelection <- adapt(Selection.from(
              IArray.from(heldOut.sorted),
              space
            ))
            split <- Split.of(trainingSelection, testSelection) match
              case Right(value) => Right(value)
              case Left(error: DesignError) =>
                Left(InferenceError.Resampling(error))
              case Left(error: CodomainMismatch) =>
                Left(InferenceError.InvalidPartition(
                  s"training/test codomains differ: ${error.left} and ${error.right}"
                ))
            heldOutSplit <- fromSplit(split)
          yield heldOutSplit

  private def adapt[A](
      value: Either[DesignError, A]
  ): Either[InferenceError, A] =
    value.left.map(InferenceError.Resampling.apply)

  private def indexSpace(size: Int): Either[InferenceError, IndexSpace] =
    IndexSpace.of(size).left.map(InferenceError.Resampling.apply)

final case class PredictiveData private (
    x: DMat,
    y: DMat
)

object PredictiveData:
  def from(x: DMat, y: DMat): Either[InferenceError, PredictiveData] =
    if x.rows != y.rows then
      Left(InferenceError.RowCountMismatch("predictive X/Y", x.rows, y.rows))
    else if x.rows < 3 then Left(InferenceError.InvalidCount("predictive rows", x.rows))
    else if x.cols < 1 then Left(InferenceError.InvalidCount("predictive X columns", x.cols))
    else if y.cols < 1 then Left(InferenceError.InvalidCount("predictive Y columns", y.cols))
    else
      FamilyPredictionMatrices.validateFinite("predictive X", x)
        .flatMap(_ => FamilyPredictionMatrices.validateFinite("predictive Y", y))
        .map(_ => PredictiveData(x, y))

final case class PredictiveGainResult(
    gain: Double,
    modelSquaredError: Double,
    baselineSquaredError: Double,
    split: HeldOutSplit,
    validity: ValidityClaim
)

final case class RrrPredictiveProtocol(
    components: ComponentCount,
    regularization: RegressionRegularization = RegressionRegularization.Ols
):
  val target: TargetSpec[TargetKind.PredictiveGain] = TargetSpec.PredictiveGain
  val nullHypothesis: NullSpec[NullKind.ResponsePermutation] = NullSpec.PermuteResponse
  val validity: ValidityClaim = ValidityClaim.Conditional

  def observe(
      data: PredictiveData,
      split: HeldOutSplit
  ): Either[InferenceError, PredictiveGainResult] =
    validateSplit(data, split).flatMap(_ => evaluate(data, split, None))

  def nullStatistic(
      data: PredictiveData,
      split: HeldOutSplit,
      seed: Seed,
      replicate: ReplicateId
  ): Either[InferenceError, Double] =
    for
      _ <- validateSplit(data, split)
      permutation <- ResamplingPlans.permutation(
        seed,
        replicate,
        split.training.domain
      )
      result <- evaluate(data, split, Some(permutation))
    yield result.gain

  private def validateSplit(
      data: PredictiveData,
      split: HeldOutSplit
  ): Either[InferenceError, Unit] =
    if split.rowCount.value == data.x.rows then Right(())
    else Left(InferenceError.RowCountMismatch(
      "predictive held-out split",
      data.x.rows,
      split.rowCount.value
    ))

  private def evaluate(
      data: PredictiveData,
      split: HeldOutSplit,
      trainingResponsePermutation: Option[Permutation]
  ): Either[InferenceError, PredictiveGainResult] =
    val trainRows = split.training.toVector
    val testRows = split.test.toVector
    val xTrain = data.x.selectRows(trainRows)
    val xTest = data.x.selectRows(testRows)
    val yTrainRaw = data.y.selectRows(trainRows)
    val yTrain = trainingResponsePermutation match
      case Some(permutation) => yTrainRaw.selectRows(permutation.toVector)
      case None              => yTrainRaw
    val yTest = data.y.selectRows(testRows)

    ReducedRankRegression
      .fit(
        MatrixView.dense(xTrain),
        MatrixView.dense(yTrain),
        components,
        regularization
      )
      .left
      .map(error =>
        InferenceError.NumericalFailure(
          "RRR predictive fit",
          error.message
        )
      )
      .flatMap(fit =>
        fit.predict(MatrixView.dense(xTest))
          .left.map(error => InferenceError.NumericalFailure("RRR held-out prediction", error.message))
      )
      .flatMap { prediction =>
        val baseline = FamilyPredictionMatrices.columnMeans(yTrain)
        val modelError = FamilyPredictionMatrices.squaredError(yTest, prediction)
        val baselineError = FamilyPredictionMatrices.squaredErrorFromRow(yTest, baseline)
        if baselineError <= Math.ulp(1.0) then
          Left(InferenceError.NumericalFailure(
            "RRR predictive gain",
            "held-out baseline has zero squared error"
          ))
        else
          Right(PredictiveGainResult(
            (baselineError - modelError) / baselineError,
            modelError,
            baselineError,
            split,
            validity
          ))
      }

private object FamilyPredictionMatrices:
  def validateFinite(role: String, matrix: DMat): Either[InferenceError, Unit] =
    val values = matrix.copyData
    var i = 0
    while i < values.length do
      if !values(i).isFinite then
        return Left(InferenceError.NonFiniteStatistic(s"$role entry $i", values(i)))
      i += 1
    Right(())

  def columnMeans(matrix: DMat): Array[Double] =
    val means = new Array[Double](matrix.cols)
    var row = 0
    while row < matrix.rows do
      var col = 0
      while col < matrix.cols do
        means(col) += matrix(row, col)
        col += 1
      row += 1
    var col = 0
    while col < matrix.cols do
      means(col) /= matrix.rows.toDouble
      col += 1
    means

  def squaredError(observed: DMat, predicted: DMat): Double =
    require(observed.rows == predicted.rows && observed.cols == predicted.cols)
    var total = 0.0
    var row = 0
    while row < observed.rows do
      var col = 0
      while col < observed.cols do
        val residual = observed(row, col) - predicted(row, col)
        total += residual * residual
        col += 1
      row += 1
    total

  def squaredErrorFromRow(observed: DMat, predicted: Array[Double]): Double =
    require(observed.cols == predicted.length)
    var total = 0.0
    var row = 0
    while row < observed.rows do
      var col = 0
      while col < observed.cols do
        val residual = observed(row, col) - predicted(col)
        total += residual * residual
        col += 1
      row += 1
    total

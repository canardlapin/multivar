package multivar.inference

import gale.backend.Backend.given
import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.LinAlgError
import gale.linalg.Matrix
import gale.linalg.Vec

/** Allocation-aware matrix and vector assembly used by inference protocols. */
private[inference] object InferenceNumerics:
  def matrixFromRowMajor(rows: Int, cols: Int, values: Array[Double]): DMat =
    require(values.length == rows * cols, "row-major data length must match matrix shape")
    val out = Matrix.newBuilder(rows, cols)
    var index = 0
    while index < values.length do
      out.writeLinear(index, values(index))
      index += 1
    out.result()

  def matrixFromRows(rows: Seq[Seq[Double]]): DMat =
    val rowCount = rows.length
    val colCount = rows.headOption.fold(0)(_.length)
    require(rows.forall(_.length == colCount), "matrix rows must have equal lengths")
    val out = Matrix.newBuilder(rowCount, colCount)
    var row = 0
    rows.foreach { values =>
      var col = 0
      values.foreach { value =>
        out(row, col) = value
        col += 1
      }
      row += 1
    }
    out.result()

  def vectorFromArray(values: Array[Double]): DVec =
    val out = Vec.newBuilder(values.length)
    var index = 0
    while index < values.length do
      out(index) = values(index)
      index += 1
    out.result()

  def vectorFromSeq(values: Seq[Double]): DVec =
    val out = Vec.newBuilder(values.length)
    var index = 0
    values.foreach { value =>
      out(index) = value
      index += 1
    }
    out.result()

  def multiply(left: DMat, right: DMat): DMat =
    left * right

  def transposeMultiply(left: DMat, right: DMat): DMat =
    left.t * right

  def selectRows(matrix: DMat, indices: IndexedSeq[Int]): DMat =
    val out = Matrix.newBuilder(indices.length, matrix.cols)
    var row = 0
    while row < indices.length do
      var col = 0
      while col < matrix.cols do
        out(row, col) = matrix(indices(row), col)
        col += 1
      row += 1
    out.result()

extension (matrix: DMat)
  private[inference] def copyData: Array[Double] =
    matrix.valuesRowMajor.toArray

  private[inference] def toRows: Vector[Vector[Double]] =
    Vector.tabulate(matrix.rows) { row =>
      Vector.tabulate(matrix.cols) { col => matrix(row, col) }
    }

  private[inference] def transpose: DMat =
    matrix.t

  private[inference] def selectRows(indices: IndexedSeq[Int]): DMat =
    InferenceNumerics.selectRows(matrix, indices)

extension (vector: DVec)
  private[inference] def toVector: Vector[Double] =
    vector.toSeq.toVector

extension (error: LinAlgError)
  private[inference] def message: String =
    error.getMessage

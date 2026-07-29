package multivar
package family.paired

import multivar.core.*

import gale.linalg.DMat
import gale.linalg.DVec

enum TaskPart:
  case Focus
  case Remainder

/** Exact in-sample quadratic loss accounting for one fitted task.
  *
  * Gains are reductions in the declared training objective, not held-out utility
  * or scientific worth.
  */
final class TaskValue private[paired] (
    val baselineLoss: Double,
    val fullLoss: Double,
    val fittedLoss: Double,
    val availableGain: Double,
    val retainedGain: Double
):
  def unretainedGain: Double =
    math.max(0.0, availableGain - retainedGain)

  def availableFraction: Double =
    safeRatio(availableGain, baselineLoss)

  def retainedFraction: Double =
    safeRatio(retainedGain, baselineLoss)

  def rankEfficiency: Double =
    safeRatio(retainedGain, availableGain)

  private def safeRatio(numerator: Double, denominator: Double): Double =
    if denominator > 0.0 then numerator / denominator else 0.0

/** Per-part gain summary. */
final class TaskPartValue private[paired] (
    val availableGain: Double,
    val retainedGain: Double,
    val shareOfAvailable: Double
)

/** One orthogonal part of a resolved task. */
final class TaskComponentBlock private[paired] (
    val part: TaskPart,
    private val singular: DVec,
    private val source: DMat,
    private val target: DMat,
    private val trainingSourceScores: DMat,
    val value: TaskPartValue
):
  def singularValues: DVec =
    singular

  /** Predictor-side directions in working coordinates. */
  def sourceWeights: DMat =
    source

  /** Response-side patterns in working coordinates. */
  def targetPatterns: DMat =
    target

  /** Training scores with singular-value scale included. */
  def sourceScores: DMat =
    trainingSourceScores

  def rank: Int =
    singular.length

  def componentValues: DVec =
    val out = new Array[Double](singular.length)
    var i = 0
    while i < singular.length do
      out(i) = singular(i) * singular(i)
      i += 1
    GaleNumerics.vectorFromArray(out)

object TaskValue:
  private[paired] def apply(
      baselineLoss: Double,
      fullLoss: Double,
      fittedLoss: Double,
      availableGain: Double,
      retainedGain: Double
  ): TaskValue =
    new TaskValue(baselineLoss, fullLoss, fittedLoss, availableGain, retainedGain)

object TaskPartValue:
  private[paired] def apply(
      availableGain: Double,
      retainedGain: Double,
      totalAvailable: Double
  ): TaskPartValue =
    val share =
      if totalAvailable > 0.0 then availableGain / totalAvailable else 0.0
    new TaskPartValue(availableGain, retainedGain, share)

object TaskComponentBlock:
  private[paired] def apply(
      part: TaskPart,
      singular: DVec,
      source: DMat,
      target: DMat,
      trainingSourceScores: DMat,
      value: TaskPartValue
  ): TaskComponentBlock =
    new TaskComponentBlock(part, singular, source, target, trainingSourceScores, value)

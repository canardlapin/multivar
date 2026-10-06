package multivar.inference

import gale.linalg.Matrix
import gale.spectral.{MatrixEnclosure, RealInterval}
import munit.FunSuite

class TaskRankPlanJvmSuite extends FunSuite:
  private def accepted[E, A](value: Either[E, A]): A =
    value.fold(error => fail(error.toString), identity)

  private val group = TaskRankGroup(
    "group",
    Vector(
      TaskRankParticipant("p1", Matrix.dense(1, 1)(1.0)),
      TaskRankParticipant("p2", Matrix.dense(1, 1)(2.0))
    ),
    accepted(MatrixEnclosure.exact(Matrix.eye(1))),
    accepted(RealInterval.exact(1.0))
  )
  private val coordinates = TaskRankCoordinates("fixture/v1", "identity", "JVM constructor boundary", Vector("f"))

  private def forged(draws: Vector[TaskRankSigns]): TaskRankPlan =
    val constructor = classOf[TaskRankPlan].getDeclaredConstructors
      .find(_.getParameterCount == 3)
      .getOrElse(fail("TaskRankPlan JVM constructor was not found"))
    constructor.newInstance(Vector(2), draws, TaskRankSampling.ExhaustiveFullSignGroup)
      .asInstanceOf[TaskRankPlan]

  private def result(plan: TaskRankPlan): Either[InferenceError, TaskRankReceipt] =
    TaskMeanRank.run(Vector(group), coordinates, plan, TaskRankMethod.FrobeniusFeasibility, 1, 20)

  test("JVM constructor bypass cannot forge incomplete or duplicate exhaustive orbits") {
    val empty = forged(Vector.empty)
    val incomplete = forged(Vector(TaskRankSigns("0", Vector(Vector(-1, -1)))))
    val duplicate = forged(Vector(
      TaskRankSigns("0", Vector(Vector(-1, -1))),
      TaskRankSigns("1", Vector(Vector(-1, -1))),
      TaskRankSigns("2", Vector(Vector(1, -1))),
      TaskRankSigns("3", Vector(Vector(1, 1)))
    ))
    Vector(empty, incomplete, duplicate).foreach { plan =>
      result(plan) match
        case Left(InferenceError.InvalidReplicatePlan(_)) => ()
        case other => fail(s"expected forged exhaustive plan refusal, got $other")
    }
  }

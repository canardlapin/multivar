package multivar.inference

import gale.linalg.{DMat, Matrix}
import gale.spectral.{MatrixEnclosure, RealInterval}
import munit.FunSuite
import resample4s.kernel.Seed

class TaskMeanRankSuite extends FunSuite:
  private def accepted[E, A](value: Either[E, A]): A = value.fold(e => fail(e.toString), identity)
  private def group(id: String, profiles: Vector[DMat], weight: Double = 1.0): TaskRankGroup =
    TaskRankGroup(id, profiles.zipWithIndex.map((x, i) => TaskRankParticipant(s"$id-$i", x)),
      accepted(MatrixEnclosure.exact(Matrix.eye(profiles.head.rows))), accepted(RealInterval.exact(weight)))
  private def coordinates(features: Int): TaskRankCoordinates =
    TaskRankCoordinates("synthetic/v1", "fixed identity row map", "analytic exact-binary fixtures", Vector.tabulate(features)(i => s"feature-$i"))
  private def run(groups: Vector[TaskRankGroup], plan: TaskRankPlan, method: TaskRankMethod,
      a: BigInt = BigInt(1), b: BigInt = BigInt(20)): TaskRankReceipt =
    accepted(TaskMeanRank.run(groups, coordinates(groups.head.participants.head.profile.cols), plan, method, a, b))

  test("noiseless two-group floor distinguishes the proved procedures") {
    val groups = Vector(group("a", Vector.fill(6)(Matrix.eye(2))), group("b", Vector.fill(6)(Matrix.eye(2))))
    val plan = accepted(TaskRankPlan.exhaustive(Vector(6, 6), 4096))
    val frobenius = run(groups, plan, TaskRankMethod.FrobeniusFeasibility)
    val operator = run(groups, plan, TaskRankMethod.OperatorFirstThenBound)
    assertEquals(frobenius.counts(1), BigInt(4))
    assertEquals(operator.counts(1), BigInt(252))
    assertEquals(frobenius.denominator, BigInt(4096))
    assertEquals(frobenius.lowerRankBound, 2)
    assertEquals(operator.lowerRankBound, 1)
    Vector(frobenius, operator).foreach { receipt =>
      assertEquals(receipt.draws.map(_.signs), plan.draws.map(_.groups))
      assert(receipt.counts.sliding(2).forall(pair => pair(0) <= pair(1)))
      assertEquals(receipt.counts.last, receipt.denominator)
    }
  }

  test("fixed B adds identity and decides exact rational boundaries without display rounding") {
    val profiles = Vector.fill(2)(Matrix.dense(1, 1)(1.0))
    val groups = Vector(group("a", profiles))
    val plan = accepted(TaskRankPlan.replayFixed(Vector(2), Vector(TaskRankSigns("0", Vector(Vector(-1, 1))))))
    TaskRankMethod.values.foreach { method =>
      val tie = run(groups, plan, method, 1, 2)
      assertEquals(tie.counts, Vector(BigInt(1), BigInt(2)))
      assertEquals(tie.denominator, BigInt(2))
      assertEquals(tie.lowerRankBound, 1)
      // This alpha rounds to .5 for display, but is strictly below 1/2.
      val denominator = BigInt(1) << 60
      val below = run(groups, plan, method, denominator / 2 - 1, denominator)
      assertEquals(below.lowerRankBound, 0)
    }
  }

  test("weighted rectangular rank-two fixture retains informative later-root exclusions") {
    val tall = Matrix.dense(3, 5)(3, 4, 0, 0, 0, 4, -3, 0, 0, 0, 0, 0, 0, 0, 0)
    val wide = Matrix.dense(2, 5)(2, 0, 0, 0, 0, 0, 1, 0, 0, 0)
    val groups = Vector(group("a", Vector.fill(3)(tall), .5), group("b", Vector.fill(4)(wide), 2.0))
    val plan = accepted(TaskRankPlan.replayFixed(Vector(3, 4), Vector(TaskRankSigns("mixed", Vector(Vector(1, 1, -1), Vector(1, -1, 1, -1))))))
    val receipt = run(groups, plan, TaskRankMethod.FrobeniusFeasibility)
    assertEquals(receipt.counts(1), BigInt(1))
    assertEquals(receipt.counts(2), BigInt(2))
    assertEquals(receipt.groups.map(_.participantIds.size), Vector(3, 4))
  }

  test("singleton and full-sign ties always count") {
    val groups = Vector(group("single", Vector(Matrix.dense(1, 2)(3, -2))))
    val plan = accepted(TaskRankPlan.exhaustive(Vector(1), 2))
    TaskRankMethod.values.foreach { method =>
      val result = run(groups, plan, method)
      assertEquals(result.counts, Vector(BigInt(2), BigInt(2)))
      assertEquals(result.lowerRankBound, 0)
    }
  }

  test("native fixed sampling replays exact signs and remains a fixed budget") {
    val first = accepted(TaskRankPlan.fixedUniform(Vector(3, 2), 31, Seed.fromLong(713L)))
    val second = accepted(TaskRankPlan.fixedUniform(Vector(3, 2), 31, Seed.fromLong(713L)))
    assertEquals(first.draws, second.draws)
    assertEquals(first.draws.size, 31)
    assert(first.draws.forall(_.groups.map(_.size) == Vector(3, 2)))
    assert(first.draws.flatMap(_.groups.flatten).forall(s => s == -1 || s == 1))
    val groups = Vector(group("a", Vector.fill(3)(Matrix.eye(2))), group("b", Vector.fill(2)(Matrix.eye(2))))
    val replay = accepted(TaskRankPlan.replayFixed(Vector(3, 2), first.draws))
    assertEquals(run(groups, first, TaskRankMethod.FrobeniusFeasibility).counts,
      run(groups, replay, TaskRankMethod.FrobeniusFeasibility).counts)
  }

  test("malformed identities, incomplete signs, unsafe scales and cancellation refuse") {
    val groups = Vector(group("a", Vector.fill(2)(Matrix.eye(2))))
    val plan = accepted(TaskRankPlan.exhaustive(Vector(2), 4))
    assert(TaskRankPlan.exhaustive(Vector(2), 3).isLeft)
    assert(TaskRankPlan.replayFixed(Vector(2), Vector(TaskRankSigns("x", Vector(Vector(1))))).isLeft)
    assert(TaskRankPlan.replayFixed(Vector(2), Vector.fill(2)(TaskRankSigns("x", Vector(Vector(1, -1))))).isLeft)
    val duplicate = groups.head.copy(participants = Vector.fill(2)(groups.head.participants.head))
    assert(TaskMeanRank.run(Vector(duplicate), coordinates(2), plan, TaskRankMethod.FrobeniusFeasibility, 1, 20).isLeft)
    var checked = 0
    val cancelled = TaskMeanRank.run(groups, coordinates(2), plan, TaskRankMethod.FrobeniusFeasibility, 1, 20,
      () => { checked += 1; checked >= 3 })
    assert(cancelled.isLeft)
    val huge = Vector(group("huge", Vector.fill(2)(Matrix.dense(1, 1)(Double.MaxValue))))
    assert(TaskMeanRank.run(huge, coordinates(1), plan, TaskRankMethod.FrobeniusFeasibility, 1, 20).isLeft)
  }

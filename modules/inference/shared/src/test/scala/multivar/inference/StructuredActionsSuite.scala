package multivar.inference

import multivar.core.{RowProjector, RowWhitening}

import gale.linalg.DMat
import resample4s.kernel.{Permutation, Seed}

class StructuredActionsSuite extends munit.FunSuite:

  private def accepted[A](value: Either[InferenceError, A]): A =
    value.fold(error => fail(error.message), identity)

  private def acceptedMultivar[A](
      value: Either[multivar.core.MultivarError, A]
  ): A =
    value.fold(error => fail(error.message), identity)

  private def values(permutation: Permutation): Vector[Int] =
    permutation.toVector

  private def permutation(values: Int*): Permutation =
    Permutation
      .from(IArray.from(values))
      .fold(error => fail(error.message), identity)

  test("nested partitions prove that each inner unit stays inside one stratum") {
    val rows = accepted(RowCount(6))
    val clusters = accepted(ClusterPartition.from(
      rows,
      Vector(Vector(0, 1), Vector(2), Vector(3, 4), Vector(5))
    ))
    val strata = accepted(StrataPartition.from(
      rows,
      Vector(Vector(0, 1, 2), Vector(3, 4, 5))
    ))

    assert(NestedPartitions.from(clusters.value, strata.value).isRight)

    val crossing = accepted(ClusterPartition.from(
      rows,
      Vector(Vector(0, 3), Vector(1, 2), Vector(4, 5))
    ))
    NestedPartitions.from(crossing.value, strata.value) match
      case Left(_: InferenceError.Resampling) => ()
      case other => fail(s"expected a nesting failure, got $other")
  }

  test("within-block and within-stratum actions never cross declared groups") {
    val rows = accepted(RowCount(8))
    val blocks = accepted(RowPartition.from(
      rows,
      Vector(Vector(0, 1), Vector(2, 3, 4), Vector(5, 6, 7))
    ))
    val strata = accepted(StrataPartition.from(
      rows,
      Vector(Vector(0, 2, 4, 6), Vector(1, 3, 5, 7))
    ))
    val seed = Seed.fromLong(9123L)
    val replicate = accepted(ReplicateId(7))
    val blockDraw = accepted(PermutationAction.withinBlocks(blocks).draw(seed, replicate))
    val strataDraw = accepted(PermutationAction.withinStrata(strata).draw(seed, replicate))

    def groupMap(partition: RowPartition): Array[Int] = partition.groupIndexByRow
    def staysWithin(permutation: Permutation, partition: RowPartition): Boolean =
      val groups = groupMap(partition)
      permutation.toVector.indices.forall { target =>
        groups(target) == groups(permutation.toVector(target))
      }

    assert(staysWithin(blockDraw, blocks))
    assert(staysWithin(strataDraw, strata.value))
    assertEquals(values(blockDraw).sorted, Vector.range(0, 8))
    assertEquals(values(strataDraw).sorted, Vector.range(0, 8))
  }

  test("whole-cluster permutation preserves within-cluster order and refuses unequal sizes") {
    val rows = accepted(RowCount(6))
    val equal = accepted(ClusterPartition.from(
      rows,
      Vector(Vector(0, 1), Vector(2, 3), Vector(4, 5))
    ))
    val action = accepted(PermutationAction.wholeClusters(equal))
    val draw = accepted(action.draw(Seed.fromLong(17L), accepted(ReplicateId(2))))

    val pairs = values(draw).grouped(2).toVector
    assert(pairs.forall(pair => pair(1) == pair(0) + 1))
    assertEquals(values(draw).sorted, Vector.range(0, 6))

    val unequal = accepted(ClusterPartition.from(
      rows,
      Vector(Vector(0), Vector(1, 2), Vector(3, 4, 5))
    ))
    PermutationAction.wholeClusters(unequal) match
      case Left(_: InferenceError.UnsupportedProblem) => ()
      case other => fail(s"expected unsupported unequal clusters, got $other")
  }

  test("row permutations apply a closed-form source-row fixture") {
    val input = InferenceNumerics.matrixFromRows(Vector(
      Vector(10.0, 100.0),
      Vector(20.0, 200.0),
      Vector(30.0, 300.0)
    ))
    val rowOrder = permutation(2, 0, 1)
    val out = accepted(rowOrder.applyTo(input))

    assertEquals(out.toRows, Vector(
      Vector(30.0, 300.0),
      Vector(10.0, 100.0),
      Vector(20.0, 200.0)
    ))
  }

  test("nuisance residual randomization preserves fitted block means") {
    val rows = accepted(RowCount(4))
    val blocks = accepted(RowPartition.from(rows, Vector(Vector(0, 1), Vector(2, 3))))
    val reference = accepted(ConditioningRef("block-means-v1"))
    val design = ResamplingDesign.nuisanceAdjustedWithinBlocks(
      blocks,
      reference,
      WhiteningRequirement.NotRequired
    )
    val nuisance = acceptedMultivar(RowProjector.fromMatrix(InferenceNumerics.matrixFromRows(Vector(
      Vector(0.5, 0.5, 0.0, 0.0),
      Vector(0.5, 0.5, 0.0, 0.0),
      Vector(0.0, 0.0, 0.5, 0.5),
      Vector(0.0, 0.0, 0.5, 0.5)
    ))))
    val action = accepted(ResidualPermutationAction.from(
      design,
      reference,
      nuisance,
      None
    ))
    val input = InferenceNumerics.matrixFromRows(Vector(
      Vector(10.0),
      Vector(12.0),
      Vector(20.0),
      Vector(24.0)
    ))
    val out = accepted(action.draw(input, Seed.fromLong(44L), accepted(ReplicateId(3))))

    assertEqualsDouble(out(0, 0) + out(1, 0), 22.0, 1e-12)
    assertEqualsDouble(out(2, 0) + out(3, 0), 44.0, 1e-12)
    assertEquals(
      Vector(out(0, 0), out(1, 0)).sorted,
      Vector(10.0, 12.0)
    )
    assertEquals(
      Vector(out(2, 0), out(3, 0)).sorted,
      Vector(20.0, 24.0)
    )
  }

  test("replicate-local actions are deterministic and independent of draw order") {
    val rows = accepted(RowCount(12))
    val action = PermutationAction.unrestricted(rows)
    val root = Seed.fromLong(1729L)
    val zero = accepted(ReplicateId(0))
    val seven = accepted(ReplicateId(7))

    val sevenFirst = accepted(action.draw(root, seven))
    val zeroAfter = accepted(action.draw(root, zero))
    val zeroFirst = accepted(action.draw(root, zero))
    val sevenAfter = accepted(action.draw(root, seven))

    assertEquals(zeroAfter, zeroFirst)
    assertEquals(sevenFirst, sevenAfter)
    assertNotEquals(zeroFirst, sevenFirst)
    assertEquals(values(zeroFirst).sorted, Vector.range(0, 12))
    assertEquals(values(sevenFirst).sorted, Vector.range(0, 12))
  }

  test("partitions canonicalize caller label identities") {
    val rows = accepted(RowCount(6))
    val first = accepted(RowPartition.from(
      rows,
      Vector(Vector(0, 3), Vector(1, 2, 4), Vector(5))
    ))
    val relabeled = accepted(RowPartition.from(
      rows,
      Vector(Vector(5), Vector(1, 2, 4), Vector(0, 3))
    ))

    assertEquals(first.labels, relabeled.labels)
    assertEquals(first.groupIndexByRow.toVector, relabeled.groupIndexByRow.toVector)
  }

  test("required whitening and conditioning identity fail closed") {
    val rows = accepted(RowCount(2))
    val reference = accepted(ConditioningRef("nuisance-v1"))
    val wrongReference = accepted(ConditioningRef("nuisance-v2"))
    val design = ResamplingDesign.nuisanceAdjusted(
      rows,
      reference,
      WhiteningRequirement.Required
    )
    val nuisance = acceptedMultivar(RowProjector.zero(2))

    ResidualPermutationAction.from(design, reference, nuisance, None) match
      case Left(_: InferenceError.UnsupportedProblem) => ()
      case other => fail(s"expected a missing-whitening rejection, got $other")

    val whitening = acceptedMultivar(RowWhitening.identity(2))
    ResidualPermutationAction.from(design, wrongReference, nuisance, Some(whitening)) match
      case Left(_: InferenceError.UnsupportedProblem) => ()
      case other => fail(s"expected a conditioning-reference rejection, got $other")
  }

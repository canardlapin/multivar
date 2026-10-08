package multivar.inference

import gale.backend.Backend.given
import gale.linalg.DMat
import resample4s.kernel.{Permutation, Seed}

class CanonicalRankGeometrySuite extends munit.FunSuite:
  private def right[A](value: Either[?, A]): A = value.fold(e => fail(e.toString), identity)
  private val input = CanonicalRankV2Fixtures.input
  private val x = DMat.tabulate(6,2)((i,j) => input(i,j))
  private val y = DMat.tabulate(6,3)((i,j) => input(i,j+2))
  private def close(a: Double, b: Double): Unit = assertEqualsDouble(a,b,1e-10 * (1.0 + math.abs(b)))

  test("all 720 independent R score-completed actions survive feature coordinate changes"):
    val transforms = Vector(
      DMat.eye(3),
      DMat.dense(3,3,Vector(1.0,0.0,2.0,0.0,1.0,0.0,0.0,0.0,1.0)),
      DMat.dense(3,3,Vector(2.0,0.0,0.0,0.0,.5,0.0,0.0,0.0,3.0)),
      DMat.dense(3,3,Vector(0.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,0.0)),
      DMat.dense(3,3,Vector(.8,-.6,0.0,.6,.8,0.0,0.0,0.0,1.0)))
    val leftChange = DMat.dense(2,2,Vector(1.0,.7,0.0,2.0))
    transforms.foreach: transform =>
      val problem = right(StepwiseCanonicalRank.from(x * leftChange, y * transform))
      assertEquals(problem.method, CanonicalRankMethod.ScoreOrthogonalPermutationV2)
      problem.correlations.zip(Vector(.8,.2)).foreach((a,b) => close(a,b))
      problem.observedWilks.zip(CanonicalRankV2Fixtures.observed).foreach((a,b) => close(a,b))
      for variables <- Vector(problem.leftVariables,problem.rightVariables) do
        val gram = variables.t * variables
        for i <- 0 until gram.rows; j <- 0 until gram.cols do close(gram(i,j),if i==j then 1.0 else 0.0)
      val counts = Array(0,0)
      CanonicalRankV2Fixtures.actions.foreach: (rows,expected) =>
        val values = right(problem.nullStatistics(right(Permutation.from(IArray.from(rows)))))
        values.zip(expected).foreach((a,b) => close(a,b))
        for k <- 0 to 1 do
          if values(k) >= problem.observedWilks(k)-1e-12 then counts(k) += 1
      assertEquals(counts.toVector,Vector(460,642))

  test("block interchange with inverse actions and equal-dimensional changes preserve tails"):
    val original = right(StepwiseCanonicalRank.from(x,y))
    val swapped = right(StepwiseCanonicalRank.from(y,x))
    val yy = DMat.tabulate(6,2)((i,j) => y(i,j))
    val equal = right(StepwiseCanonicalRank.from(x,yy))
    val changed = right(StepwiseCanonicalRank.from(x * DMat.dense(2,2,Vector(1.0,2.0,0.0,1.0)),
      yy * DMat.dense(2,2,Vector(.8,-.6,.6,.8))))
    CanonicalRankV2Fixtures.actions.foreach: (rows,_) =>
      val forward = right(Permutation.from(IArray.from(rows)))
      val inverse = right(Permutation.from(IArray.from(Vector.tabulate(6)(i => rows.indexOf(i)))))
      right(original.nullStatistics(forward)).zip(right(swapped.nullStatistics(inverse))).foreach((a,b) => close(a,b))
      right(equal.nullStatistics(forward)).zip(right(changed.nullStatistics(forward))).foreach((a,b) => close(a,b))

  test("unidentified tied-root cuts refuse without losing observed zero-root hypotheses"):
    def walsh(i: Int, j: Int) = if Integer.bitCount(i & j) % 2 == 0 then 1.0 else -1.0
    val left = DMat.tabulate(8,2)((i,j) => walsh(i,j+1))
    for rho <- Vector(0.0,.6,.60000000001) do
      val rightSide = DMat.tabulate(8,2)((i,j) => rho*walsh(i,j+1) + math.sqrt(1-rho*rho)*walsh(i,j+3))
      val problem = right(StepwiseCanonicalRank.from(left,rightSide))
      assertEquals(problem.correlations.size,2)
      assertEquals(problem.observedWilks.size,2)
      val refused = FixedCanonicalRank.run(problem,PermutationAction.unrestricted(right(RowCount(8))),
        Seed.fromLong(42L),right(MonteCarloDraws(9)),right(Alpha(.05)))
      assert(refused.left.toOption.exists(_.isInstanceOf[InferenceError.UnidentifiedCanonicalTail]))
      assert(CanonicalRankSpectrum.from(left,rightSide).isRight)

  test("nearly singular and nonfinite candidates refuse before numerical interpretation"):
    val singular = DMat.tabulate(6,2)((i,j) => if j==0 then x(i,0) else x(i,0)+1e-15*x(i,1))
    assert(StepwiseCanonicalRank.from(singular,y).isLeft)
    assert(CanonicalRankSpectrum.from(singular,y).isLeft)
    assert(CanonicalRankSpectrum.from(DMat.tabulate(6,2)((_,_) => Double.NaN),y).isLeft)

package multivar.inference

import gale.backend.Backend.given
import gale.linalg.DMat
import resample4s.kernel.{Seed, StreamDomain}

class GaussianCanonicalRankSuite extends munit.FunSuite:
  private def right[A](value: Either[?, A]): A = value.fold(e => fail(e.toString), identity)
  private val input = CanonicalRankV2Fixtures.input
  private val x = DMat.tabulate(6,2)((i,j) => input(i,j))
  private val y = DMat.tabulate(6,3)((i,j) => input(i,j+2))
  private val seed = Seed.fromLong(871254L)
  private def close(a: Double, b: Double): Unit = assertEqualsDouble(a,b,1e-10 * (1.0 + math.abs(b)))

  test("Gaussian reference seed paths, normal inputs and every Wilks draw match independent Python and R"):
    val problem = right(CanonicalRankSpectrum.from(x,y))
    val result = right(GaussianCanonicalRank.run(problem,seed,right(MonteCarloDraws(39)),right(Alpha(.05))))
    assertEquals(result.method,CanonicalRankMethod.GaussianInterlacingWilksV1)
    assertEquals(result.hypothesisSeeds.map(_.value),CanonicalRankV2Fixtures.hypothesisSeeds)
    assertEquals(result.references.map(r => (r.rows,r.smallerColumns,r.largerColumns,r.nullRank)),Vector((6,2,3,0),(6,1,3,1)))
    val local = GaussianCanonicalRank.child(result.hypothesisSeeds.head,StreamDomain.custom(1012).toOption.get,0)
    val block = GaussianCanonicalRank.gaussian(6,2,GaussianCanonicalRank.child(local,StreamDomain.custom(1013).toOption.get,0))
    for i <- 0 until 6; j <- 0 until 2 do close(block(i,j),CanonicalRankV2Fixtures.firstNormalBlock(i*2+j))
    result.receipts.zipWithIndex.foreach: (receipt,k) =>
      receipt.nullValues.zip(CanonicalRankV2Fixtures.gaussian(k)).foreach((a,b) => close(a,b))
      close(receipt.observed,CanonicalRankV2Fixtures.observed(k))
      val count = CanonicalRankV2Fixtures.gaussian(k).count(_ >= CanonicalRankV2Fixtures.observed(k))
      assertEquals(receipt.exceedances,count)
      assertEquals(receipt.consumed.value,39)
      assertEqualsDouble(receipt.pValue.value,(count+1.0)/40,0.0)
      assertEquals(receipt.provenance,ReplicateProvenance.Deterministic(result.hypothesisSeeds(k),GaussianCanonicalRank.referenceAlgorithm,Seed.derivationAlgorithm))
    val closed = result.receipts.map(_.pValue.value).scanLeft(0.0)(math.max).tail
    assertEquals(result.adjustedPValues.map(_.value),closed)
    assertEquals(result.detectableRank,closed.takeWhile(_ <= .05).size)
    assertEquals(result.completedCompactFits,78L)

  test("reference streams are stable under B changes, scheduling, feature transforms and block interchange"):
    val a = right(CanonicalRankSpectrum.from(x,y))
    val b = right(CanonicalRankSpectrum.from(y * DMat.dense(3,3,Vector(1.0,0.0,2.0,0.0,1.0,0.0,0.0,0.0,1.0)),x))
    val first = right(GaussianCanonicalRank.run(a,seed,right(MonteCarloDraws(19)),right(Alpha(.05))))
    val second = right(GaussianCanonicalRank.run(b,seed,right(MonteCarloDraws(39)),right(Alpha(.05))))
    for k <- 0 to 1 do
      first.receipts(k).nullValues.zip(second.receipts(k).nullValues).foreach((u,v) => close(u,v))
      close(first.receipts(k).observed,second.receipts(k).observed)
      Vector(18,0,7).foreach: index =>
        val draw = right(first.references(k).draw(first.hypothesisSeeds(k),right(ReplicateId(index))))
        close(draw.observedWilks.head,first.receipts(k).nullValues(index))
    assertNotEquals(first.hypothesisSeeds(0),first.hypothesisSeeds(1))

  test("one-column Gaussian reference follows the analytic uncentered Beta(1,2) law"):
    // m=6, a=1, b=2: r^2~Beta(b/2,(m-b)/2), E[W]=1/2,
    // P(W >= -log(1-.8^2))=(1-.8^2)^2=.1296. Extra centering changes this law.
    val left = DMat.tabulate(6,1)((i,_) => x(i,0))
    val rightSide = DMat.tabulate(6,2)((i,j) => y(i,j))
    val problem = right(CanonicalRankSpectrum.from(left,rightSide))
    val result = right(GaussianCanonicalRank.run(problem,Seed.fromLong(635827L),right(MonteCarloDraws(2047)),right(Alpha(.05))))
    val receipt = result.receipts.head
    assertEqualsDouble(receipt.nullValues.sum/2047,.5,.04)
    assertEqualsDouble(receipt.pValue.value,.1296,.03)
    assertEqualsDouble(receipt.observed,-math.log(.36),1e-12)

  test("tied and zero roots remain available without fitting unidentified tail directions"):
    val left = DMat.tabulate(8,2)((i,j) => if i==j then 1.0 else 0.0)
    for rho <- Vector(0.0,.6) do
      val rightSide = DMat.tabulate(8,2)((i,j) => if i==j then rho else if i==j+2 then math.sqrt(1-rho*rho) else 0.0)
      val result = right(GaussianCanonicalRank.run(right(CanonicalRankSpectrum.from(left,rightSide)),seed,right(MonteCarloDraws(19)),right(Alpha(.05))))
      assertEquals(result.receipts.size,2)
      if rho==0.0 then result.receipts.foreach(r => assertEqualsDouble(r.pValue.value,1.0,0.0))

  test("interlacing bounds every observed tail by a reduced-block full Wilks statistic"):
    val whole = right(CanonicalRankSpectrum.from(x,y))
    val reduced = right(CanonicalRankSpectrum.from(DMat.tabulate(6,1)((i,_) => x(i,1)),y))
    assert(whole.correlations(1) <= reduced.correlations.head + 1e-12)
    assert(whole.observedWilks(1) <= reduced.observedWilks.head + 1e-12)

  test("reference dimensions, memory and receipt budgets refuse invalid requests"):
    val problem = right(CanonicalRankSpectrum.from(x,y))
    val count = right(MonteCarloDraws(39));val alpha = right(Alpha(.05))
    assert(GaussianCanonicalReference.forHypothesis(problem,right(ComponentIx(2))).isLeft)
    assert(GaussianCanonicalRank.run(problem,seed,count,alpha,maximumNullValues=77).isLeft)
    assert(GaussianCanonicalRank.run(problem,seed,count,alpha,maximumElements=0).isLeft)
    val reference = right(GaussianCanonicalReference.forHypothesis(problem,right(ComponentIx(1))))
    assert(reference.draw(seed,right(ReplicateId(0)),maximumElements=0).isLeft)

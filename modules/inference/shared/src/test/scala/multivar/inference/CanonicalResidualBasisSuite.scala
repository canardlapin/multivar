package multivar.inference

import gale.backend.Backend.given
import gale.linalg.{DMat, QROptions, QRPivoting}

class CanonicalResidualBasisSuite extends munit.FunSuite:
  private def ok[A](value: Either[?, A]): A = value.fold(e => fail(e.toString), identity)
  private def same(a: DMat, b: DMat): Unit =
    assertEquals((a.rows, a.cols), (b.rows, b.cols))
    for i <- 0 until a.rows; j <- 0 until a.cols do
      assertEqualsDouble(a(i, j), b(i, j), 1e-12)

  for rows <- Vector(12, 80) do
    test(s"Huh-Jhun ignores redundant nuisance reflectors at $rows rows"):
      val unique = DMat.tabulate(rows, 3): (i, j) =>
        if j == 0 then 1.0 else if j == 1 then (i % 7 - 3) * 0.37 else (i * i % 11 - 5) * 0.23
      val redundant = DMat.tabulate(rows, 4)((i, j) => if j == 0 then 1.0 else unique(i, j - 1))
      val qr = redundant.qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(1e-12)))
      val rank = qr.diagnostics.rank.get
      assertEquals(rank, 3)
      val selected = DMat.tabulate(rows, rank)((i, j) => redundant(i, qr.columnPermutation(j)))
      val selectors = DMat.tabulate(rows, rows - rank)((i, j) => if i == j + rank then 1.0 else 0.0)
      val expected = ok(selected.qr(QROptions(pivoting = QRPivoting.Disabled)).applyQ(selectors))
      val actual = ok(CanonicalResidualBasis.from(redundant))
      assertEquals(actual.nuisanceRank, rank)
      same(actual.matrix, expected)
      val projector = actual.matrix * actual.matrix.t
      val base = ok(CanonicalResidualBasis.from(unique)).matrix
      same(projector, base * base.t)

  test("full-rank Huh-Jhun retains the existing pivoted QR coordinates"):
    val z = DMat.tabulate(12, 2)((i, j) => if j == 0 then 1.0 else (i - 3) * 0.2)
    val qr = z.qr(QROptions(pivoting = QRPivoting.Column, rankTolerance = Some(1e-12)))
    val selectors = DMat.tabulate(12, 10)((i, j) => if i == j + 2 then 1.0 else 0.0)
    same(ok(CanonicalResidualBasis.from(z)).matrix, ok(qr.applyQ(selectors)))

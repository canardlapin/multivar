package multivar.inference

import gale.linalg.DMat

/** A single-copy deletion from an identified N-of-N source multiset. The source
  * domain still has N units; the empirical mean has N-1 observations.
  */
final class ProfileDeletionCounts private (
    val source: ProfileBootstrapCounts,
    val removedUnit: Int,
    val multiplicities: Vector[Int]
):
  def units: Int = source.units
  def sampleSize: Int = units - 1

  def meanProfiles(profiles: DMat, coordinatesPerUnit: Int,
      maxWorkingBytes: Long = 256L * 1024 * 1024,
      cancelled: () => Boolean = () => false): Either[InferenceError, DMat] =
    if coordinatesPerUnit < 1 || profiles.cols < 1 ||
        profiles.rows.toLong != units.toLong * coordinatesPerUnit then
      Left(InferenceError.InvalidReplicatePlan("deletion readout requires the original complete source-unit domain"))
    else for
      _ <- ProfileBootstrapMath.allocation(BigInt(32) * coordinatesPerUnit * profiles.cols, maxWorkingBytes)
      _ <- ProfileBootstrapMath.finite(profiles, "deletion mean readout", cancelled)
      mean <- ProfileBootstrapMath.aggregate(profiles, multiplicities, sampleSize, coordinatesPerUnit, cancelled)
    yield mean

object ProfileDeletionCounts:
  def fromBootstrap(source: ProfileBootstrapCounts, removedUnit: Int): Either[InferenceError, ProfileDeletionCounts] =
    if removedUnit < 0 || removedUnit >= source.units then
      Left(InferenceError.InvalidReplicatePlan("deleted unit is outside the original sampling-unit domain"))
    else if source.multiplicities(removedUnit) == 0 then
      Left(InferenceError.InvalidReplicatePlan("cannot delete a copy absent from the bootstrap multiset"))
    else Right(new ProfileDeletionCounts(source, removedUnit,
      source.multiplicities.updated(removedUnit, source.multiplicities(removedUnit) - 1)))

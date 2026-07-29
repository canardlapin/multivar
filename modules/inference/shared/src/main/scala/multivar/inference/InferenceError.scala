package multivar.inference

import resample4s.spi.DesignError

enum ExecutionPolicyMismatch:
  case FixedMonteCarloUnsupported
  case PerRungDraws
  case Alpha
  case BatchSize
  case SequentialBoundary

enum InferenceError:
  case InvalidIdentifier(role: String, value: String)
  case InvalidDescription(role: String, value: String)
  case InvalidCount(role: String, value: Int)
  case InvalidNonNegativeCount(role: String, value: Int)
  case InvalidProbability(role: String, value: Double, inclusiveZero: Boolean)
  case InvalidTolerance(role: String, value: Double)
  case InvalidComponent(value: Int)
  case EmptyComponentSet(role: String)
  case DuplicateComponent(value: Int)
  case UnorderedComponents(previous: Int, next: Int)
  case ComponentOutOfRange(value: Int, rank: Int)
  case InvalidUnit(detail: String)
  case InvalidPartition(detail: String)
  case RowCountMismatch(role: String, expected: Int, actual: Int)
  case InvalidSpectrum(detail: String)
  case InvalidReplicatePlan(detail: String)
  case Resampling(error: DesignError)
  case NumericalFailure(role: String, detail: String)
  case InvalidValidity(detail: String)
  case UnsupportedProblem(detail: String)
  case FitProtocolMismatch(
      compiled: FitDescriptor[?],
      protocol: FitDescriptor[?]
  )
  case TargetProtocolMismatch(
      compiled: TargetLabel,
      protocol: TargetLabel
  )
  case NullProtocolMismatch(
      compiled: NullLabel,
      protocol: NullLabel
  )
  case ExecutionPolicyConflict(
      mismatch: ExecutionPolicyMismatch
  )
  case FailedAssumption(id: AssumptionId, detail: String)
  case NonFiniteStatistic(role: String, value: Double)
  case RankLoss(expected: Int, actual: Int)
  case UnitBeyondRank(unit: UnitId, rank: Int)
  case ReplicateFailure(replicate: ReplicateId, cause: InferenceError)
  case BudgetExhausted(consumed: Int, allocated: Int)
  case UnsupportedEvidence(requested: RequestedEvidence)
  case UnsupportedUnitPolicy(requested: UnitPolicy)

  def message: String =
    this match
      case InvalidIdentifier(role, value) =>
        s"$role must be non-empty and trimmed, got '$value'"
      case InvalidDescription(role, value) =>
        s"$role must be non-empty and trimmed, got '$value'"
      case InvalidCount(role, value) =>
        s"$role must be positive, got $value"
      case InvalidNonNegativeCount(role, value) =>
        s"$role must be non-negative, got $value"
      case InvalidProbability(role, value, inclusiveZero) =>
        val interval = if inclusiveZero then "[0, 1]" else "(0, 1)"
        s"$role must be finite and in $interval, got $value"
      case InvalidTolerance(role, value) =>
        s"$role must be finite and non-negative, got $value"
      case InvalidComponent(value) =>
        s"component index must be non-negative, got $value"
      case EmptyComponentSet(role) =>
        s"$role must contain at least one component"
      case DuplicateComponent(value) =>
        s"component set contains duplicate index $value"
      case UnorderedComponents(previous, next) =>
        s"component set must be strictly increasing, got $previous before $next"
      case ComponentOutOfRange(value, rank) =>
        s"component index $value is outside fitted rank $rank"
      case InvalidUnit(detail) =>
        s"invalid latent unit: $detail"
      case InvalidPartition(detail) =>
        s"invalid row partition: $detail"
      case RowCountMismatch(role, expected, actual) =>
        s"$role expected $expected rows, got $actual"
      case InvalidSpectrum(detail) =>
        s"invalid ordered spectrum: $detail"
      case InvalidReplicatePlan(detail) =>
        s"invalid replicate plan: $detail"
      case Resampling(error) =>
        s"resampling failed [${error.code.value}]: ${error.message}"
      case NumericalFailure(role, detail) =>
        s"$role failed: $detail"
      case InvalidValidity(detail) =>
        s"invalid validity declaration: $detail"
      case UnsupportedProblem(detail) =>
        s"unsupported inference problem: $detail"
      case FitProtocolMismatch(compiled, protocol) =>
        s"compiled fit '${compiled.label}' does not match protocol fit '${protocol.label}'"
      case TargetProtocolMismatch(compiled, protocol) =>
        s"compiled target '${compiled.value}' does not match protocol target '${protocol.value}'"
      case NullProtocolMismatch(compiled, protocol) =>
        s"compiled null '${compiled.value}' does not match protocol null '${protocol.value}'"
      case ExecutionPolicyConflict(mismatch) =>
        mismatch match
          case ExecutionPolicyMismatch.FixedMonteCarloUnsupported =>
            "the ladder executor currently requires a sequential Monte Carlo policy"
          case ExecutionPolicyMismatch.PerRungDraws =>
            "compiled sequential maxDraws must equal the ladder per-rung allocation"
          case ExecutionPolicyMismatch.Alpha =>
            "compiled sequential alpha must equal the ladder execution alpha"
          case ExecutionPolicyMismatch.BatchSize =>
            "compiled sequential batchSize must equal the ladder execution batch size"
          case ExecutionPolicyMismatch.SequentialBoundary =>
            "compiled sequential boundary must equal the ladder execution boundary"
      case FailedAssumption(id, detail) =>
        s"assumption ${id.value} failed: $detail"
      case NonFiniteStatistic(role, value) =>
        s"$role must be finite, got $value"
      case RankLoss(expected, actual) =>
        s"replicate rank dropped from $expected to $actual"
      case UnitBeyondRank(unit, rank) =>
        s"unit ${unit.value} is beyond fitted rank $rank"
      case ReplicateFailure(replicate, cause) =>
        s"replicate ${replicate.value} failed: ${cause.message}"
      case BudgetExhausted(consumed, allocated) =>
        s"Monte Carlo budget exhausted after $consumed of $allocated draws"
      case UnsupportedEvidence(requested) =>
        s"requested evidence '$requested' is unavailable from the significance ladder"
      case UnsupportedUnitPolicy(requested) =>
        s"unit policy '$requested' forms a multi-axis subspace, but the significance ladder currently implements axis statistics only"

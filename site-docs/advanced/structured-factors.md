# Structured factors

Structured factorization is for component directions that must be sparse,
smooth, or measured under a non-Euclidean geometry. Those choices define the
estimand, so the API asks for them rather than inventing defaults.

## Define the rank-one problem

`RankOneStructuredFactorization.from` needs:

1. an `OpTable` containing the data;
2. a `StructuredFactorGeometry` for rows and another for columns;
3. `HomogeneousFactorPenalties` for each factor; and
4. a `RankOneFactorInitialization` before `solve`.

The row and column parameter identities must be distinct. Each geometry must
match the corresponding data space. The solver returns
`Either[SparseFunctionalFactorizationError, RankOneStructuredFit]`.

The fit contains the row factor, column factor, strength, reconstruction,
stopping status, achieved optimization guarantee, and a certificate with
coordinate residuals. Reaching the iteration limit is reported as an iteration
limit, not as convergence.

## Choose a multi-component estimand

There are two different ways to request more than one component:

- `GreedyStructuredFactorization` fits rank-one components in sequence and
  subtracts each fitted reconstruction.
- `JointRankKStructuredFactorization` fits a simultaneous generalized-Stiefel
  frame.

They are different estimands. Greedy deflation does not stand in for a joint
rank-\(k\) solution. The joint implementation accepts the exact smooth case;
it rejects unsupported nonsmooth joint fitting instead of routing it through
the greedy solver.

Start with rank one when testing a new geometry or penalty. Move to greedy or
joint rank-\(k\) only after deciding which estimand answers the scientific
question.

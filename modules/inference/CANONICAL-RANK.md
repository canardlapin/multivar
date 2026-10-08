# Canonical rank references

The input consists of two full-column-rank matrices with `m > p+q` rows,
already in zero-mean coordinates. `CanonicalRankSpectrum` retains all
`min(p,q)` roots, including zeros. It never centers the rows again. QR and SVD
are Gale operations; no covariance inverse or local eigensolver is used.

Two method identities intentionally express different inferential status:

| Method | Construction | Statistical status |
| --- | --- | --- |
| `canonical-rank-score-orthogonal-permutation/v2` | Complete singular vectors in whitened score coordinates, remove preceding fitted modes, refit every permuted tail | A research comparator; uniform partial-null validity is not established |
| `canonical-rank-gaussian-interlacing-wilks/v1` | Full observed Wilks tails compared with independent Gaussian complete-null reference draws | Conservative under the conditional Gaussian assumptions below |

These arithmetic results do not admit a downstream application, authenticate
independent discovery, or license dependent or non-Gaussian rows.

## Score geometry and ties

For `X=Qx Rx` and `Y=Qy Ry`, take the full thin SVD of `Qx'Qy` and complete
its left and right singular-vector bases orthogonally in their respective QR
coordinates. Multiplication by `Qx` and `Qy` yields orthonormal scores spanning
every original candidate dimension. The earlier Euclidean completion of
canonical coefficients is removed: coefficient orthogonality generally does
not imply score orthogonality. Its historical implementation remains available
at revision `ab811e257dd67f77e8c3b70cb1ea600f274429a3`.

For separated roots, invertible feature changes preserve the corresponding
tail score subspaces and therefore the permuted tail statistics. Canonical
signs and rotations within the appended complement do not affect this result.
If a cut splits adjacent roots within the fixed absolute tolerance `1e-10`,
the permutation API returns `UnidentifiedCanonicalTail` before sampling. It
does not resolve an unidentified eigenspace by an arbitrary SVD orientation.
The observed spectrum remains available, including zero roots. The Gaussian
reference uses only roots and needs no such cut-direction refusal.

`CanonicalRankGeometrySuite` checks all 720 row actions of an independently
computed example, feature shears/scales/rotations/order, block interchange,
equal and unequal dimensions, score orthogonality, and explicit refusals.
The repaired example has H2 count 642/720 in both feature coordinate systems;
the old completion gave 642/720 versus 672/720.

## Conservative Gaussian Wilks reference

Assume fixed candidate spaces and `m` iid zero-mean joint Gaussian rows with
positive-definite marginal covariance. For fixed nuisance design `Z`, this
also holds after applying `Q'` when `Q'Q=I_m`, `Q'Z=0`, and the original
conditional row covariance is `I_n`; here `m=n-rank(Z)`. No second intercept
is subtracted. The method does not establish these assumptions from data.

Let `a=min(p,q)`, `b=max(p,q)`, and consider the null
`rank(SigmaXY)<=s`, where `0<=s<a`. A population change of coordinates on the
smaller block produces `a-s` variables independent of the entire larger
block. QR with those columns first makes their cross-basis matrix `C0` a
row submatrix of the full cross-basis matrix `C`. Singular-value interlacing
gives `r_(s+j)(C)<=r_j(C0)` for `1<=j<=a-s`. Since `-log(1-r^2)` increases,

```
W_(s+1)(C) = -sum(j=1..a-s, log(1-r_(s+j)(C)^2))
           <= full_Wilks(C0).
```

The right side has the complete-null law for independent Gaussian matrices
with `m` rows and dimensions `(a-s,b)`. Population coordinates enter only
the proof; the procedure never estimates them. This extends the rootwise
argument in [Johnstone, Lemma 1 and appendix](https://arxiv.org/pdf/1009.5854)
to the unscaled Wilks sum. Reducing both dimensions, recentering residual
coordinates, or using unequal Bartlett factors is not justified by it.

For `B` independent draws from this dominating reference, the inclusive
probability `(1 + count(Wref >= Wobs))/(B+1)` is super-uniform: couple the
observed statistic below an independent complete-null statistic, then use the
ordinary Monte Carlo rank argument. If true rank is `r`, any false rejection
after prefix-maximum closure implies rejection of the first true null,
`H(r+1)`. Its probability is at most alpha. Dependence between hypotheses
does not invalidate this closure argument.

The claim concerns a rank lower bound in the fixed candidate spaces. It is
not an upper bound on total association rank, and conservative error control
does not promise useful power.

## Reproducibility and numerical boundaries

`GaussianCanonicalRank.run` assigns each hypothesis a child of its caller's
dataset-specific root seed. Each replicate and each block gets another child:
`seed-path/v1` custom domains 1011, 1012 and 1013 respectively. Ordinals and
stream contents do not depend on B or scheduling. Separate datasets must use
separate roots; reusing a random reference table is not an exact deterministic
CDF for a binomial calibration study.

The reference algorithm identity is
`gaussian-canonical-interlacing-splitmix64-boxmuller52/v1`. It uses Resample4s
`Rand` for SplitMix64 words and Box-Muller transforms of open midpoint 52-bit
uniforms, filling each matrix in row-major order. Hypothesis seeds, every
null statistic, completed B, numerical diagnostics and closure are retained.
The finite-sample argument assumes the Gaussian sampling model and exact
arithmetic; floating-point/PRNG execution is checked independently, not
presented as a machine-verified probability theorem.

The QR rank tolerance is `1e-12`; unit or numerically unit correlations with
`1-r^2<=1e-12` refuse. No failed numerical replicate is dropped or redrawn.
Budgets cover planned owned numeric cells and retained receipt cardinality,
not a measured process-heap guarantee. Rejection/admission is unavailable on
failure, and calibration must retain the failure rather than condition on
successful datasets.

`tools/r-parity/generate_canonical_rank_v2.py` supplies independently generated
normal matrices using Python integer arithmetic; base R QR/SVD computes all
reference values. Both platforms check those fixtures. A separate one-column
law check uses `m=6,a=1,b=2`, where `r^2~Beta(1,2)`, `E[W]=1/2`, and
`P(W>=-log(1-.8^2))=.1296`; this distinguishes zero-mean from re-centered rows.

The stepwise refit enum now names score-orthogonal v2 explicitly. Consumers
must record the new method identity and must not relabel historical v1 draws
as results of the repaired algorithm.

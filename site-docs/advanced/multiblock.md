# Multiblock models

Multiblock methods combine matrices whose rows or entities have a declared
relationship. That relationship is a statistical assumption, so Multivar does
not infer it from row positions.

## Choose the row relationship

Multivar keeps three cases separate:

| Row relationship | Main type | Use when |
| --- | --- | --- |
| Shared rows | `AlignedSharedScoreGlrm` | every block measures the same observations |
| Independent rows | `DirectSumStudy` | blocks have separate row spaces |
| Shared entities with different observations | `EntityAlignedStudy` | rows map to a common entity hub |

There is no implicit conversion between these study types.

## Fit covariance-style multiset association

For independent row spaces:

1. wrap each complete block as a `CompleteStudyView`;
2. build a `DirectSumStudy`, which records block offsets and column geometry;
3. declare a `MaximizeAssociation` problem over its row relationship; and
4. call `MultisetAssociation.fit(study, problem, components)`.

The fit returns eigenvalues, feature axes, metric loadings, direct-sum scores,
and scores split back into their original views. The finite feature-space
eigensystem is the explicit dense boundary; callers can reject that
materialization through `StoragePolicy`.

## Fit shared latent scores across aligned blocks

`AlignedSharedScoreGlrm.from` requires one `SharedRowBinding`, a latent space,
and a non-empty vector of `AlignedGlrmBlock` values. Each block supplies its
own observation pattern, feature layout, decoder, loss scaling, and geometry.
All blocks must use the exact shared-row binding.

The resulting program evaluates one shared row-code matrix against every
block. Its fitted encoder can estimate a new code from any declared subset of
blocks, provided those observations identify the latent coordinates.

Use the [GLRM guide](mixed-data.md) first if feature domains, observation
patterns, and fitted encoding are unfamiliar.

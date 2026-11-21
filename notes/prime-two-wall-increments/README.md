# Prime 2 increments

The proof note gives the $O(h)$ bound for the weighted squared wall increment of the prime 2 tail in Masullo’s §12.3. It takes $K\to\infty$ at fixed $X$, then $X\to\infty$ with $h=1/\log X$. Its last section computes the coefficient of $h$ and states the cancellation still needed for $o(h)$. [Proof note](proof.md) and [discussion on #2469](https://github.com/mccorvie/riemann-masullo/issues/2469#issuecomment-7067444368).

## Lean checks

In a checkout of the fork’s `prime-two-wall-increments` branch, run the following from the repository root with Elan installed. This folder pins Lean 4.33.1 and records the mathlib commit and dependencies in its own Lake manifest. The cache command downloads prebuilt mathlib files. [Elan](https://github.com/leanprover/elan) and [dependency revisions](lake-manifest.json).

```sh
cd notes/prime-two-wall-increments
lake exe cache get
lake build
for audit in AuditContinuity AuditIncrements AuditPairSieve AuditDyadicModes AuditKernelWall; do
  lake env lean "$audit.lean" || exit 1
done
lake env leanchecker -v PrimeTwoContinuity PrimeTwoIncrements PrimeTwoPairSieve PrimeTwoDyadicModes PrimeTwoKernelWall PrimeTwoRatios PrimeTwoEta
```

The theorem `integrated_overlap_bound` bounds the double integral of the interval overlap by $hCM^2$ for a finite measure of mass $M$, a weight bounded in absolute value by $C$ and a complex factor of modulus at most one. In the proof note the measure is the total variation of $\mu_X$ and the complex factor is the product of its phase at one point and the conjugate phase at the other. Theorems `weighted_pair_limit` and `weighted_pair_zero_iff` check the limiting overlap sum for an absolutely summable family of real jumps and characterise when that sum vanishes. [Lean source](PrimeTwoIncrements.lean) and [theorem and axiom checks](AuditIncrements.lean).

The other checks cover the one-sided limits of an exterior sum, finite Fourier orthogonality for ratio pairs, the dyadic band identities, the three-wall profile and the finite matrix commutator. [Continuity](PrimeTwoContinuity.lean), [ratio pairs](PrimeTwoPairSieve.lean), [dyadic bands](PrimeTwoDyadicModes.lean) and [three walls](PrimeTwoKernelWall.lean).

The analytic work identifying $\mu_X$, bounding its total variation uniformly in $X$ and passing from the heat measure to the zeta line sums is in the proof note. Those steps specify how the Lean lemmas apply to the original tail. [Area estimate](prime-two-denominators.md#every-wall), [two-scale limit](prime-two-increments.md#two-scales) and [remaining line sums](prime-two-kernel-wall.md).

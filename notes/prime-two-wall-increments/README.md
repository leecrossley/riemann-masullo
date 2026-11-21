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
lake env leanchecker -v PrimeTwoContinuity PrimeTwoIncrements PrimeTwoWallEnergy PrimeTwoPairSieve PrimeTwoDyadicModes PrimeTwoKernelWall PrimeTwoRatios PrimeTwoEta
```

The theorem `wall_energy_identity` expands the weighted square of the complex wall integral and proves the change of order that gives the interval overlap. The theorem `density_wall_energy_bound` applies this to an integrable complex density and gives $hCM^2$, where $C$ bounds the weight and $M$ is the integral of the absolute density. For the tail, `tail_energy_isBigO` proves $O(h)$ from a uniform bound on $M$ and a bound $Bh$ on the change contributed by the removed denominators wherever the weight is positive. [Wall energy proofs](PrimeTwoWallEnergy.lean) and [theorem and axiom checks](AuditIncrements.lean).

For an absolutely summable family of real jumps, `weighted_pair_limit` computes the limiting overlap sum. The theorem `weighted_pair_zero_iff` characterises when it vanishes. [Jump lemmas](PrimeTwoIncrements.lean).

The other checks cover the one-sided limits of an exterior sum, finite Fourier orthogonality for ratio pairs, the dyadic band identities, the three-wall profile and the finite matrix commutator. [Continuity](PrimeTwoContinuity.lean), [ratio pairs](PrimeTwoPairSieve.lean), [dyadic bands](PrimeTwoDyadicModes.lean) and [three walls](PrimeTwoKernelWall.lean).

The proof note supplies the inputs for Masullo’s tail after $K\to\infty$ by identifying the heat density and bounding its total variation and the contribution of the removed block. It also identifies the coefficient of $h$ with the zeta line sums. [Area estimate](prime-two-denominators.md#every-wall), [denominator estimate](prime-two-increments.md#distinct-ratios), [two-scale limit](prime-two-increments.md#two-scales) and [remaining line sums](prime-two-kernel-wall.md).

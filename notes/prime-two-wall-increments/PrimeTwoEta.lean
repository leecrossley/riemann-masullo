import Std

namespace PrimeTwoEta

open Lean.Grind

variable {R : Type u} [CommRing R]

def geometricSum (v : R) : Nat → R
  | 0 => 0
  | n + 1 => geometricSum v n + v ^ n

theorem geometric_identity (v : R) (n : Nat) :
    (1 - v) * geometricSum v n = 1 - v ^ n := by
  induction n with
  | zero =>
    simp only [geometricSum, Semiring.pow_zero, Semiring.mul_zero]
    grind
  | succ n ih =>
    simp only [geometricSum, Semiring.pow_succ]
    grind

def reciprocal (a c zeta : R) (n : Nat) : R :=
  a * c * geometricSum (1 - c * (a * zeta)) n

theorem reciprocal_error (a c zeta : R) (n : Nat) :
    1 - zeta * reciprocal a c zeta n = (1 - c * (a * zeta)) ^ n := by
  unfold reciprocal
  have h := geometric_identity (1 - c * (a * zeta)) n
  grind

theorem residue_multiplier_at_zero (a c zeta : R) (n : Nat) (h : zeta = 0) :
    (1 - c * (a * zeta)) ^ n = 1 := by
  rw [← reciprocal_error, h, Semiring.zero_mul]
  grind

def factorial : Nat → Nat
  | 0 => 1
  | n + 1 => (n + 1) * factorial n

set_option maxRecDepth 65536 in
set_option maxHeartbeats 10000000 in
theorem tail_integer_bound :
    52 * 10 ^ 19935 * 5120 ^ 24000 * (factorial 24000) ^ 2 < (3 * 10 ^ 12) ^ 24000 := by
  decide +kernel

end PrimeTwoEta

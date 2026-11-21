import Std

namespace PrimeTwoRatios

open Lean.Grind

attribute [local instance] Semiring.natCast

variable {R : Type u} [CommRing R]

def sumTo (f : Nat → R) : Nat → R
  | 0 => 0
  | n + 1 => sumTo f n + f (n + 1)

def convolution (f g : Nat → R) (n : Nat) : R :=
  sumTo (fun d => if n % d = 0 then f d * g (n / d) else 0) n

def convolutionPower (f : Nat → R) : Nat → Nat → R
  | 0 => fun n => if n = 1 then 1 else 0
  | j + 1 => convolution f (convolutionPower f j)

private theorem convolution_one (f g : Nat → R) :
    convolution f g 1 = f 1 * g 1 := by
  change 0 + f 1 * g 1 = _
  grind

private theorem convolution_three (f g : Nat → R) :
    convolution f g 3 = f 1 * g 3 + f 3 * g 1 := by
  change 0 + f 1 * g 3 + 0 + f 3 * g 1 = _
  grind

private theorem convolution_seven (f g : Nat → R) :
    convolution f g 7 = f 1 * g 7 + f 7 * g 1 := by
  change 0 + f 1 * g 7 + 0 + 0 + 0 + 0 + 0 + f 7 * g 1 = _
  grind

private theorem power_one (f : Nat → R) (h : f 1 = 1) (j : Nat) :
    convolutionPower f j 1 = 1 := by
  induction j with
  | zero => rfl
  | succ j ih =>
    rw [convolutionPower, convolution_one, h, ih, Semiring.one_mul]

private theorem power_three (f : Nat → R) (h : f 1 = 1) (j : Nat) :
    convolutionPower f j 3 = (j : R) * f 3 := by
  induction j with
  | zero =>
    change 0 = (Nat.cast 0 : R) * f 3
    rw [Semiring.natCast_zero, Semiring.zero_mul]
  | succ j ih =>
    rw [convolutionPower, convolution_three, h, ih, power_one f h, Semiring.natCast_succ]
    grind

theorem seven_three_moment (eta derivative : Nat → R) (logSeven : R)
    (one : eta 1 = 1) (three : eta 3 = 1)
    (derivOne : derivative 1 = 0) (derivSeven : derivative 7 = -logSeven) (j : Nat) :
    convolution derivative (convolutionPower eta j) 7 * convolutionPower eta (j + 1) 3 =
      -(((j + 1 : Nat) : R) * logSeven) := by
  rw [convolution_seven, derivOne, derivSeven, power_one eta one, power_three eta one, three]
  grind

def inPulseWindow (n d : Nat) : Prop :=
  (3 * d < 2 * n ∧ 2 * n < 5 * d) ∨ (2 * d < 5 * n ∧ 3 * n < 2 * d)

theorem supported_integer_or_reciprocal {n d : Nat} (h : inPulseWindow n d)
    (reduced : d = 1 ∨ n = 1) :
    (n = 2 ∧ d = 1) ∨ (n = 1 ∧ d = 2) := by
  unfold inPulseWindow at h
  omega

theorem aggregate_lower_bound_arithmetic :
    544 * 600 < (600 - 544 : Nat) * 10 ^ 20000 := by
  decide +kernel

theorem regularised_pole_split (e derivative d c heat phi : R) (h : e * phi = heat) :
    heat * (d + c * (derivative - d * e)) + (1 - c * e) * (derivative * phi - d * heat) =
      derivative * phi := by
  rw [← h]
  grind

theorem dyadic_density_cutoff (v : R) (j : Nat) :
    -v + sumTo (fun k => (1 - v) * v ^ k) j = -(v ^ (j + 1)) := by
  induction j with
  | zero =>
    simp only [sumTo, Semiring.pow_succ, Semiring.pow_zero, Semiring.one_mul]
    grind
  | succ j ih =>
    change -v + (sumTo _ j + (1 - v) * v ^ (j + 1)) = -(v ^ ((j + 1) + 1))
    rw [Semiring.pow_succ v (j + 1)]
    grind

theorem dyadic_density_slope (a v : R) (j : Nat) :
    a * v + sumTo (fun k => (a * v - (k : R) * a * (1 - v)) * v ^ k) j =
      ((j + 1 : Nat) : R) * a * v ^ (j + 1) := by
  induction j with
  | zero =>
    simp only [sumTo, Semiring.pow_succ, Semiring.pow_zero, Semiring.one_mul,
      Semiring.natCast_succ, Semiring.natCast_zero]
    grind
  | succ j ih =>
    change a * v + (sumTo _ j + (a * v - ((j + 1 : Nat) : R) * a * (1 - v)) *
      v ^ (j + 1)) = (((j + 1) + 1 : Nat) : R) * a * v ^ ((j + 1) + 1)
    rw [Semiring.pow_succ v (j + 1), Semiring.natCast_succ (j + 1)]
    grind

end PrimeTwoRatios

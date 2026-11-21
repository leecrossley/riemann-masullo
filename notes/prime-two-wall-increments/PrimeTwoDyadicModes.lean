import Mathlib.Analysis.Normed.Group.InfiniteSum
import Mathlib.Analysis.Complex.Basic
import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring

open Finset
open scoped ComplexConjugate

namespace PrimeTwoDyadicModes

def phase (z : ℂ) (n : ℕ) : ℂ := z ^ padicValNat 2 n

def unsignedSum (z : ℂ) (n : ℕ) : ℂ :=
  ∑ j ∈ range n, phase z (j + 1)

def alternatingSum (z : ℂ) (n : ℕ) : ℂ :=
  ∑ j ∈ range n, (-1 : ℂ) ^ j * phase z (j + 1)

private theorem phase_odd (z : ℂ) (j : ℕ) : phase z (2 * j + 1) = 1 := by
  have h : ¬2 ∣ 2 * j + 1 := by omega
  simp [phase, padicValNat.eq_zero_of_not_dvd h]

private theorem phase_double (z : ℂ) (j : ℕ) :
    phase z (2 * (j + 1)) = z * phase z (j + 1) := by
  simp [phase, padicValNat.mul (p := 2) (by decide : 2 ≠ 0) (Nat.succ_ne_zero j),
    pow_add]

theorem unsigned_double (z : ℂ) (n : ℕ) :
    unsignedSum z (2 * n) = n + z * unsignedSum z n := by
  induction n with
  | zero => simp [unsignedSum]
  | succ n ih =>
    change unsignedSum z (2 * n + 2) = _
    have step : unsignedSum z (2 * n + 2) =
        unsignedSum z (2 * n) + phase z (2 * n + 1) + phase z (2 * (n + 1)) := by
      simp [unsignedSum, sum_range_succ, Nat.mul_add]
    rw [step, ih, phase_odd, phase_double]
    simp [unsignedSum, sum_range_succ, mul_add]
    ring

theorem alternating_double (z : ℂ) (n : ℕ) :
    alternatingSum z (2 * n) = n - z * unsignedSum z n := by
  induction n with
  | zero => simp [alternatingSum, unsignedSum]
  | succ n ih =>
    change alternatingSum z (2 * n + 2) = _
    have step : alternatingSum z (2 * n + 2) = alternatingSum z (2 * n) +
        (-1 : ℂ) ^ (2 * n) * phase z (2 * n + 1) +
        (-1 : ℂ) ^ (2 * n + 1) * phase z (2 * (n + 1)) := by
      simp [alternatingSum, sum_range_succ, Nat.mul_add]
    rw [step, ih, phase_odd, phase_double]
    simp [unsignedSum, sum_range_succ, pow_add, pow_mul, mul_add]
    ring

theorem unsigned_block (z : ℂ) (m : ℕ) :
    (2 - z) * unsignedSum z (2 ^ m) = 2 ^ m + (1 - z) * z ^ m := by
  induction m with
  | zero =>
    simp [unsignedSum, phase]
    ring
  | succ m ih =>
    rw [pow_succ' (2 : ℕ), unsigned_double, Nat.cast_pow, Nat.cast_ofNat, pow_succ, pow_succ]
    linear_combination z * ih

theorem alternating_block (z : ℂ) (m : ℕ) :
    (2 - z) * alternatingSum z (2 ^ (m + 1)) =
      (1 - z) * (2 ^ (m + 1) - z ^ (m + 1)) := by
  rw [pow_succ' (2 : ℕ), alternating_double, Nat.cast_pow, Nat.cast_ofNat, pow_succ, pow_succ]
  linear_combination -z * unsigned_block z m

theorem split_square (right left low : ℂ) :
    (right - left - low) * conj (right - left - low) =
      low * conj low + (right * conj right - right * conj low - low * conj right) +
      (left * conj left - right * conj left - left * conj right + left * conj low +
        low * conj left) := by
  simp only [map_sub]
  ring

theorem band_sum_error (u v : Finset ℤ) (huv : Disjoint u v) (a b t : ℤ → ℂ)
    (ha : Summable fun k => ‖a k‖) (hb : Summable fun k => ‖b k‖)
    (ht : ∀ k ∈ u, t k = 0) :
    ‖(∑ k ∈ u, (a k + b k + t k)) - ∑' k, a k‖ ≤
      (∑' k : {j // j ∉ v}, ‖b k‖) + ∑' k : {j // j ∉ u}, ‖a k‖ := by
  have other : (∑ k ∈ u, ‖b k‖) ≤ ∑' k : {j // j ∉ v}, ‖b k‖ := by
    have bound := hb.sum_le_tsum (u ∪ v) (fun k _ => norm_nonneg (b k))
    rw [sum_union huv, ← hb.sum_add_tsum_subtype_compl v] at bound
    linarith
  have split : (∑ k ∈ u, (a k + b k + t k)) - ∑' k, a k =
      (∑ k ∈ u, b k) - ∑' k : {j // j ∉ u}, a k := by
    rw [sum_add_distrib, sum_add_distrib, sum_eq_zero ht, add_zero,
      ← ha.of_norm.sum_add_tsum_subtype_compl u]
    abel
  rw [split]
  exact (norm_sub_le _ _).trans (add_le_add ((norm_sum_le u b).trans other)
    (norm_tsum_le_tsum_norm (ha.subtype _)))

theorem half_jump (before jump : ℂ) :
    (‖before + jump‖ ^ 2 - ‖before‖ ^ 2) / 2 =
      (conj before * jump).re + ‖jump‖ ^ 2 / 2 := by
  rw [add_comm before jump]
  simp only [Complex.sq_norm, Complex.normSq_add]
  rw [mul_comm jump (conj before)]
  ring

end PrimeTwoDyadicModes

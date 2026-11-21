import Mathlib.Analysis.Fourier.AddCircleMulti
import Mathlib.Tactic.Linarith

open Finset
open scoped ComplexConjugate InnerProductSpace

namespace PrimeTwoPairSieve

theorem matching_exponents {ι P : Type*} [Fintype P] (s : Finset ι)
    (v : ι → P → ℤ) (a : ι → ℂ) :
    (‖∑ i ∈ s, a i • UnitAddTorus.mFourierLp 2 (v i)‖ : ℂ) ^ 2 =
      ∑ i ∈ s, ∑ j ∈ s, if v i = v j then conj (a i) * a j else 0 := by
  classical
  have hnorm := inner_self_eq_norm_sq_to_K (𝕜 := ℂ)
    (∑ i ∈ s, a i • UnitAddTorus.mFourierLp 2 (v i))
  refine hnorm.symm.trans ?_
  rw [sum_inner]
  simp only [inner_sum, inner_smul_left, inner_smul_right]
  apply sum_congr rfl
  intro i _
  apply sum_congr rfl
  intro j _
  rw [(orthonormal_iff_ite.mp (UnitAddTorus.orthonormal_mFourier (d := P))) (v i) (v j)]
  split <;> simp_all [mul_comm]

theorem one_quotient_bound {G : Type*} [Group G] [DecidableEq G]
    (s : Finset G) (a : G → ℂ) (u : G) :
    (∑ r ∈ s.filter (fun r => u * r ∈ s), ‖a r‖ * ‖a (u * r)‖) ≤
      ∑ r ∈ s, ‖a r‖ ^ 2 := by
  let t := s.filter (fun r => u * r ∈ s)
  have image_subset : t.image (fun r => u * r) ⊆ s := by
    intro r hr
    obtain ⟨r', hr', rfl⟩ := mem_image.mp hr
    exact (mem_filter.mp hr').2
  have left : (∑ r ∈ t, ‖a r‖ ^ 2) ≤ ∑ r ∈ s, ‖a r‖ ^ 2 :=
    sum_le_sum_of_subset_of_nonneg (filter_subset _ _) (fun _ _ _ => sq_nonneg _)
  have right : (∑ r ∈ t, ‖a (u * r)‖ ^ 2) ≤ ∑ r ∈ s, ‖a r‖ ^ 2 := by
    calc
      _ = ∑ r ∈ t.image (fun r => u * r), ‖a r‖ ^ 2 :=
        (sum_image (mul_right_injective u).injOn).symm
      _ ≤ _ := sum_le_sum_of_subset_of_nonneg image_subset (fun _ _ _ => sq_nonneg _)
  have twice : 2 * (∑ r ∈ t, ‖a r‖ * ‖a (u * r)‖) ≤
      (∑ r ∈ t, ‖a r‖ ^ 2) + ∑ r ∈ t, ‖a (u * r)‖ ^ 2 := by
    rw [mul_sum, ← sum_add_distrib]
    apply sum_le_sum
    intro r _
    nlinarith [sq_nonneg (‖a r‖ - ‖a (u * r)‖)]
  change (∑ r ∈ t, ‖a r‖ * ‖a (u * r)‖) ≤ _
  linarith

theorem quotient_set_bound {G : Type*} [Group G] [DecidableEq G]
    (s U : Finset G) (a : G → ℂ) :
    (∑ u ∈ U, ∑ r ∈ s.filter (fun r => u * r ∈ s), ‖a r‖ * ‖a (u * r)‖) ≤
      U.card * ∑ r ∈ s, ‖a r‖ ^ 2 := by
  calc
    _ ≤ ∑ _u ∈ U, ∑ r ∈ s, ‖a r‖ ^ 2 :=
      sum_le_sum fun u _ => one_quotient_bound s a u
    _ = _ := by simp

end PrimeTwoPairSieve

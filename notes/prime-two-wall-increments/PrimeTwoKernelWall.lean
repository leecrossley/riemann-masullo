import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
import Mathlib.Data.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic.Ring

open Matrix

namespace PrimeTwoKernelWall

noncomputable def profile (u δ t : ℝ) : ℝ := Real.sinh (t * u) / Real.sinh (δ * u)

theorem profile_endpoints {u δ : ℝ} (hu : 0 < u) (hδ : 0 < δ) :
    profile u δ 0 = 0 ∧ profile u δ δ = 1 := by
  have h : Real.sinh (δ * u) ≠ 0 := Real.sinh_ne_zero.mpr (mul_pos hδ hu).ne'
  simp [profile, h]

theorem profile_derivative (u δ t : ℝ) :
    HasDerivAt (profile u δ) (Real.cosh (t * u) * u / Real.sinh (δ * u)) t := by
  exact (((Real.hasDerivAt_sinh (t * u)).comp t ((hasDerivAt_id t).mul_const u)).div_const
    (Real.sinh (δ * u))) |>.congr_deriv (by simp)

theorem profile_curvature (u δ t : ℝ) :
    HasDerivAt (fun s => Real.cosh (s * u) * u / Real.sinh (δ * u))
      (u ^ 2 * profile u δ t) t := by
  have h := ((((Real.hasDerivAt_cosh (t * u)).comp t ((hasDerivAt_id t).mul_const u)).mul_const u).div_const
    (Real.sinh (δ * u)))
  apply h.congr_deriv
  simp [profile]
  ring

theorem compressed_commutator {ι : Type*} [Fintype ι] (a n : Matrix ι ι ℂ)
    (b : Matrix ι Unit ℂ) (l : Matrix Unit ι ℂ) (d e : Matrix Unit Unit ℂ) :
    ((fromBlocks a b 0 d * fromBlocks n 0 l e -
      fromBlocks n 0 l e * fromBlocks a b 0 d).toBlocks₁₁) =
      a * n - n * a + b * l := by
  classical
  rw [fromBlocks_multiply, fromBlocks_multiply]
  ext i j
  simp [toBlocks₁₁, fromBlocks]
  ring

theorem top_jet_trace (m : ℕ) (a n : Matrix (Fin (m + 1)) (Fin (m + 1)) ℂ)
    (b : Matrix (Fin (m + 1)) Unit ℂ) :
    (a * n - n * a + b * single () (Fin.last m) (m + 1 : ℂ)).trace =
      (m + 1 : ℂ) * b (Fin.last m) () := by
  classical
  rw [trace_add, trace_sub, trace_mul_comm a n, sub_self, zero_add, trace_mul_comm b,
    trace_single_mul]
  rfl

end PrimeTwoKernelWall

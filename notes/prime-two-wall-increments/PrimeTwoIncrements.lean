import Mathlib.Analysis.Normed.Group.Tannery
import Mathlib.Analysis.Normed.Ring.InfiniteSum
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

open Filter MeasureTheory Set Topology

namespace PrimeTwoIncrements

noncomputable def overlap (ω : ℝ → ℝ) (h a b : ℝ) : ℝ :=
  let l := max a b - h / 2
  let u := min a b + h / 2
  if l ≤ u then (∫ x in l..u, ω x) / h else 0

theorem overlap_bound {ω : ℝ → ℝ} {C : ℝ} (hC : ∀ x, ‖ω x‖ ≤ C)
    {h : ℝ} (hh : 0 < h) (a b : ℝ) : ‖overlap ω h a b‖ ≤ C := by
  have hC₀ := (norm_nonneg (ω 0)).trans (hC 0)
  dsimp only [overlap]
  split
  · next hlu =>
    rw [norm_div, Real.norm_eq_abs h, abs_of_pos hh, div_le_iff₀ hh]
    refine (intervalIntegral.norm_integral_le_of_norm_le_const (fun x _ => hC x)).trans ?_
    rw [abs_of_nonneg (sub_nonneg.mpr hlu)]
    apply mul_le_mul_of_nonneg_left _ hC₀
    have hab : min a b ≤ max a b := (min_le_left a b).trans (le_max_left a b)
    linarith
  · simpa using hC₀

theorem integrated_overlap_bound {α : Type*} [MeasurableSpace α] (ν : Measure α)
    [IsFiniteMeasure ν] (position : α → ℝ) (phase : α → α → ℂ)
    (hphase : ∀ x y, ‖phase x y‖ ≤ 1) {ω : ℝ → ℝ} {C : ℝ} (hC : ∀ s, ‖ω s‖ ≤ C)
    {h : ℝ} (hh : 0 < h) :
    ‖h • (∫ x, ∫ y, overlap ω h (position x) (position y) • phase x y ∂ν ∂ν)‖ ≤
      h * C * (ν.real univ) ^ 2 := by
  have pointwise (x y : α) :
      ‖overlap ω h (position x) (position y) • phase x y‖ ≤ C := by
    rw [norm_smul]
    calc
      _ ≤ ‖overlap ω h (position x) (position y)‖ * 1 :=
        mul_le_mul_of_nonneg_left (hphase x y) (norm_nonneg _)
      _ ≤ C := by simpa using overlap_bound hC hh (position x) (position y)
  have inner (x : α) :
      ‖∫ y, overlap ω h (position x) (position y) • phase x y ∂ν‖ ≤ C * ν.real univ :=
    norm_integral_le_of_norm_le_const (Eventually.of_forall (pointwise x))
  rw [norm_smul, Real.norm_eq_abs, abs_of_pos hh]
  calc
    _ ≤ h * ((C * ν.real univ) * ν.real univ) :=
      mul_le_mul_of_nonneg_left
        (norm_integral_le_of_norm_le_const (Eventually.of_forall inner)) hh.le
    _ = _ := by ring

theorem overlap_off {ω : ℝ → ℝ} {a b : ℝ} (hab : a ≠ b) :
    Tendsto (fun h => overlap ω h a b) (𝓝[>] 0) (𝓝 0) := by
  have hgap : 0 < max a b - min a b := by
    rcases lt_or_gt_of_ne hab with h | h <;> simp [h.le, h]
  apply tendsto_const_nhds.congr'
  filter_upwards [mem_nhdsWithin_of_mem_nhds (Iio_mem_nhds hgap)] with h hsmall
  have hlt : h < max a b - min a b := hsmall
  have hempty : ¬ max a b - h / 2 ≤ min a b + h / 2 := by linarith
  simp [overlap, hempty]

theorem overlap_on {ω : ℝ → ℝ} (hω : Continuous ω) (a : ℝ) :
    Tendsto (fun h => overlap ω h a a) (𝓝[>] 0) (𝓝 (ω a)) := by
  rw [Metric.tendsto_nhds]
  intro ε hε
  obtain ⟨δ, hδ, hball⟩ := (Metric.continuousAt_iff.mp (hω.continuousAt (x := a)))
    (ε / 2) (half_pos hε)
  filter_upwards [self_mem_nhdsWithin, mem_nhdsWithin_of_mem_nhds (Iio_mem_nhds hδ)] with h hp hs
  have hh : 0 < h := hp
  have hsmall : h < δ := hs
  let l := a - h / 2
  let u := a + h / 2
  have hlu : l ≤ u := by dsimp [l, u]; linarith
  have hlen : u - l = h := by dsimp [l, u]; ring
  have hbound : ∀ x ∈ uIoc l u, ‖ω x - ω a‖ ≤ ε / 2 := by
    intro x hx
    rw [uIoc_of_le hlu] at hx
    apply le_of_lt
    apply hball
    rw [Real.dist_eq, abs_lt]
    dsimp [l, u] at hx
    constructor <;> linarith [hx.1, hx.2]
  have hmean : (∫ x in l..u, ω x - ω a) = (∫ x in l..u, ω x) - h * ω a := by
    rw [intervalIntegral.integral_sub (hω.intervalIntegrable _ _) intervalIntegrable_const,
      intervalIntegral.integral_const, hlen, smul_eq_mul]
  have estimate := intervalIntegral.norm_integral_le_of_norm_le_const hbound
  rw [hlen, abs_of_pos hh] at estimate
  have hdiff : (∫ x in l..u, ω x) / h - ω a = (∫ x in l..u, ω x - ω a) / h := by
    rw [hmean, sub_div, mul_div_cancel_left₀ _ hh.ne']
  simp only [overlap, max_self, min_self]
  change dist (if l ≤ u then (∫ x in l..u, ω x) / h else 0) (ω a) < ε
  rw [if_pos hlu, dist_eq_norm, hdiff, norm_div, Real.norm_eq_abs h, abs_of_pos hh]
  exact ((div_le_iff₀ hh).2 estimate).trans_lt (half_lt_self hε)

theorem weighted_pair_limit {ι : Type*} (p : ι → ℝ) (hp : Function.Injective p)
    (w : ι → ℝ) (hw : Summable (fun i => ‖w i‖)) {ω : ℝ → ℝ} (hω : Continuous ω)
    {C : ℝ} (hC : ∀ x, ‖ω x‖ ≤ C) :
    Tendsto (fun h => ∑' t : ι × ι, w t.1 * w t.2 * overlap ω h (p t.1) (p t.2))
      (𝓝[>] 0) (𝓝 (∑' i, w i ^ 2 * ω (p i))) := by
  classical
  let v (t : ι × ι) := if t.1 = t.2 then w t.1 ^ 2 * ω (p t.1) else 0
  let bound (t : ι × ι) := ‖w t.1 * w t.2‖ * C
  have hb : Summable bound := (hw.mul_norm hw).mul_right C
  have major {h : ℝ} (hh : 0 < h) (t : ι × ι) :
      ‖w t.1 * w t.2 * overlap ω h (p t.1) (p t.2)‖ ≤ bound t := by
    rw [norm_mul]
    exact mul_le_mul_of_nonneg_left (overlap_bound hC hh _ _) (norm_nonneg _)
  have limit (t : ι × ι) : Tendsto
      (fun h => w t.1 * w t.2 * overlap ω h (p t.1) (p t.2)) (𝓝[>] 0) (𝓝 (v t)) := by
    obtain ⟨i, j⟩ := t
    by_cases hij : i = j
    · subst j
      simpa [v, pow_two] using (overlap_on hω (p i)).const_mul (w i * w i)
    · simpa [v, hij] using (overlap_off (ω := ω) (hp.ne hij)).const_mul (w i * w j)
  have positive : ∀ᶠ h : ℝ in 𝓝[>] 0, 0 < h := self_mem_nhdsWithin
  have hv : Summable v := hb.of_norm_bounded fun t =>
    le_of_tendsto (tendsto_norm.comp (limit t)) (positive.mono fun _ hh => major hh t)
  have regroup : (∑' t, v t) = ∑' i, w i ^ 2 * ω (p i) := by
    rw [hv.tsum_prod]
    apply tsum_congr
    intro i
    rw [tsum_eq_single i (fun j hji => by simp [v, Ne.symm hji])]
    simp [v]
  rw [← regroup]
  exact tendsto_tsum_of_dominated_convergence hb limit
    (positive.mono fun _ hh t => major hh t)

theorem weighted_pair_zero_iff {ι : Type*} (p : ι → ℝ) (hp : Function.Injective p)
    (w : ι → ℝ) (hw : Summable (fun i => ‖w i‖)) {ω : ℝ → ℝ} (hω : Continuous ω)
    (hpos : ∀ x, 0 ≤ ω x) {C : ℝ} (hC : ∀ x, ‖ω x‖ ≤ C) :
    Tendsto (fun h => ∑' t : ι × ι, w t.1 * w t.2 * overlap ω h (p t.1) (p t.2))
      (𝓝[>] 0) (𝓝 0) ↔ ∀ i, 0 < ω (p i) → w i = 0 := by
  have diag : Function.Injective (fun i : ι => (i, i)) := fun _ _ heq => congrArg Prod.fst heq
  have bounds : Summable (fun i => ‖w i * w i‖ * C) :=
    ((hw.mul_norm hw).comp_injective diag).mul_right C
  have pos (i : ι) : 0 ≤ w i ^ 2 * ω (p i) := mul_nonneg (sq_nonneg _) (hpos _)
  have sum : Summable (fun i => w i ^ 2 * ω (p i)) := by
    apply bounds.of_norm_bounded
    intro i
    rw [norm_mul, pow_two]
    exact mul_le_mul_of_nonneg_left (hC _) (norm_nonneg _)
  have limit := weighted_pair_limit p hp w hw hω hC
  constructor
  · intro hz i hi
    have heq := tendsto_nhds_unique limit hz
    have hle := sum.le_tsum i (fun j _ => pos j)
    rw [heq] at hle
    exact sq_eq_zero_iff.mp ((mul_eq_zero.mp (le_antisymm hle (pos i))).resolve_right hi.ne')
  · intro hz
    have zero (i : ι) : w i ^ 2 * ω (p i) = 0 := by
      rcases (hpos (p i)).eq_or_lt with hi | hi
      · rw [← hi, mul_zero]
      · simp [hz i hi]
    simpa only [zero, tsum_zero] using limit

end PrimeTwoIncrements

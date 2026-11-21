import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.Normed.Group.Tannery
import Mathlib.Topology.Order.LeftRight

open Filter Set Topology

namespace PrimeTwoContinuity

variable {ι : Type*}

noncomputable def step (p : ℝ) (w : ℂ) (s : ℝ) : ℂ :=
  if p < s then w else if p = s then w / 2 else 0

noncomputable def wall (p : ι → ℝ) (w : ι → ℂ) (s : ℝ) : ℂ :=
  ∑' i, step (p i) (w i) s

noncomputable def below (p : ι → ℝ) (w : ι → ℂ) (b : ℝ) : ℂ :=
  ∑' i, if p i < b then w i else 0

noncomputable def onWall (p : ι → ℝ) (w : ι → ℂ) (b : ℝ) : ℂ :=
  ∑' i, if p i = b then w i else 0

noncomputable def exteriorTail (p : ι → ℝ) (w : ι → ℂ) (pole : ℂ) (s : ℝ) : ℂ :=
  pole - ∑' i, step s (w i) (p i)

private theorem step_bound (p s : ℝ) (w : ℂ) : ‖step p w s‖ ≤ ‖w‖ := by
  rcases lt_trichotomy p s with h | h | h
  · simp [step, h]
  · simp [step, h]
  · simp [step, not_lt_of_ge h.le, ne_of_gt h]

private theorem selected_summable (w : ι → ℂ) (hw : Summable (fun i => ‖w i‖))
    (P : ι → Prop) [DecidablePred P] : Summable (fun i => if P i then w i else 0) := by
  apply Summable.of_norm_bounded hw
  intro i
  by_cases h : P i <;> simp [h]

private theorem fibre_sum (p : ι → ℝ) (w : ι → ℂ) (x : ℝ) :
    (∑' j : p ⁻¹' {x}, w j) = onWall p w x := by
  classical
  rw [tsum_subtype]
  unfold onWall
  apply tsum_congr
  intro i
  by_cases hi : p i = x <;> simp [Set.indicator, hi]

private theorem step_left (p b : ℝ) (w : ℂ) :
    Tendsto (step p w) (𝓝[<] b) (𝓝 (if p < b then w else 0)) := by
  by_cases h : p < b
  · rw [if_pos h]
    apply tendsto_const_nhds.congr'
    filter_upwards [mem_nhdsWithin_of_mem_nhds (Ioi_mem_nhds h)] with s hs
    simp [step, (show p < s from hs)]
  · rw [if_neg h]
    apply tendsto_const_nhds.congr'
    filter_upwards [self_mem_nhdsWithin] with s hs
    have hsp : s < p := lt_of_lt_of_le (show s < b from hs) (le_of_not_gt h)
    simp [step, not_lt_of_ge hsp.le, ne_of_gt hsp]

private theorem step_right (p b : ℝ) (w : ℂ) :
    Tendsto (step p w) (𝓝[>] b) (𝓝 (if p ≤ b then w else 0)) := by
  by_cases h : p ≤ b
  · rw [if_pos h]
    apply tendsto_const_nhds.congr'
    filter_upwards [self_mem_nhdsWithin] with s hs
    simp [step, lt_of_le_of_lt h (show b < s from hs)]
  · rw [if_neg h]
    apply tendsto_const_nhds.congr'
    filter_upwards [mem_nhdsWithin_of_mem_nhds (Iio_mem_nhds (lt_of_not_ge h))] with s hs
    have hsp : s < p := hs
    simp [step, not_lt_of_ge hsp.le, ne_of_gt hsp]

variable (p : ι → ℝ) (w : ι → ℂ) (hw : Summable (fun i => ‖w i‖)) (b : ℝ)

include hw

theorem left_limit : Tendsto (wall p w) (𝓝[<] b) (𝓝 (below p w b)) := by
  exact tendsto_tsum_of_dominated_convergence hw (fun i => step_left (p i) b (w i))
    (.of_forall fun s i => step_bound (p i) s (w i))

theorem right_limit : Tendsto (wall p w) (𝓝[>] b) (𝓝 (below p w b + onWall p w b)) := by
  unfold below onWall
  rw [← (selected_summable w hw (fun i => p i < b)).tsum_add
    (selected_summable w hw (fun i => p i = b))]
  apply tendsto_tsum_of_dominated_convergence hw _ (.of_forall fun s i => step_bound (p i) s (w i))
  intro i
  have h := step_right (p i) b (w i)
  rcases lt_trichotomy (p i) b with hi | hi | hi
  · simpa [hi, hi.le, ne_of_lt hi] using h
  · simpa [hi] using h
  · simpa [not_lt_of_ge hi.le, ne_of_gt hi, not_le_of_gt hi] using h

theorem midpoint_value : wall p w b = below p w b + onWall p w b / 2 := by
  unfold wall below onWall
  rw [← tsum_div_const, ← (selected_summable w hw (fun i => p i < b)).tsum_add
    ((selected_summable w hw (fun i => p i = b)).div_const 2)]
  apply tsum_congr
  intro i
  rcases lt_trichotomy (p i) b with hi | hi | hi
  · simp [step, hi, ne_of_lt hi]
  · simp [step, hi]
  · simp [step, not_lt_of_ge hi.le, ne_of_gt hi]

theorem continuousAt_iff_onWall_zero : ContinuousAt (wall p w) b ↔ onWall p w b = 0 := by
  have left := left_limit p w hw b
  have right := right_limit p w hw b
  constructor
  · intro hc
    have hl := tendsto_nhds_unique left (hc.continuousWithinAt (s := Iio b))
    have hr := tendsto_nhds_unique right (hc.continuousWithinAt (s := Ioi b))
    simpa using hr.trans hl.symm
  · intro h
    have at_b : wall p w b = below p w b := by
      rw [midpoint_value p w hw b, h]
      simp
    apply continuousAt_iff_continuous_left'_right'.2
    constructor
    · change Tendsto _ _ (𝓝 (wall p w b))
      rw [at_b]
      exact left
    · change Tendsto _ _ (𝓝 (wall p w b))
      rw [at_b]
      simpa [h] using right

theorem continuousOn_iff_onWall_zero (U : Set ℝ) (hU : IsOpen U) :
    ContinuousOn (wall p w) U ↔ ∀ s ∈ U, onWall p w s = 0 := by
  constructor
  · intro hc s hs
    exact (continuousAt_iff_onWall_zero p w hw s).1 ((hc s hs).continuousAt (hU.mem_nhds hs))
  · intro h s hs
    exact ((continuousAt_iff_onWall_zero p w hw s).2 (h s hs)).continuousWithinAt

theorem wall_eq_zero_of_onWall_zero (h : ∀ x, onWall p w x = 0) : wall p w b = 0 := by
  classical
  rw [midpoint_value p w hw b, h b, zero_div, add_zero]
  unfold below
  have grouped := ((selected_summable w hw (fun i => p i < b)).hasSum.tsum_fiberwise p).tsum_eq
  have zero_group (x : ℝ) :
      (∑' j : p ⁻¹' {x}, if p j < b then w j else 0) = 0 := by
    by_cases hx : x < b
    · calc
        _ = ∑' j : p ⁻¹' {x}, w j := tsum_congr fun j => by
          have hj : p j = x := j.property
          rw [hj, if_pos hx]
        _ = onWall p w x := fibre_sum p w x
        _ = 0 := h x
    · have hz (j : p ⁻¹' {x}) : (if p j < b then w j else 0) = 0 := by
        have hj : p j = x := j.property
        rw [hj, if_neg hx]
      simp_rw [hz]
      exact tsum_zero
  simp_rw [zero_group] at grouped
  simpa using grouped.symm

theorem exteriorTail_eq (pole : ℂ) :
    exteriorTail p w pole b = (pole - ∑' i, w i) + wall p w b := by
  have outer := Summable.of_norm_bounded hw (fun i => step_bound b (p i) (w i))
  have inner := Summable.of_norm_bounded hw (fun i => step_bound (p i) b (w i))
  have both : (∑' i, step b (w i) (p i)) + wall p w b = ∑' i, w i := by
    unfold wall
    rw [← outer.tsum_add inner]
    apply tsum_congr
    intro i
    rcases lt_trichotomy (p i) b with hi | hi | hi
    · simp [step, hi, not_lt_of_ge hi.le, ne_of_gt hi]
    · simp [step, hi, add_halves]
    · simp [step, hi, not_lt_of_ge hi.le, ne_of_gt hi]
  unfold exteriorTail
  rw [← both]
  abel

theorem exteriorTail_continuousAt_iff (pole : ℂ) :
    ContinuousAt (exteriorTail p w pole) b ↔ onWall p w b = 0 := by
  have eq_fun : exteriorTail p w pole = fun s => (pole - ∑' i, w i) + wall p w s :=
    funext fun s => exteriorTail_eq p w hw s pole
  rw [eq_fun, ← continuousAt_iff_onWall_zero p w hw b]
  constructor
  · intro h
    have hsub := h.sub (continuousAt_const (y := pole - ∑' i, w i))
    change ContinuousAt (fun s => ((pole - ∑' i, w i) + wall p w s) - (pole - ∑' i, w i)) b at hsub
    simpa using hsub
  · exact fun h => continuousAt_const.add h

theorem exteriorTail_eq_pole (pole : ℂ) (h : ∀ x, onWall p w x = 0) :
    exteriorTail p w pole b = pole := by
  have total := (hw.of_norm.hasSum.tsum_fiberwise p).tsum_eq
  simp_rw [fibre_sum p w, h] at total
  rw [exteriorTail_eq p w hw b pole, wall_eq_zero_of_onWall_zero p w hw b h]
  simpa using congrArg (fun z : ℂ => pole - z) total.symm

end PrimeTwoContinuity

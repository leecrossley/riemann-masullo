import PrimeTwoIncrements
import Mathlib.Analysis.Asymptotics.Defs
import Mathlib.MeasureTheory.Integral.Prod
import Mathlib.MeasureTheory.Integral.Bochner.ContinuousLinearMap

open MeasureTheory Set Filter
open scoped ComplexConjugate

namespace PrimeTwoWallEnergy

variable {α : Type*} [MeasurableSpace α]

def wallWindow (position : α → ℝ) (h s : ℝ) : Set α :=
  {x | s ∈ Ioo (position x - h / 2) (position x + h / 2)}

private theorem measurableSet_wallWindow {position : α → ℝ} (hp : Measurable position)
    (h s : ℝ) : MeasurableSet (wallWindow position h s) :=
  (measurableSet_lt (hp.sub measurable_const) measurable_const).inter
    (measurableSet_lt measurable_const (hp.add measurable_const))

private theorem common_interval (ω : ℝ → ℝ) {h : ℝ} (hh : 0 < h) (a b : ℝ) :
    ∫ s in Ioo (a - h / 2) (a + h / 2) ∩ Ioo (b - h / 2) (b + h / 2), ω s =
      h * PrimeTwoIncrements.overlap ω h a b := by
  rw [Ioo_inter_Ioo, max_sub_sub_right, min_add_add_right]
  dsimp only [PrimeTwoIncrements.overlap]
  split
  · next hle =>
    rw [intervalIntegral.integral_of_le hle, integral_Ioc_eq_integral_Ioo]
    field_simp
  · next hle =>
    rw [Ioo_eq_empty_of_le (le_of_not_ge hle), setIntegral_empty, mul_zero]

theorem wall_energy_identity (ν : Measure α) [IsFiniteMeasure ν]
    {position : α → ℝ} (hp : Measurable position) {phase : α → ℂ} (hm : Measurable phase)
    (hphase : ∀ x, ‖phase x‖ ≤ 1) {ω : ℝ → ℝ} (hω : Integrable ω)
    {h : ℝ} (hh : 0 < h) :
    Integrable (fun s => ω s * ‖∫ x in wallWindow position h s, phase x ∂ν‖ ^ 2) ∧
      ((∫ s, ω s * ‖∫ x in wallWindow position h s, phase x ∂ν‖ ^ 2 : ℝ) : ℂ) =
        h • (∫ x, ∫ y, PrimeTwoIncrements.overlap ω h (position x) (position y) •
          (phase x * conj (phase y)) ∂ν ∂ν) := by
  classical
  let window := wallWindow position h
  let both : Set (ℝ × (α × α)) := {v | v.2.1 ∈ window v.1 ∧ v.2.2 ∈ window v.1}
  let pair (t : α × α) := phase t.1 * conj (phase t.2)
  let kernel (s : ℝ) (t : α × α) : ℂ :=
    both.indicator (fun v => ω v.1 • pair v.2) (s, t)
  have hwindow (s : ℝ) : MeasurableSet (window s) :=
    measurableSet_wallWindow hp h s
  have slab {point : ℝ × (α × α) → α} (hq : Measurable point) :
      MeasurableSet {v | point v ∈ window v.1} :=
    (measurableSet_lt ((hp.comp hq).sub measurable_const) measurable_fst).inter
      (measurableSet_lt measurable_fst ((hp.comp hq).add measurable_const))
  have hboth : MeasurableSet both :=
    (slab (point := fun v => v.2.1) (measurable_fst.comp measurable_snd)).inter
      (slab (point := fun v => v.2.2) (measurable_snd.comp measurable_snd))
  have hpair : Integrable pair (ν.prod ν) := by
    apply (integrable_const (1 : ℝ)).mono'
    · exact ((hm.comp measurable_fst).mul
        (Complex.continuous_conj.measurable.comp (hm.comp measurable_snd))).aestronglyMeasurable
    · apply Eventually.of_forall
      intro t
      simpa [pair, norm_mul] using
        (mul_le_mul (hphase t.1) (hphase t.2) (norm_nonneg _) zero_le_one)
  have hkernel : Integrable (Function.uncurry kernel) (volume.prod (ν.prod ν)) :=
    (hω.smul_prod hpair).indicator hboth
  have atWall (s : ℝ) : ∫ t, kernel s t ∂ν.prod ν =
      ((ω s * ‖∫ x in window s, phase x ∂ν‖ ^ 2 : ℝ) : ℂ) := by
    change ∫ t, (window s ×ˢ window s).indicator (fun t => ω s • pair t) t ∂ν.prod ν = _
    rw [integral_indicator ((hwindow s).prod (hwindow s)), integral_smul]
    change ω s • (∫ t in window s ×ˢ window s, phase t.1 * conj (phase t.2) ∂ν.prod ν) = _
    rw [setIntegral_prod_mul (μ := ν) (ν := ν) phase (fun y => conj (phase y)) (window s) (window s),
      integral_conj, Complex.mul_conj']
    simp [Complex.real_smul]
  have atPair (t : α × α) : ∫ s, kernel s t =
      h • (PrimeTwoIncrements.overlap ω h (position t.1) (position t.2) • pair t) := by
    let interval (x : α) := Ioo (position x - h / 2) (position x + h / 2)
    change ∫ s, (interval t.1 ∩ interval t.2).indicator (fun s => ω s • pair t) s = _
    rw [integral_indicator (measurableSet_Ioo.inter measurableSet_Ioo), integral_smul_const,
      common_interval ω hh, mul_smul]
  have energyIntegrable := hkernel.integral_prod_left.congr (Eventually.of_forall atWall)
  refine ⟨?_, ?_⟩
  · exact energyIntegrable.re
  · rw [← integral_complex_ofReal]
    calc
      _ = ∫ s, ∫ t, kernel s t ∂ν.prod ν :=
        integral_congr_ae (Eventually.of_forall fun s => (atWall s).symm)
      _ = ∫ t, (∫ s, kernel s t) ∂ν.prod ν := integral_integral_swap hkernel
      _ = ∫ x, ∫ y, (∫ s, kernel s (x, y)) ∂ν ∂ν := integral_prod _ hkernel.integral_prod_right
      _ = _ := by simp_rw [atPair, integral_smul, pair]

theorem wall_energy_bound (ν : Measure α) [IsFiniteMeasure ν]
    {position : α → ℝ} (hp : Measurable position) {phase : α → ℂ} (hm : Measurable phase)
    (hphase : ∀ x, ‖phase x‖ ≤ 1) {ω : ℝ → ℝ} (hω : Integrable ω)
    {C : ℝ} (hC : ∀ s, ‖ω s‖ ≤ C) {h : ℝ} (hh : 0 < h) :
    |∫ s, ω s * ‖∫ x in wallWindow position h s, phase x ∂ν‖ ^ 2| ≤
      h * C * (ν.real univ) ^ 2 := by
  have hpair (x y : α) : ‖phase x * conj (phase y)‖ ≤ 1 := by
    simpa [norm_mul] using (mul_le_mul (hphase x) (hphase y) (norm_nonneg _) zero_le_one)
  have bound := PrimeTwoIncrements.integrated_overlap_bound ν position _ hpair hC hh
  rw [← (wall_energy_identity ν hp hm hphase hω hh).2] at bound
  simpa only [Complex.norm_real, Real.norm_eq_abs] using bound

theorem density_wall_energy_bound (μ : Measure α) {position : α → ℝ} (hp : Measurable position)
    {density : α → ℂ} (hm : Measurable density) (hd : Integrable density μ)
    {ω : ℝ → ℝ} (hω : Integrable ω) {C : ℝ} (hC : ∀ s, ‖ω s‖ ≤ C) {h : ℝ} (hh : 0 < h) :
    Integrable (fun s => ω s * ‖∫ x in wallWindow position h s, density x ∂μ‖ ^ 2) ∧
      |∫ s, ω s * ‖∫ x in wallWindow position h s, density x ∂μ‖ ^ 2| ≤
        h * C * (∫ x, ‖density x‖ ∂μ) ^ 2 := by
  let ν := μ.withDensity (fun x => ‖density x‖ₑ)
  let phase (x : α) := ‖density x‖⁻¹ • density x
  have mass : ν univ = ∫⁻ x, ‖density x‖ₑ ∂μ := by
    dsimp only [ν]
    rw [withDensity_apply _ MeasurableSet.univ, Measure.restrict_univ]
  let : IsFiniteMeasure ν := ⟨mass ▸ hasFiniteIntegral_iff_enorm.mp hd.2⟩
  have hmPhase : Measurable phase := hm.norm.inv.smul hm
  have normPhase (x : α) : ‖phase x‖ ≤ 1 := by
    by_cases hx : density x = 0
    · simp [phase, hx]
    · simp [phase, hx]
  have recover (x : α) : ‖density x‖ₑ.toReal • phase x = density x := by
    by_cases hx : density x = 0
    · simp [phase, hx]
    · simp [phase, hx]
  have wall (s : ℝ) : (∫ x in wallWindow position h s, phase x ∂ν) =
      ∫ x in wallWindow position h s, density x ∂μ := by
    dsimp only [ν]
    rw [setIntegral_withDensity_eq_setIntegral_toReal_smul hm.enorm
      (Eventually.of_forall fun _ => enorm_lt_top) phase (measurableSet_wallWindow hp h s)]
    simp_rw [recover]
  have total : ν.real univ = ∫ x, ‖density x‖ ∂μ := by
    rw [measureReal_def, mass, integral_norm_eq_lintegral_enorm hd.aestronglyMeasurable]
  exact ⟨by simpa only [wall] using (wall_energy_identity ν hp hmPhase normPhase hω hh).1,
    by simpa only [wall, total] using wall_energy_bound ν hp hmPhase normPhase hω hC hh⟩

theorem tail_energy_bound (μ : Measure α) {position : α → ℝ} (hp : Measurable position)
    {density : α → ℂ} (hm : Measurable density) (hd : Integrable density μ)
    {ω : ℝ → ℝ} (hω : Integrable ω) (hpos : ∀ s, 0 ≤ ω s) {C : ℝ} (hC : ∀ s, ‖ω s‖ ≤ C)
    {h : ℝ} (hh : 0 < h) {tail : ℝ → ℂ} (ht : Measurable tail) {B : ℝ}
    (herr : ∀ s, 0 < ω s → ‖tail s - ∫ x in wallWindow position h s, density x ∂μ‖ ≤ B * h) :
    Integrable (fun s => ω s * ‖tail s‖ ^ 2) ∧
      ∫ s, ω s * ‖tail s‖ ^ 2 ≤
        2 * (h * C * (∫ x, ‖density x‖ ∂μ) ^ 2) + 2 * (B * h) ^ 2 * ∫ s, ω s := by
  let full (s : ℝ) := ∫ x in wallWindow position h s, density x ∂μ
  obtain ⟨hi, hb⟩ := density_wall_energy_bound μ hp hm hd hω hC hh
  change Integrable (fun s => ω s * ‖full s‖ ^ 2) at hi
  have major (s : ℝ) : ω s * ‖tail s‖ ^ 2 ≤
      2 * (ω s * ‖full s‖ ^ 2) + (2 * (B * h) ^ 2) * ω s := by
    rcases (hpos s).eq_or_lt with hs | hs
    · simp [← hs]
    have hn : ‖tail s‖ ≤ ‖full s‖ + B * h := by
      calc
        ‖tail s‖ = ‖(tail s - full s) + full s‖ := by rw [sub_add_cancel]
        _ ≤ ‖tail s - full s‖ + ‖full s‖ := norm_add_le _ _
        _ ≤ ‖full s‖ + B * h := by linarith [herr s hs]
    have hsq : (‖full s‖ + B * h) ^ 2 ≤ 2 * (‖full s‖ ^ 2 + (B * h) ^ 2) := by
      nlinarith [sq_nonneg (‖full s‖ - B * h)]
    calc
      _ ≤ ω s * (‖full s‖ + B * h) ^ 2 :=
        mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (norm_nonneg _) hn 2) (hpos s)
      _ ≤ ω s * (2 * (‖full s‖ ^ 2 + (B * h) ^ 2)) := mul_le_mul_of_nonneg_left hsq (hpos s)
      _ = _ := by ring
  have hmajor := (hi.const_mul 2).add (hω.const_mul (2 * (B * h) ^ 2))
  have htail : Integrable (fun s => ω s * ‖tail s‖ ^ 2) := by
    apply hmajor.mono' (hω.aestronglyMeasurable.mul (ht.norm.pow_const 2).aestronglyMeasurable)
    exact Eventually.of_forall fun s => by
      simpa only [Pi.mul_apply, Pi.add_apply, Real.norm_eq_abs,
        abs_of_nonneg (mul_nonneg (hpos s) (sq_nonneg _))] using major s
  have hle := integral_mono htail hmajor major
  simp only [Pi.add_apply] at hle
  rw [integral_add (hi.const_mul 2) (hω.const_mul _), integral_const_mul, integral_const_mul] at hle
  exact ⟨htail, hle.trans (add_le_add
    (mul_le_mul_of_nonneg_left ((le_abs_self _).trans hb) (by norm_num : (0 : ℝ) ≤ 2)) le_rfl)⟩

theorem tail_energy_isBigO {ι : Type*} (l : Filter ι) (μ : Measure α)
    {position : α → ℝ} (hp : Measurable position)
    {ω : ℝ → ℝ} (hω : Integrable ω) (hpos : ∀ s, 0 ≤ ω s) {C : ℝ} (hC : ∀ s, ‖ω s‖ ≤ C)
    {h : ι → ℝ} (hh : ∀ᶠ i in l, h i ∈ Ioc 0 1) {density : ι → α → ℂ} {M : ℝ}
    (hd : ∀ᶠ i in l, Measurable (density i) ∧ Integrable (density i) μ ∧
      (∫ x, ‖density i x‖ ∂μ) ≤ M)
    {tail : ι → ℝ → ℂ} {B : ℝ}
    (ht : ∀ᶠ i in l, Measurable (tail i) ∧ ∀ s, 0 < ω s →
      ‖tail i s - ∫ x in wallWindow position (h i) s, density i x ∂μ‖ ≤ B * h i) :
    (fun i => ∫ s, ω s * ‖tail i s‖ ^ 2) =O[l] h := by
  apply Asymptotics.isBigO_iff.mpr
  refine ⟨2 * C * M ^ 2 + 2 * B ^ 2 * (∫ s, ω s), ?_⟩
  filter_upwards [hh, hd, ht] with i hs hden htail
  obtain ⟨hm, hi, hM⟩ := hden
  obtain ⟨hmt, herr⟩ := htail
  have hmass := pow_le_pow_left₀ (integral_nonneg fun x => norm_nonneg (density i x)) hM 2
  have hwidth : h i ^ 2 ≤ h i := by
    simpa only [pow_two, mul_one] using mul_le_mul_of_nonneg_left hs.2 hs.1.le
  have hC₀ := (norm_nonneg (ω 0)).trans (hC 0)
  have hweight : 0 ≤ ∫ s, ω s := integral_nonneg (f := ω) hpos
  have henergy : 0 ≤ ∫ s, ω s * ‖tail i s‖ ^ 2 :=
    integral_nonneg fun s => mul_nonneg (hpos s) (sq_nonneg _)
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg henergy, abs_of_pos hs.1]
  calc
    _ ≤ 2 * (h i * C * (∫ x, ‖density i x‖ ∂μ) ^ 2) + 2 * (B * h i) ^ 2 * ∫ s, ω s :=
      (tail_energy_bound μ hp hm hi hω hpos hC hs.1 hmt herr).2
    _ = (2 * (h i * C)) * (∫ x, ‖density i x‖ ∂μ) ^ 2 + (2 * B ^ 2 * ∫ s, ω s) * h i ^ 2 := by ring
    _ ≤ (2 * (h i * C)) * M ^ 2 + (2 * B ^ 2 * ∫ s, ω s) * h i :=
      add_le_add (mul_le_mul_of_nonneg_left hmass (mul_nonneg zero_le_two (mul_nonneg hs.1.le hC₀)))
        (mul_le_mul_of_nonneg_left hwidth (mul_nonneg (mul_nonneg zero_le_two (sq_nonneg _)) hweight))
    _ = _ := by ring

end PrimeTwoWallEnergy

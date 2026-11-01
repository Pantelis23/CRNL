/-
§112 — the merge threshold is rail geometry.

A one-species bistable element with rails `r₁ < r₂ < r₃` (low rail, saddle, high rail). Merging `k`
tanks of which `j` have failed to the low rail conserves counts, so the merged concentration is the
count-weighted mean `x(k, j) = ((k − j) r₃ + j r₁) / k`. FINDINGS §112 states, and its table rests on:

  * the merged tank sits below the saddle exactly when `j > f k`, with `f = (r₃ − r₂)/(r₃ − r₁)`;
  * hence the least number of failures that flips the merge is `m(k) = ⌊f k⌋ + 1`;
  * unanimity (`m = k`) holds exactly when `f ≥ 1 − 1/k`;
  * one failure flips it (`m = 1`) exactly when `f < 1/k`.

All four are proved here, with no assumption beyond `r₁ < r₂ < r₃` and `k ≥ 1`.
-/
import Mathlib

namespace CRNL.MergeLaw

variable {r₁ r₂ r₃ : ℝ}

/-- The merged concentration of `k` tanks, `j` of them at the low rail. -/
noncomputable def xMerged (r₁ r₃ : ℝ) (k j : ℕ) : ℝ := ((k - j) * r₃ + j * r₁) / k

/-- The rail-geometry fraction `f = (r₃ − r₂)/(r₃ − r₁)`. -/
noncomputable def frac (r₁ r₂ r₃ : ℝ) : ℝ := (r₃ - r₂) / (r₃ - r₁)

/-- The least number of low tanks that drags the merge below the saddle. -/
noncomputable def threshold (r₁ r₂ r₃ : ℝ) (k : ℕ) : ℕ := ⌊frac r₁ r₂ r₃ * k⌋₊ + 1

theorem frac_pos (h₁₂ : r₁ < r₂) (h₂₃ : r₂ < r₃) : 0 < frac r₁ r₂ r₃ :=
  div_pos (by linarith) (by linarith)

theorem frac_lt_one (h₁₂ : r₁ < r₂) (h₂₃ : r₂ < r₃) : frac r₁ r₂ r₃ < 1 := by
  rw [frac, div_lt_one (by linarith)]
  linarith

/-- The merged tank is below the saddle exactly when more than `f k` tanks have failed. -/
theorem below_saddle_iff (h₁₂ : r₁ < r₂) (h₂₃ : r₂ < r₃) {k : ℕ} (hk : 0 < k) (j : ℕ) :
    xMerged r₁ r₃ k j < r₂ ↔ frac r₁ r₂ r₃ * k < j := by
  have hk' : (0 : ℝ) < k := by exact_mod_cast hk
  have hd : (0 : ℝ) < r₃ - r₁ := by linarith
  rw [xMerged, div_lt_iff₀ hk', frac, div_mul_eq_mul_div, div_lt_iff₀ hd]
  constructor <;> intro h <;> nlinarith

/-- §112's merge law: the threshold `⌊f k⌋ + 1` is exactly the least flipping count. -/
theorem below_saddle_iff_threshold_le (h₁₂ : r₁ < r₂) (h₂₃ : r₂ < r₃) {k : ℕ} (hk : 0 < k)
    (j : ℕ) : xMerged r₁ r₃ k j < r₂ ↔ threshold r₁ r₂ r₃ k ≤ j := by
  have hf : 0 ≤ frac r₁ r₂ r₃ * k := mul_nonneg (frac_pos h₁₂ h₂₃).le (Nat.cast_nonneg k)
  rw [below_saddle_iff h₁₂ h₂₃ hk, threshold, Nat.add_one_le_iff, Nat.floor_lt hf]

/-- Unanimity regime: every tank must fail before the merge flips iff `f ≥ 1 − 1/k`. -/
theorem threshold_eq_k_iff (h₁₂ : r₁ < r₂) (h₂₃ : r₂ < r₃) {k : ℕ} (hk : 0 < k) :
    threshold r₁ r₂ r₃ k = k ↔ 1 - 1 / (k : ℝ) ≤ frac r₁ r₂ r₃ := by
  have hk' : (0 : ℝ) < k := by exact_mod_cast hk
  have hf0 : 0 ≤ frac r₁ r₂ r₃ * k := mul_nonneg (frac_pos h₁₂ h₂₃).le hk'.le
  have hflt : frac r₁ r₂ r₃ * k < k := by
    have := frac_lt_one h₁₂ h₂₃
    nlinarith
  rw [threshold]
  constructor
  · intro h
    have hfl : ⌊frac r₁ r₂ r₃ * k⌋₊ = k - 1 := by omega
    have : ((k - 1 : ℕ) : ℝ) ≤ frac r₁ r₂ r₃ * k := by
      rw [← hfl]; exact Nat.floor_le hf0
    rw [Nat.cast_sub (by omega), Nat.cast_one] at this
    rw [sub_le_comm, le_div_iff₀ hk']
    nlinarith
  · intro h
    have hge : ((k - 1 : ℕ) : ℝ) ≤ frac r₁ r₂ r₃ * k := by
      rw [Nat.cast_sub (by omega), Nat.cast_one]
      have := mul_le_mul_of_nonneg_right h hk'.le
      rwa [sub_mul, one_div, inv_mul_cancel₀ hk'.ne', one_mul] at this
    have h1 : k - 1 ≤ ⌊frac r₁ r₂ r₃ * k⌋₊ := Nat.le_floor hge
    have h2 : ⌊frac r₁ r₂ r₃ * k⌋₊ < k := (Nat.floor_lt hf0).mpr hflt
    omega

/-- One-failure regime: a single low tank flips the merge iff `f < 1/k`. -/
theorem threshold_eq_one_iff (h₁₂ : r₁ < r₂) (h₂₃ : r₂ < r₃) {k : ℕ} (hk : 0 < k) :
    threshold r₁ r₂ r₃ k = 1 ↔ frac r₁ r₂ r₃ < 1 / (k : ℝ) := by
  have hk' : (0 : ℝ) < k := by exact_mod_cast hk
  have hf0 : 0 ≤ frac r₁ r₂ r₃ * k := mul_nonneg (frac_pos h₁₂ h₂₃).le hk'.le
  rw [threshold, Nat.add_eq_right, Nat.floor_eq_zero, lt_div_iff₀ hk']

end CRNL.MergeLaw

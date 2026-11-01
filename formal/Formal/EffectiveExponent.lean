/-
§84 / rule 21 — the effective exponent is the closed form's own correction term.

FINDINGS §84 derived the escape action near the pitchfork in closed form,
  `A(γ) = (1 + 2γ)(1 − 2γ)² / (2(1 − γ)(1 + γ)²)`,
and §85.1 then observed that with `ε = ½ − γ` this is `(32/9) ε² (1 − (5/3)ε + O(ε²))`, so the local
exponent `ν_eff(ε) = d ln A / d ln ε` is `2 − (5/3)ε + O(ε²)` — the drift every window fit had read
as a disagreement. The record states it with `O(·)` and checks it numerically. Here it is exact:

  * `A(½ − ε) = 4ε²(1 − ε) / ((½ + ε)(3/2 − ε)²)`;
  * `ν_eff(ε) = 2 − ε/(1 − ε) − ε/(½ + ε) + 2ε/(3/2 − ε)` on `0 < ε < 1`, as a derivative;
  * `|ν_eff(ε) − (2 − (5/3)ε)| ≤ 7 ε²` on `0 < ε ≤ ¼` — the `O(ε²)` with its constant;
  * `ν_eff(1/10) = 235/126 = 1.86507…`, the value §85.1 printed as 1.8651.
-/
import Mathlib

namespace CRNL.EffectiveExponent

open Real

/-- §84's closed-form escape action. -/
noncomputable def A (γ : ℝ) : ℝ := (1 + 2 * γ) * (1 - 2 * γ) ^ 2 / (2 * (1 - γ) * (1 + γ) ^ 2)

/-- The same action in `ε = ½ − γ`. -/
noncomputable def Aε (ε : ℝ) : ℝ := 4 * ε ^ 2 * (1 - ε) / ((1 / 2 + ε) * (3 / 2 - ε) ^ 2)

theorem A_half_sub (ε : ℝ) :
    A (1 / 2 - ε) = Aε ε := by
  unfold A Aε
  have e₁ : 2 * (1 - (1 / 2 - ε)) = 2 * (1 / 2 + ε) := by ring
  have e₂ : (1 + (1 / 2 - ε)) = 3 / 2 - ε := by ring
  rw [e₁, e₂]
  field_simp
  ring

/-- The local exponent `ν_eff = ε · d ln A / dε`. -/
noncomputable def ν (ε : ℝ) : ℝ := ε * deriv (fun e => log (Aε e)) ε

/-- `ln A` as a sum of logs, valid where every factor is positive. -/
noncomputable def L (e : ℝ) : ℝ :=
  log 4 + 2 * log e + log (1 - e) - log (1 / 2 + e) - 2 * log (3 / 2 - e)

theorem log_Aε_eq_L {e : ℝ} (h0 : 0 < e) (h1 : e < 1) : log (Aε e) = L e := by
  have a : 0 < 1 - e := by linarith
  have b : 0 < 1 / 2 + e := by linarith
  have c : 0 < 3 / 2 - e := by linarith
  unfold Aε L
  rw [log_div (by positivity) (by positivity), log_mul (by positivity) a.ne',
    log_mul (by positivity) (by positivity), log_mul b.ne' (by positivity), log_pow, log_pow]
  push_cast
  ring

theorem hasDerivAt_L {e : ℝ} (h0 : 0 < e) (h1 : e < 1) :
    HasDerivAt L (2 / e - 1 / (1 - e) - 1 / (1 / 2 + e) + 2 / (3 / 2 - e)) e := by
  have a : (1 - e) ≠ 0 := by linarith
  have b : (1 / 2 + e) ≠ 0 := by linarith
  have c : (3 / 2 - e) ≠ 0 := by linarith
  have d1 := ((hasDerivAt_id e).log h0.ne').const_mul 2
  have d2 := (((hasDerivAt_id e).const_sub 1).log a)
  have d3 := (((hasDerivAt_id e).const_add (1 / 2)).log b)
  have d4 := (((hasDerivAt_id e).const_sub (3 / 2)).log c).const_mul 2
  have := ((((hasDerivAt_const e (log 4)).add d1).add d2).sub d3).sub d4
  convert this using 1
  · ext y; simp [L]
  · simp; field_simp; ring

/-- §85.1's exponent, exactly. -/
theorem ν_eq {ε : ℝ} (h0 : 0 < ε) (h1 : ε < 1) :
    ν ε = 2 - ε / (1 - ε) - ε / (1 / 2 + ε) + 2 * ε / (3 / 2 - ε) := by
  have hev : (fun e => log (Aε e)) =ᶠ[nhds ε] L := by
    filter_upwards [Ioo_mem_nhds h0 h1] with e he
    exact log_Aε_eq_L he.1 he.2
  have hd : deriv (fun e => log (Aε e)) ε = _ := hev.deriv_eq.trans (hasDerivAt_L h0 h1).deriv
  rw [ν, hd]
  have a : (1 - ε) ≠ 0 := by linarith
  have b : (1 / 2 + ε) ≠ 0 := by linarith
  have c : (3 / 2 - ε) ≠ 0 := by linarith
  field_simp

/-- The remainder after the leading correction, in closed form. -/
theorem ν_sub_linear {ε : ℝ} (h0 : 0 < ε) (h1 : ε < 1) :
    ν ε - (2 - 5 / 3 * ε) =
      ε ^ 2 * (-1 / (1 - ε) + 4 / (1 + 2 * ε) + 8 / (9 - 6 * ε)) := by
  rw [ν_eq h0 h1]
  have a : (1 - ε) ≠ 0 := by linarith
  have b : (1 / 2 + ε) ≠ 0 := by linarith
  have b' : (1 + 2 * ε) ≠ 0 := by linarith
  have c : (3 / 2 - ε) ≠ 0 := by linarith
  have c' : (9 - 6 * ε) ≠ 0 := by linarith
  grind

/-- `ν_eff = 2 − (5/3)ε + O(ε²)`, with the constant: the error is at most `7 ε²` for `ε ≤ ¼`. -/
theorem ν_close {ε : ℝ} (h0 : 0 < ε) (h1 : ε ≤ 1 / 4) :
    |ν ε - (2 - 5 / 3 * ε)| ≤ 7 * ε ^ 2 := by
  rw [ν_sub_linear h0 (by linarith), abs_mul, abs_of_nonneg (sq_nonneg ε), mul_comm]
  refine mul_le_mul_of_nonneg_right ?_ (sq_nonneg ε)
  have a : 0 < 1 - ε := by linarith
  have b : 0 < 1 + 2 * ε := by linarith
  have c : 0 < 9 - 6 * ε := by linarith
  have t1 : 1 / (1 - ε) ≤ 4 / 3 := by rw [div_le_div_iff₀ a (by norm_num)]; linarith
  have t2 : 4 / (1 + 2 * ε) ≤ 4 := by rw [div_le_iff₀ b]; linarith
  have t3 : 8 / (9 - 6 * ε) ≤ 16 / 15 := by rw [div_le_div_iff₀ c (by norm_num)]; linarith
  have p1 : 0 ≤ 1 / (1 - ε) := by positivity
  have p2 : 0 ≤ 4 / (1 + 2 * ε) := by positivity
  have p3 : 0 ≤ 8 / (9 - 6 * ε) := by positivity
  rw [neg_div, abs_le]
  constructor <;> linarith

/-- The number §85.1 printed: `ν_eff(0.1) = 1.8651`. Exactly `235/126`. -/
theorem ν_tenth : ν (1 / 10) = 235 / 126 := by
  rw [ν_eq (by norm_num) (by norm_num)]
  norm_num

end CRNL.EffectiveExponent

/-
The bistability edge is well defined: the bistable inputs form an interval.

`weak_upstream.edge()` (and §105's `bistability_edge`) scans the upstream level and returns the least
one at which stage 2's rate cubic
    F_h(x) = h k₁ₐ x² + k₂ᵦ − k₁ᵣ x³ − k₂ᵣ x          (h = the Hill factor of the input)
has three positive roots. That is only "the edge" if those inputs form an interval. Proved here, for
any positive rate constants:

  * every real root of `F_h` is positive (for `x ≤ 0` every term is `≥ 0` and `k₂ᵦ > 0`), so
    "three positive roots" is "three real roots";
  * the cubic's discriminant, as a function of `h`, is
        Δ(h) = −4 k₁ₐ³ k₂ᵦ h³ + k₁ₐ² k₂ᵣ² h² + 18 k₁ᵣ k₁ₐ k₂ᵣ k₂ᵦ h − (4 k₁ᵣ k₂ᵣ³ + 27 k₁ᵣ² k₂ᵦ²)
    (`disc_eq`, checked against the general formula);
  * `Δ(h)/h` is concave on `h > 0`, so `{h > 0 | Δ h > 0}` is convex — an interval.

Not proved here: the classical fact that a real cubic has three distinct real roots iff its
discriminant is positive. The interval statement is about the discriminant.
-/
import Mathlib

namespace CRNL.EdgeInterval

open Set

variable {k1a k1r k2b k2r : ℝ}

/-- Stage 2's rate cubic at Hill level `h`. -/
def F (k1a k1r k2b k2r h x : ℝ) : ℝ := h * k1a * x ^ 2 + k2b - k1r * x ^ 3 - k2r * x

/-- Every real root is positive: `F` is strictly positive on `x ≤ 0`. -/
theorem F_pos_of_nonpos (ha : 0 < k1a) (hr : 0 < k1r) (hb : 0 < k2b) (hrr : 0 < k2r) {h x : ℝ}
    (hh : 0 ≤ h) (hx : x ≤ 0) : 0 < F k1a k1r k2b k2r h x := by
  unfold F
  have t1 : 0 ≤ h * k1a * x ^ 2 := by positivity
  have t2 : 0 ≤ -(k1r * x ^ 3) := by
    have : x ^ 3 ≤ 0 := by
      have := pow_le_pow_left₀ (neg_nonneg.mpr hx) le_rfl 3
      nlinarith [sq_nonneg x, mul_nonneg (neg_nonneg.mpr hx) (sq_nonneg x)]
    nlinarith
  have t3 : 0 ≤ -(k2r * x) := by nlinarith
  linarith

/-- The general cubic discriminant `18abcd − 4b³d + b²c² − 4ac³ − 27a²d²`. -/
def cubicDisc (a b c d : ℝ) : ℝ :=
  18 * a * b * c * d - 4 * b ^ 3 * d + b ^ 2 * c ^ 2 - 4 * a * c ^ 3 - 27 * a ^ 2 * d ^ 2

/-- The discriminant of `F_h` as a polynomial in `h`. -/
def Δ (k1a k1r k2b k2r h : ℝ) : ℝ :=
  -4 * k1a ^ 3 * k2b * h ^ 3 + k1a ^ 2 * k2r ^ 2 * h ^ 2 + 18 * k1r * k1a * k2r * k2b * h -
    (4 * k1r * k2r ^ 3 + 27 * k1r ^ 2 * k2b ^ 2)

/-- `F_h = a x³ + b x² + c x + d` with `a = −k₁ᵣ, b = h k₁ₐ, c = −k₂ᵣ, d = k₂ᵦ`; its discriminant is `Δ h`. -/
theorem disc_eq (h : ℝ) : cubicDisc (-k1r) (h * k1a) (-k2r) k2b = Δ k1a k1r k2b k2r h := by
  unfold cubicDisc Δ; ring

/-- `Δ(h)/h` for `h > 0`. -/
noncomputable def q (k1a k1r k2b k2r h : ℝ) : ℝ :=
  -4 * k1a ^ 3 * k2b * h ^ 2 + k1a ^ 2 * k2r ^ 2 * h + 18 * k1r * k1a * k2r * k2b -
    (4 * k1r * k2r ^ 3 + 27 * k1r ^ 2 * k2b ^ 2) * h⁻¹

theorem Δ_eq_mul_q {h : ℝ} (hh : 0 < h) : Δ k1a k1r k2b k2r h = h * q k1a k1r k2b k2r h := by
  unfold Δ q; field_simp

theorem q_concave (ha : 0 < k1a) (hr : 0 < k1r) (hb : 0 < k2b) (hrr : 0 < k2r) :
    ConcaveOn ℝ (Ioi 0) (q k1a k1r k2b k2r) := by
  have c2 : ConcaveOn ℝ (Ioi (0:ℝ)) fun h => -(4 * k1a ^ 3 * k2b) * h ^ 2 :=
    ((even_two.convexOn_pow).subset (subset_univ _) (convex_Ioi 0)).neg.smul
      (by positivity : (0:ℝ) ≤ 4 * k1a ^ 3 * k2b) |>.congr fun h _ => by simp
  have c1 : ConcaveOn ℝ (Ioi (0:ℝ)) fun h => k1a ^ 2 * k2r ^ 2 * h + 18 * k1r * k1a * k2r * k2b :=
    ((concaveOn_id (convex_Ioi (0:ℝ))).smul (by positivity : (0:ℝ) ≤ k1a ^ 2 * k2r ^ 2)).add_const
      (18 * k1r * k1a * k2r * k2b) |>.congr fun h _ => by simp [smul_eq_mul]
  have c0 : ConcaveOn ℝ (Ioi (0:ℝ)) fun h => -(4 * k1r * k2r ^ 3 + 27 * k1r ^ 2 * k2b ^ 2) * h⁻¹ :=
    ((convexOn_zpow (𝕜 := ℝ) (-1)).neg.smul
      (by positivity : (0:ℝ) ≤ 4 * k1r * k2r ^ 3 + 27 * k1r ^ 2 * k2b ^ 2))
      |>.congr fun h _ => by simp; ring
  refine ((c2.add c1).add c0).congr fun h _ => ?_
  simp only [q, Pi.add_apply]; ring

/-- The edge theorem: the Hill levels at which the discriminant is positive form an interval. -/
theorem bistable_levels_convex (ha : 0 < k1a) (hr : 0 < k1r) (hb : 0 < k2b) (hrr : 0 < k2r) :
    Convex ℝ {h ∈ Ioi (0:ℝ) | 0 < Δ k1a k1r k2b k2r h} := by
  have hc := (q_concave ha hr hb hrr).convex_gt 0
  convert hc using 1
  ext h
  simp only [mem_ofPred_eq, mem_Ioi]
  constructor
  · rintro ⟨hh, hΔ⟩
    exact ⟨hh, by rw [Δ_eq_mul_q hh] at hΔ; exact pos_of_mul_pos_right hΔ hh.le⟩
  · rintro ⟨hh, hq⟩
    exact ⟨hh, by rw [Δ_eq_mul_q hh]; exact mul_pos hh hq⟩

end CRNL.EdgeInterval

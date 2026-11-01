/-
The Schlögl element's construction from its rails (`experiments/cascade_schlogl.schlogl_consts`).

Every element in CRNL is built by choosing three rails `r₁ < r₂ < r₃` and setting the four rate
constants to the Vieta coefficients, `k₁ₐ = k e₁, k₁ᵣ = k, k₂ᵦ = k e₃, k₂ᵣ = k e₂`. The code assumes,
and every section uses, that this makes the mass-action drift vanish exactly at the rails, with the
outer rails stable and the middle one a saddle. Proved here.
-/
import Mathlib

namespace CRNL.Schlogl

/-- The rate constants `schlogl_consts(r1, r2, r3, k1r)`, in the code's order. -/
structure Consts where
  k1a : ℝ
  k1r : ℝ
  k2b : ℝ
  k2r : ℝ

/-- `schlogl_consts`: `(k1r·e₁, k1r, k1r·e₃, k1r·e₂)`. -/
def consts (r₁ r₂ r₃ k : ℝ) : Consts :=
  ⟨k * (r₁ + r₂ + r₃), k, k * (r₁ * r₂ * r₃), k * (r₁ * r₂ + r₁ * r₃ + r₂ * r₃)⟩

/-- The mass-action drift of one undriven element, `λ(x) − μ(x)` in concentration units. -/
def drift (c : Consts) (x : ℝ) : ℝ := c.k1a * x ^ 2 + c.k2b - c.k1r * x ^ 3 - c.k2r * x

/-- Its derivative. -/
def drift' (c : Consts) (x : ℝ) : ℝ := 2 * c.k1a * x - 3 * c.k1r * x ^ 2 - c.k2r

theorem hasDerivAt_drift (c : Consts) (x : ℝ) : HasDerivAt (drift c) (drift' c x) x := by
  have := (((hasDerivAt_pow 2 x).const_mul c.k1a).add_const c.k2b).sub
    ((hasDerivAt_pow 3 x).const_mul c.k1r) |>.sub ((hasDerivAt_id x).const_mul c.k2r)
  convert this using 1
  · ext y; simp [drift]
  · simp [drift']; ring

/-- The construction factorises the drift: it is `−k (x − r₁)(x − r₂)(x − r₃)`. -/
theorem drift_eq (r₁ r₂ r₃ k x : ℝ) :
    drift (consts r₁ r₂ r₃ k) x = -k * (x - r₁) * (x - r₂) * (x - r₃) := by
  simp [drift, consts]; ring

theorem drift_rails (r₁ r₂ r₃ k : ℝ) :
    drift (consts r₁ r₂ r₃ k) r₁ = 0 ∧ drift (consts r₁ r₂ r₃ k) r₂ = 0 ∧
      drift (consts r₁ r₂ r₃ k) r₃ = 0 := by
  simp [drift_eq]

/-- The slope at each rail, in closed form. -/
theorem drift'_rails (r₁ r₂ r₃ k : ℝ) :
    drift' (consts r₁ r₂ r₃ k) r₁ = -k * (r₁ - r₂) * (r₁ - r₃) ∧
      drift' (consts r₁ r₂ r₃ k) r₂ = k * (r₂ - r₁) * (r₃ - r₂) ∧
      drift' (consts r₁ r₂ r₃ k) r₃ = -k * (r₃ - r₁) * (r₃ - r₂) := by
  simp [drift', consts]; refine ⟨by ring, by ring, by ring⟩

/-- Ordered rails and a positive rate give two stable rails and an unstable saddle between them. -/
theorem stability {r₁ r₂ r₃ k : ℝ} (hk : 0 < k) (h₁₂ : r₁ < r₂) (h₂₃ : r₂ < r₃) :
    drift' (consts r₁ r₂ r₃ k) r₁ < 0 ∧ 0 < drift' (consts r₁ r₂ r₃ k) r₂ ∧
      drift' (consts r₁ r₂ r₃ k) r₃ < 0 := by
  obtain ⟨h1, h2, h3⟩ := drift'_rails r₁ r₂ r₃ k
  rw [h1, h2, h3]
  have a : 0 < r₂ - r₁ := by linarith
  have b : 0 < r₃ - r₂ := by linarith
  have c : 0 < r₃ - r₁ := by linarith
  refine ⟨?_, ?_, ?_⟩
  · have : -k * (r₁ - r₂) * (r₁ - r₃) = -(k * (r₂ - r₁) * (r₃ - r₁)) := by ring
    rw [this]; exact neg_neg_of_pos (by positivity)
  · positivity
  · have : -k * (r₃ - r₁) * (r₃ - r₂) = -(k * (r₃ - r₁) * (r₃ - r₂)) := by ring
    rw [this]; exact neg_neg_of_pos (by positivity)

/-- All four constructed rate constants are positive whenever the rails are positive — the check
`Elem.__init__` makes numerically ("rails give a non-positive rate constant"). -/
theorem consts_pos {r₁ r₂ r₃ k : ℝ} (hk : 0 < k) (h₁ : 0 < r₁) (h₂ : 0 < r₂) (h₃ : 0 < r₃) :
    let c := consts r₁ r₂ r₃ k
    0 < c.k1a ∧ 0 < c.k1r ∧ 0 < c.k2b ∧ 0 < c.k2r := by
  simp only [consts]
  refine ⟨by positivity, hk, by positivity, by positivity⟩

end CRNL.Schlogl

/-
§122 — the committor of a birth-death chain, the closed form behind φ.

§122 found that φ = w_d/(k₁t) — the stage-1 factor of the near-miss floor — equals the probability
that stage 1, one step below its saddle, reaches its high rail before its low rail. That probability
is computed (`asymptotic_floor.committor`) by the closed form below, on states shifted so that the
low rail is `0` and the high rail is `N`:

    ρ j = ∏_{i=1}^{j} μ i / λ i,      h n = (∑_{j<n} ρ j) / (∑_{j<N} ρ j).

Proved here: `h 0 = 0`, `h N = 1`, and `h` solves the backward equation
`λ n (h (n+1) − h n) + μ n (h (n−1) − h n) = 0` at every interior state, given only positive birth
rates. (The identification of φ with `h` is a measurement, §122; this is only the formula.)
-/
import Mathlib

namespace CRNL.Committor

open Finset

variable (lam mu : ℕ → ℝ)

/-- `ρ j = ∏_{i=1}^{j} μ i / λ i`. -/
noncomputable def ρ (j : ℕ) : ℝ := ∏ i ∈ range j, mu (i + 1) / lam (i + 1)

/-- Partial sums `∑_{j<n} ρ j`. -/
noncomputable def S (n : ℕ) : ℝ := ∑ j ∈ range n, ρ lam mu j

/-- The committor to `N` before `0`. -/
noncomputable def h (N n : ℕ) : ℝ := S lam mu n / S lam mu N

theorem ρ_succ (j : ℕ) : ρ lam mu (j + 1) = ρ lam mu j * (mu (j + 1) / lam (j + 1)) := by
  unfold ρ; rw [prod_range_succ]

theorem ρ_pos (hlam : ∀ n, 0 < lam n) (hmu : ∀ n, 0 < mu n) (j : ℕ) : 0 < ρ lam mu j :=
  prod_pos fun _ _ => div_pos (hmu _) (hlam _)

theorem h_zero (N : ℕ) : h lam mu N 0 = 0 := by simp [h, S]

theorem h_top (hlam : ∀ n, 0 < lam n) (hmu : ∀ n, 0 < mu n) {N : ℕ} (hN : 0 < N) :
    h lam mu N N = 1 := by
  have : 0 < S lam mu N := by
    obtain ⟨k, rfl⟩ : ∃ k, N = k + 1 := ⟨N - 1, by omega⟩
    exact sum_pos (fun j _ => ρ_pos lam mu hlam hmu j) ⟨0, by simp⟩
  exact div_self this.ne'

/-- The backward equation at every interior state `0 < n`. -/
theorem backward (hlam : ∀ n, 0 < lam n) (N n : ℕ) (hn : 0 < n) :
    lam n * (h lam mu N (n + 1) - h lam mu N n) + mu n * (h lam mu N (n - 1) - h lam mu N n) = 0 := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
  simp only [h, Nat.add_sub_cancel, S, sum_range_succ]
  rw [ρ_succ]
  have hl : lam (k + 1) ≠ 0 := (hlam _).ne'
  rcases eq_or_ne (∑ j ∈ range N, ρ lam mu j) 0 with h0 | h0
  · simp [h0]
  · field_simp
    ring

end CRNL.Committor

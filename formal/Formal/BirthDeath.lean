/-
The birth-death product formula behind `Elem.stage1_stationary` and `pi_low`.

A chain on `{0, …, N}` (a reflecting box: no birth out of `N`, no death out of `0`) with birth rates
`λ n` (n → n+1) and death rates `μ n` (n → n−1). The code takes as its exact stationary law

    π n ∝ ∏_{i < n} λ i / μ (i+1)

computed in log space. Proved here: positive death rates make this unnormalised `π` satisfy detailed
balance on every edge and global balance `π Q = 0` at every state of the box, boundaries included.
-/
import Mathlib

namespace CRNL.BirthDeath

open Finset

variable (lam mu : ℕ → ℝ)

/-- The unnormalised product-formula weight. -/
noncomputable def w (n : ℕ) : ℝ := ∏ i ∈ range n, lam i / mu (i + 1)

theorem w_succ (n : ℕ) : w lam mu (n + 1) = w lam mu n * (lam n / mu (n + 1)) := by
  unfold w; rw [prod_range_succ]

/-- Detailed balance across every edge `n ↔ n+1`. -/
theorem detailed_balance (hmu : ∀ n, mu (n + 1) ≠ 0) (n : ℕ) :
    w lam mu n * lam n = w lam mu (n + 1) * mu (n + 1) := by
  rw [w_succ]; field_simp [hmu n]

/-- Global balance `(π Q)ₙ = 0` in the box `{0,…,N}`, stated as inflow = outflow with no reference to
detailed balance: births stop at `N`, deaths stop at `0`. -/
theorem global_balance (hmu : ∀ n, mu (n + 1) ≠ 0) (N n : ℕ) :
    (if 0 < n then w lam mu (n - 1) * lam (n - 1) else 0) +
        (if n < N then w lam mu (n + 1) * mu (n + 1) else 0) =
      w lam mu n * ((if n < N then lam n else 0) + (if 0 < n then mu n else 0)) := by
  have db := detailed_balance lam mu hmu
  rcases Nat.eq_zero_or_pos n with rfl | hpos
  · simp only [lt_irrefl, ite_false, zero_add, add_zero]
    split_ifs <;> simp [db 0]
  · obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
    simp only [hpos, ite_true, Nat.add_sub_cancel]
    split_ifs
    · rw [db k, ← db (k + 1)]; ring
    · rw [db k]; ring

end CRNL.BirthDeath

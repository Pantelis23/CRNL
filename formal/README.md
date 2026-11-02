# formal/ — machine-checked proofs of CRNL's exact claims

Lean 4 + Mathlib. Only claims that are *exact* belong here: algebra, closed forms, combinatorics.
Measured results (everything from the CME, every fitted law) do not, and nothing here upgrades one.

| file | what is proved | where it is used |
|---|---|---|
| `Formal/MergeLaw.lean` | merged tank below the saddle ⇔ `j > f k`; least flipping count `⌊f k⌋ + 1`; unanimity ⇔ `f ≥ 1 − 1/k`; one failure flips ⇔ `f < 1/k` | FINDINGS §112 and its regime table |
| `Formal/Schlogl.lean` | `schlogl_consts` factorises the drift as `−k(x−r₁)(x−r₂)(x−r₃)`; outer rails stable, saddle unstable; slopes in closed form; constants positive for positive rails | every element in the repo |
| `Formal/EffectiveExponent.lean` | `A(½−ε)` closed form; `ν_eff` exactly, as a derivative; `|ν_eff − (2 − 5ε/3)| ≤ 7ε²` on `(0, ¼]`; `ν_eff(0.1) = 235/126` | FINDINGS §84/§85.1, rule 21 |
| `Formal/BirthDeath.lean` | the product formula `π n ∝ ∏ λᵢ/μᵢ₊₁` satisfies detailed balance on every edge and `(πQ)ₙ = 0` (inflow = outflow, stated without reference to detailed balance) in the reflecting box; needs only `μ ≠ 0` | `Elem.stage1_stationary` |
| `Formal/Committor.lean` | the birth-death committor closed form is 0 and 1 at the rails and solves the backward equation at every interior state; needs only `λ > 0` | `asymptotic_floor.committor` (§122: φ) |
| `Formal/EdgeInterval.lean` | every real root of stage 2's cubic is positive; its discriminant as a cubic in the Hill level; `Δ(h)/h` concave, so the bistable Hill levels form an interval | `weak_upstream.edge()` (§124.3) |

**Build:** `cd formal && ~/.elan/bin/lake build` (first time: `lake exe cache get` for Mathlib).
**Audit:** `lake env lean scripts/axioms.lean`. Every theorem must depend only on
`[propext, Classical.choice, Quot.sound]` — `sorryAx` anywhere means something was assumed.

`tests/test_formal.py` runs both, and separately checks that the Python helpers each theorem is
about (`Element.m_of_k`, `x_merged`, `schlogl_consts`, `tau`) agree with the proved statements.
A proof about a definition the code does not implement proves nothing about the code.

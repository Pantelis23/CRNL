/-
§124.11 — the saddle Hessian from a Lyapunov equation, and why the transport integrand vanishes there.

At a fixed point of the rate equations p = 0, and the Riccati equation for Z = ∂²W reduces to the
algebraic equation `Aᵀ Z + Z A + Z D Z = 0` (A the drift Jacobian, D the noise matrix).
`prefactor.saddle_lyapunov` takes Z at the saddle as `C⁻¹`, where `A C + C Aᵀ + D = 0`. Proved here,
for any square matrices over a field with Z invertible:

  * `riccati_iff_lyapunov`: the algebraic Riccati equation holds iff `C = Z⁻¹` solves the Lyapunov
    equation;
  * `trace_identity`: the Riccati equation forces `2 tr A + tr (D Z) = 0`. This is the statement that
    the transport integrand `−tr A − ½ tr(D Z)` vanishes at the fixed point, so ln K converges at both
    ends of the optimal path.

No symmetry of Z or D is assumed; none is needed.
-/
import Mathlib

namespace CRNL.RiccatiLyapunov

open Matrix

variable {n : Type*} [Fintype n] [DecidableEq n] {K : Type*} [Field K]

theorem riccati_iff_lyapunov (A D Z : Matrix n n K) (hZ : IsUnit Z.det) :
    Aᵀ * Z + Z * A + Z * D * Z = 0 ↔ A * Z⁻¹ + Z⁻¹ * Aᵀ + D = 0 := by
  have hl : Z⁻¹ * Z = 1 := nonsing_inv_mul Z hZ
  have hr : Z * Z⁻¹ = 1 := mul_nonsing_inv Z hZ
  have cl : ∀ X : Matrix n n K, Z⁻¹ * (Z * X) = X := fun X => by
    rw [← Matrix.mul_assoc, hl, Matrix.one_mul]
  have cr : ∀ X : Matrix n n K, X * (Z * Z⁻¹) = X := fun X => by rw [hr, Matrix.mul_one]
  have cl' : ∀ X : Matrix n n K, Z * (Z⁻¹ * X) = X := fun X => by
    rw [← Matrix.mul_assoc, hr, Matrix.one_mul]
  have cr' : ∀ X : Matrix n n K, X * (Z⁻¹ * Z) = X := fun X => by rw [hl, Matrix.mul_one]
  constructor
  · intro h
    have e := congrArg (fun M => Z⁻¹ * M * Z⁻¹) h
    simp only [Matrix.mul_add, Matrix.add_mul, Matrix.mul_zero, Matrix.zero_mul, Matrix.mul_assoc,
      cl, cr] at e
    rw [← e]; abel
  · intro h
    have e := congrArg (fun M => Z * M * Z) h
    simp only [Matrix.mul_add, Matrix.add_mul, Matrix.mul_zero, Matrix.zero_mul, Matrix.mul_assoc,
      cl', cr'] at e
    rw [← e]; simp only [Matrix.mul_assoc]; abel

theorem trace_identity (A D Z : Matrix n n K) (hZ : IsUnit Z.det)
    (h : Aᵀ * Z + Z * A + Z * D * Z = 0) : 2 * trace A + trace (D * Z) = 0 := by
  have hl : Z⁻¹ * Z = 1 := nonsing_inv_mul Z hZ
  have hr : Z * Z⁻¹ = 1 := mul_nonsing_inv Z hZ
  have h2 := congrArg (fun M => trace (Z⁻¹ * M)) h
  simp only [Matrix.mul_add, trace_add, Matrix.mul_zero, trace_zero] at h2
  have e1 : trace (Z⁻¹ * (Aᵀ * Z)) = trace A := by
    rw [← Matrix.mul_assoc, trace_mul_comm, ← Matrix.mul_assoc, hr, Matrix.one_mul, trace_transpose]
  have e2 : trace (Z⁻¹ * (Z * A)) = trace A := by rw [← Matrix.mul_assoc, hl, Matrix.one_mul]
  have e3 : trace (Z⁻¹ * (Z * D * Z)) = trace (D * Z) := by
    rw [Matrix.mul_assoc Z D Z, ← Matrix.mul_assoc, hl, Matrix.one_mul]
  rw [e1, e2, e3] at h2
  linear_combination h2

end CRNL.RiccatiLyapunov

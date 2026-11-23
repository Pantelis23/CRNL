# What krsolve needs from KernRift (§134.23–134.24)

This is measured on `solver.kr`'s dense kernel. 95% of the top fronts' time is in `partial_lu`, and that runs at about
2.3 GFLOP/s per core on a Ryzen 9 7900X and about 49 GFLOP/s overall on the Core Ultra 9 275HX laptop. The laptop's
ceilings are roughly 200 GFLOP/s for ideal scalar code and about 1 TFLOP/s for AVX2 + FMA. The items are ordered by
expected gain for krsolve.

1. **Keep loop-carried f64 locals in registers.** A variable declared before a loop and updated in it (the GEMM
   accumulators `c00..c11`) is stored to and reloaded from the stack on every iteration (`0xd0(%rsp)`…), so each
   `c = c + a*b` waits on a store-to-load round trip. Locals declared *inside* the loop body already stay in xmm.
   - *Done when:* `gemm_core`'s k-loop has no `(%rsp)` traffic for `c00..c11`.
   - *Expected:* about 3–4× on the top fronts.
2. **Inline small leaf functions**, or add an indexed f64 load/store intrinsic. `ldf(a, i)` and `stf(a, i, v)` compile
   to real calls, and a 4-way unroll written through `ldf` ran 3.6× slower.
   - *Done when:* the `ldf` form of the kernel is no slower than the explicit-pointer form.
3. **Load f64 directly into xmm.** Use `movsd (mem),%xmm` instead of `mov (mem),%rax` followed by `movq %rax,%xmm`.
4. **FMA.** Either an explicit `fma(a, b, c)` builtin emitting `vfmadd231sd`, or opt-in contraction. The explicit form
   is preferred, because FMA changes rounding and CRNL validates against stored cells.
5. **SIMD: `f64x4` (AVX2)**, with load/store, add/mul/FMA and broadcast. With it, a 4×4 register-blocked GEMM could
   reach a large fraction of peak.
   - *Expected:* a further 4–8× on dense fronts.
   - *Done when:* a GEMM micro-benchmark reaches ≥ 30% of one core's AVX2 FMA peak.
6. **Optional: hybrid-CPU awareness.** A query for which CPUs are P-cores and E-cores, plus thread affinity
   (`sched_setaffinity`). §134.24's dynamic chunking covers most of this already.

Items 1–3 are code-generator fixes that help any KernRift program with hot float loops. Items 4–5 are features.

## The GPU path: deferred until after the codegen fixes (decided 2026-11-23)

**Hardware:**
- Desktop: RX 7700 XT / 7800 XT, RDNA 3, MLRift `amdgpu-native`.
- Laptop: RTX 5080 Laptop, 16 GB, sm_120, MLRift `nvptx` under WSL2.

**Precision must be FP64.** FP32 cannot hold the pinned Q-process: the pin is 2^900 and solutions reach about 1e−300.
Consumer FP64 is about 1/32 (RDNA 3) and 1/64 (Blackwell) of FP32, so roughly 0.5–1.2 TFLOP/s on spec sheets. A
vectorised CPU kernel on the laptop would plausibly reach 300–500 GFLOP/s, so the GPU's edge over a *good* CPU
kernel is small. Offload would also only touch the dense top-front LU, which caps the whole-solve speed-up at about
3–4×.

**Registered go/no-go, written before any GPU number exists.** This needs MLRift FP64 mul, add and FMA in `@kernel`
on both targets. Then benchmark an FP64 GEMM on each GPU against the best CPU kernel on the same machine, after the
codegen fixes.
- If the GPU is at least 5× faster, build the top-front offload.
- If not, drop it and record the numbers.

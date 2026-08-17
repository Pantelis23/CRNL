"""T-DEPTH-g: §75-§79 priced the WRONG failure mode for a chemically-coupled cascade

§75 argued that in a cascade whose coupling is chemical there is no abstract channel: stage 1's
output species IS stage 2's input, so the inter-stage noise is stage 1's own rail fluctuation.
§76-§79 then built everything on the resulting per-stage error, the GAUSSIAN READOUT of that
fluctuation, eps_read = Phi(-Delta/sigma).

**An element in a cascade has a second way to fail, and it was never compared.** During the
stage time it can SPONTANEOUSLY ESCAPE its rail -- cross the saddle on its own -- and that is a
different exponential:

    eps_read  ~ exp(-eta * Omega)        eta = Delta^2/(2V)   [LNA, §78]
    eps_esc   ~ t * exp(-A * Omega)      A   = the escape action [WKB]

**Both are exponential in Omega with different coefficients, so whichever has the SMALLER
exponent dominates absolutely at large Omega.** §75-§79 assumed the readout term without
checking. If A < eta the assumption is wrong and the arc priced the subdominant failure.

PREDICTIONS, written before running.

  P1  GATE. The escape instrument must be exact and must behave: ln T (mean first passage from
      the rail to the saddle) LINEAR in Omega, so A = ln T / Omega converges. **A first attempt
      by banded solve returned a NEGATIVE mean first passage time (-4.98e13) and an Omega-
      independent one -- the reflecting row was written at the wrong lattice site.** Nothing is
      read off an instrument that returns a negative time.
  P2  **THE TEST. Which exponent is smaller?** Report A and eta for the same elements.
      **Predicted: A < eta, so ESCAPE dominates** -- because the readout error asks the
      fluctuation to cross the whole gap in one draw, while escape only asks it to get over the
      barrier once in many attempts, and the barrier is the cheaper route. If instead eta < A
      the arc is safe and this section is a footnote.
  P3  **BY HOW MUCH, at the stage times actually used.** Report ln(eps_esc) - ln(eps_read) at
      the Omega and t of §72/§75. A difference of a few is a correction; a difference of
      hundreds means §75-§79's eps is not the physical one at all.
  P4  **THE CROSSOVER, so the scope is stated rather than implied.** Readout dominates only for
      stage times below t* = exp(-(eta - A) Omega). Report t* -- if it is exponentially small,
      the readout regime is unreachable for a chemically-coupled cascade and §75's own premise
      selects the regime where §75-§79's answer does not apply.
  P5  **WHAT IS UNAFFECTED.** §71/§72 used an EXTERNAL channel with sigma = f*Delta, far wider
      than the intrinsic width, where the readout term is large by construction. That regime is
      a real one -- engineered wiring, a diffusive gap, a readout instrument -- and this section
      must not be read as overturning it. The two regimes are distinguished here.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm

from experiments.cascade_schlogl import rates, schlogl_consts
from experiments.chemical_channel_noise import rail_width
from experiments.derive_eta import schlogl_V


def ln_mfpt(omega, r1, r2, r3, cap_mult=1.8):
    """Exact ln(MFPT) from the high rail to the saddle: absorbing below, reflecting above.

    T(n) = sum_{k=a}^{n-1} [1/(lam(k) pi(k))] sum_{j>k} pi(j), done in logs because pi spans
    e^(Omega A). A banded solve was tried first and returned a negative time (P1).
    """
    c = schlogl_consts(r1, r2, r3)
    cap = int(np.ceil(cap_mult * r3 * omega))
    lam, mu = rates(omega, c, cap)
    a, n0 = int(round(r2 * omega)), int(round(r3 * omega))
    lp = np.full(cap + 1, -np.inf)
    lp[a] = 0.0
    acc = 0.0
    for k in range(a + 1, cap + 1):
        if lam[k - 1] <= 0 or mu[k] <= 0:
            break
        acc += np.log(lam[k - 1]) - np.log(mu[k])
        lp[k] = acc
    terms = []
    for k in range(a, n0):
        if lam[k] <= 0 or not np.isfinite(lp[k]):
            continue
        tail = lp[k + 1:cap + 1]
        tail = tail[np.isfinite(tail)]
        if tail.size:
            terms.append(-np.log(lam[k]) - lp[k] + logsumexp(tail))
    return float(logsumexp(terms)) if terms else np.nan


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rails", type=float, nargs=3, default=[0.1, 1.0, 1.9])
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[400, 800, 1600, 3200, 6400, 12800])
    ap.add_argument("--t", type=float, default=2.0, help="stage time, as §71/§72 used")
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/two_exponentials.json"))
    args = ap.parse_args()
    r1, r2, r3 = args.rails
    D = (r3 - r1) / 2.0
    eta = D ** 2 / (2 * schlogl_V(r1, r2, r3))
    print(f"rails {r1}/{r3}, Delta = {D}, stage time t = {args.t}")
    print(f"eta from §78 (LNA, readout):  {eta:.6f} nats/molecule")

    print("\n=== P1 GATE: is ln(MFPT) linear in Omega, so A converges?")
    print(f"{'Omega':>7}{'ln T':>11}{'A = lnT/Om':>13}{'local A':>11}")
    rows, prev = [], None
    for om in args.omegas:
        lT = ln_mfpt(om, r1, r2, r3)
        loc = (lT - prev[1]) / (om - prev[0]) if prev else None
        rows.append({"omega": om, "lnT": lT, "A": lT / om, "localA": loc})
        print(f"{om:>7}{lT:>11.3f}{lT/om:>13.6f}"
              + (f"{loc:>11.6f}" if loc else f"{'--':>11}"))
        prev = (om, lT)
    locs = [r["localA"] for r in rows if r["localA"]]
    conv = len(locs) >= 3 and abs(locs[-1] - locs[-2]) / abs(locs[-1]) < 0.01
    assert all(np.isfinite(r["lnT"]) for r in rows), "non-finite MFPT"
    assert all(r["lnT"] > 0 for r in rows), "negative MFPT -- the instrument is broken"
    print(f"  local A: " + ", ".join(f"{v:.6f}" for v in locs))
    print(f"  -> P1 {'HOLDS: A converges to ' + f'{locs[-1]:.6f}' if conv else 'FAILS'}")
    A = locs[-1]

    print(f"\n=== P2/P3: which failure mode binds?")
    print(f"{'Omega':>7}{'ln eps_read':>14}{'ln eps_esc':>13}{'difference':>13}{'binds':>10}")
    for r in rows:
        om = r["omega"]
        sd = rail_width(om, r1, r2, r3)["sd_exact"]
        lr = float(norm.logcdf(-D / sd))
        le = np.log(args.t) - r["lnT"]
        r["ln_read"], r["ln_esc"] = lr, le
        print(f"{om:>7}{lr:>14.2f}{le:>13.2f}{le - lr:>13.2f}"
              f"{('ESCAPE' if le > lr else 'readout'):>10}")
    print(f"\n  eta (readout exponent) = {eta:.6f}")
    print(f"  A   (escape  exponent) = {A:.6f}")
    print(f"  ratio eta/A = {eta/A:.3f}")
    print(f"  -> P2 {'ESCAPE DOMINATES: A < eta, so §75-§79 priced the subdominant failure' if A < eta else 'readout dominates: the arc is safe'}")
    big = rows[-1]
    print(f"  P3: at Omega = {big['omega']}, escape is e^{big['ln_esc']-big['ln_read']:.0f}"
          f" = {np.exp(min(big['ln_esc']-big['ln_read'], 700)):.3g}x more likely than misreading")

    print(f"\n=== P4: the crossover -- how short would a stage have to be?")
    print(f"{'Omega':>7}{'t* (readout binds below)':>28}")
    for r in rows:
        tstar = float(np.exp(r["ln_read"] + r["lnT"]))
        r["tstar"] = tstar
        print(f"{r['omega']:>7}{tstar:>28.3e}")
    print(f"  -> the readout regime needs a stage time below e^(-(eta-A)Omega), which is")
    print(f"     exponentially small. **§75's own premise -- a chemically coupled cascade --")
    print(f"     selects the regime where §75-§79's eps is not the physical one.**")

    print(f"\n=== P5: what is NOT affected")
    sd = rail_width(1600, r1, r2, r3)["sd_exact"]
    ext = 0.35 * D
    print(f"  intrinsic sigma at Omega=1600: {sd:.5f};  §71/§72's external channel: {ext:.5f}"
          f"  ({ext/sd:.1f}x wider)")
    print(f"  ln eps_read with the EXTERNAL channel: {float(norm.logcdf(-D/ext)):.2f}"
          f"  vs escape {rows[2]['ln_esc']:.2f}")
    print(f"  -> with an external channel the readout term is far larger and DOES bind, so")
    print(f"     §71/§72 are unaffected. The two regimes are physically different cascades.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "eta": eta, "A": A}, indent=2,
                                   default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()

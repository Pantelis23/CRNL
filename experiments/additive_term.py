"""T15-b: does a network with an ADDITIVE drift term break the categorical zero?

§29 and §30 established that AM's categorical `s-only` zero is a theorem, because the
drift of every signal coordinate is exactly proportional to that coordinate. §29 then
made a prediction about networks outside that class:

    "A network whose `b_delta` carries an additive term would give a small but
     NONZERO `s-only`, not a zero."

That prediction has never been tested, and until it is, "the identity is why the zero
happens" is an explanation with no independent confirmation -- rule 17's exact
situation. `am_asymmetric` supplies the test. Working the algebra out by hand from its
propensities, with s = n_X + n_Y and delta = n_X - n_Y:

    d(delta)/dt = delta * (k/Omega) * [ n_B - gamma*(s - 1) ]
                       + (k*beta/Omega) * [ n_B*s - gamma*((s^2 + delta^2)/2 - s) ]

The first term is §29's identity exactly. The second is an ADDITIVE term that exists
for every beta != 0 and vanishes identically at beta = 0. Verified against the
network's own propensities at gamma in {0.10, 0.30, 0.45} x beta in {0, 0.05, 0.20,
0.50}: worst relative residual 3.5e-14, and the additive term's median size relative to
the multiplicative one runs 0.00 / 0.08 / 0.40 / 0.93 across those betas. So beta tunes
the very thing §29 says the zero depends on, in a network already in the repo.

THE OBSERVABLE IS ARMED SO THE ADDITIVE TERM WORKS AGAINST IT. beta > 0 makes X the
better catalyst, and the additive term is positive over most of the simplex, so it
pushes delta UP. The champion is therefore **Y**, the disfavoured symbol: the additive
term erodes Y's lead, which is the direction that can produce errors.

WHY THIS TEST IS IMMUNE TO THE CONFOUND THAT KILLED §30.1. beta also moves the basin
boundary and so lowers the barrier, and §30.2 is a fresh reminder that a magnitude
trend along such an axis proves nothing. **But the claim here is CATEGORICAL, not a
magnitude**: at beta = 0, `s-only` is exactly 0 while `full` fails at a measurable
rate, and no change in barrier height can turn an exact zero into a nonzero -- a
theorem does not weaken, it either applies or does not. The magnitude trend in beta is
reported but explicitly NOT read as mechanism.

TWO START RULES, so the barrier confound is broken by construction rather than argued
away (the lesson of §30.2, applied before the run instead of after):

  * `fixed`   -- delta_0/Omega = -eps*delta_star(gamma) for every beta. The start is
                 identical across beta, so the arms are paired on one state, but the
                 barrier falls as the saddle slides toward the start.
  * `matched` -- delta_0/Omega = saddle(gamma, beta) - eps*delta_star(gamma), so the
                 start keeps a FIXED distance from the basin boundary and the barrier
                 is approximately held. Identical to `fixed` at beta = 0 by
                 construction.

If `s-only` turns nonzero under BOTH rules, the additive term is doing it. If it turns
nonzero only under `fixed`, the barrier is doing it and §29's prediction is unsupported.

HOW THE CELL WAS CHOSEN, disclosed because it was chosen by looking. The first pass
used eps = 0.25, which puts the start 15 counts clear of the saddle; the additive term
there is 5-10x the multiplicative one at the start state and still nothing crossed,
because 15 counts is simply a long way. A probe at eps = 0.06 (3-6 counts of margin,
4,000 trials) found crossings, and that is the eps used below. **The predictions were
fixed before either run and are unchanged** -- what the probe selected is the cell, the
same way every experiment here selects an Omega range. The probe's own numbers are
reported in FINDINGS alongside the confirming run rather than quietly replaced by it.

Note that a small margin does not weaken P2: at beta = 0 the theorem forces `s-only`
to exactly 0 for ANY margin, even one count. That is what makes a near-boundary cell
the sharpest place to look rather than a compromised one.

PREDICTIONS, written before running:

  P1  The additive identity above holds to < 1e-12 relative. (Already checked to
      3.5e-14 while designing the cells; re-run here so the number is in the record
      alongside the result that uses it.)
  P2  At beta = 0, `s-only` is EXACTLY 0 with zero sign flips, reproducing §29 on a
      network built by a different constructor. `am_asymmetric(gamma, 0)` is
      `am_reversible(gamma)` reaction for reaction, so a nonzero here would mean the
      harness, not the chemistry, and would invalidate the whole run.
  P3  THE TEST. At beta > 0, `s-only` is NONZERO under both start rules, and delta
      changes sign in trajectories that could not have crossed at beta = 0. §29's
      prediction stands or falls here.
  P4  THE CONTROL THAT DECIDES ADMISSIBILITY. The deterministic ODE from each start
      must still reach the Y attractor. `am_asymmetric`'s own docstring warns that the
      tilt produces a SYSTEMATIC bias distinct from the random error the wall protects
      against, and that "both show up as a wrong answer". A cell whose ODE already
      fails is measuring the bias, is not a restoration failure at all, and is reported
      as inadmissible rather than scored.
  P5  `full` reproduces the exact CME (via `splitting_probability`, which takes the
      favoured predicate rather than assuming n[0] > n[1]) at every cell, or nothing
      else in the row is admissible.
  P6  `delta-only` continues to carry the answer at 2-18% as in §24.1, at beta = 0 at
      least. I have no prediction for how it behaves at beta > 0 and am not making one.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import propensities_batch
from crnl.cme import splitting_probability
from crnl.deterministic import integrate
from crnl.networks.am_asymmetric import am_asymmetric, basin_boundary
from crnl.networks.am_reversible import delta_star
from crnl.vectorized import compile_network
from experiments.noise_placement import project

ARMS = ("full", "delta-only", "s-only")


def additive_residual(gamma: float, beta: float, omega: int, rng, n: int) -> float:
    comp = compile_network(am_asymmetric(gamma, beta), float(omega))
    worst = 0.0
    for _ in range(n):
        c = np.sort(rng.integers(0, omega + 1, size=2))
        st = np.array([c[0], c[1] - c[0], omega - c[1]], dtype=np.int64)
        a = propensities_batch(comp, st[None, :].astype(float))[0]
        b = comp.S @ a
        s = float(st[0] + st[1]); d = float(st[0] - st[1]); nb = float(st[2])
        mult = d * (nb - gamma * (s - 1.0)) / omega
        add = beta * (nb * s - gamma * ((s * s + d * d) / 2.0 - s)) / omega
        scale = max(abs(mult), abs(add), 1e-300)
        worst = max(worst, abs(float(b[0] - b[1]) - mult - add) / scale)
    return worst


def make_start(gamma: float, beta: float, omega: int, eps: float, rule: str):
    """Start with Y ahead. `fixed` ignores beta; `matched` tracks the saddle."""
    off = eps * delta_star(gamma)
    frac = -off if rule == "fixed" else basin_boundary(gamma, beta) - off
    if not np.isfinite(frac):
        raise ValueError("no saddle (monostable)")
    nb = int(round(omega * gamma / (1.0 + gamma)))
    rest = omega - nb
    d0 = int(round(frac * omega))
    if (rest - d0) % 2:
        d0 -= 1
    nx, ny = (rest + d0) // 2, (rest - d0) // 2
    if min(nx, ny, nb) < 1:
        raise ValueError(f"degenerate start {(nx, ny, nb)}")
    return np.array([nx, ny, nb], dtype=np.int64)


def ode_is_safe(gamma: float, beta: float, n0: np.ndarray, omega: int) -> tuple:
    """P4: does the noiseless flow from this start still reach the Y attractor?"""
    net = am_asymmetric(gamma, beta)
    traj = integrate(net, np.asarray(n0, float) / omega, t_span=(0.0, 400.0))
    fin = traj.final()
    x, y = float(fin[0]), float(fin[1])
    return (y > x), x - y


def run_arm(comp, n0, rng, *, dt, thr, trials, t_max, mode, max_steps=300_000):
    """CLE with projected noise. Champion is Y: an error is n_X > n_Y at absorption."""
    S = comp.S.astype(float)
    m = np.tile(np.asarray(n0, dtype=float), (trials, 1))
    live = np.ones(trials, bool)
    t = np.zeros(trials)
    d0 = float(n0[0] - n0[1])
    sgn = np.sign(d0)
    worst = np.full(trials, abs(d0))
    flipped = np.zeros(trials, bool)
    rejected = steps = 0
    var_full = var_kept = 0.0
    for _ in range(max_steps):
        if not live.any():
            break
        idx = np.where(live)[0]
        a = propensities_batch(comp, m[idx])
        mean = a * dt
        xi = np.sqrt(np.clip(mean, 0.0, None)) * rng.standard_normal(mean.shape)
        nz_full = xi @ S.T
        nz = project(nz_full, mode)
        var_full += float((nz_full ** 2).sum())
        var_kept += float((nz ** 2).sum())
        cand = m[idx] + (mean @ S.T) + nz
        ok = (cand >= 0.0).all(axis=1)
        rejected += int((~ok).sum())
        steps += len(idx)
        upd = idx[ok]
        m[upd] = cand[ok]
        t[upd] += dt
        cur = (m[upd][:, 0] - m[upd][:, 1]) * sgn
        worst[upd] = np.minimum(worst[upd], cur)
        flipped[upd] |= (cur <= 0.0)
        absorbed = np.abs(m[idx, 0] - m[idx, 1]) >= thr
        live[idx[absorbed | (t[idx] >= t_max)]] = False
    fin = np.abs(m[:, 0] - m[:, 1]) >= thr
    nok = int(fin.sum())
    wrong = int((m[fin][:, 0] > m[fin][:, 1]).sum()) if nok else 0
    return {"p": wrong / nok if nok else float("nan"), "wrong": wrong, "n_ok": nok,
            "unfinished": int(trials - nok), "rejected": rejected,
            "reject_frac": rejected / steps if steps else float("nan"),
            "flipped": int(flipped.sum()),
            "closest": float((worst / abs(d0)).min()),
            "variance_kept_frac": var_kept / var_full if var_full > 0 else float("nan")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.25)
    ap.add_argument("--omega", type=int, default=80)
    ap.add_argument("--betas", type=float, nargs="+",
                    default=[0.0, 0.02, 0.05, 0.10, 0.15])
    ap.add_argument("--eps-frac", type=float, default=0.25)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--rules", type=str, nargs="+", default=["fixed", "matched"])
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--trials", type=int, default=60000)
    ap.add_argument("--seed", type=int, default=20260803)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/additive_term.json"))
    args = ap.parse_args()

    t0 = time.time()
    g, om = args.gamma, args.omega
    ds = delta_star(g)
    thr = max(2, int(round(args.theta * ds * om)))
    rng = np.random.default_rng(args.seed)

    print("=== P1: d(delta)/dt = delta*[n_B - g(s-1)]/Om + beta*[n_B*s "
          "- g((s^2+d^2)/2 - s)]/Om")
    ident = {}
    for b in args.betas:
        r = additive_residual(g, b, 97, rng, 60)
        ident[str(b)] = r
        print(f"  beta = {b:<5} worst relative residual {r:.3e}")
    print(f"  -> P1 {'HOLDS' if max(ident.values()) < 1e-12 else 'FAILS'}")

    print(f"\ngamma={g} Omega={om} delta*={ds:.4f} threshold={thr} "
          f"eps={args.eps_frac} trials={args.trials}   champion is Y")
    rows = []
    for rule in args.rules:
        print(f"\n--- start rule: {rule}")
        print(f"{'beta':>6}{'saddle':>9}{'start':>16}{'ODE':>6}{'CME':>9}"
              f"{'full':>9}{'delta':>9}{'s-only':>9}{'s flips':>9}"
              f"{'s closest':>11}{'unfin':>7}")
        for b in args.betas:
            try:
                n0 = make_start(g, b, om, args.eps_frac, rule)
            except ValueError as e:
                print(f"{b:>6.2f}   SKIPPED ({e})")
                continue
            safe, xy_end = ode_is_safe(g, b, n0, om)
            net = am_asymmetric(g, b)
            res = splitting_probability(
                net, om, float(om), n0,
                lambda s: abs(int(s[0]) - int(s[1])) >= thr,
                lambda s: int(s[1]) > int(s[0]))
            exact = 1.0 - res["split"] if res["valid"] else float("nan")
            comp = compile_network(net, float(om))
            got = {}
            for mode in ARMS:
                r = np.random.default_rng(args.seed + 77 + int(1000 * b))
                got[mode] = run_arm(comp, n0, r, dt=args.dt, thr=thr,
                                    trials=args.trials, t_max=6000.0, mode=mode)
            so = got["s-only"]
            print(f"{b:>6.2f}{basin_boundary(g, b):>9.4f}{str(n0.tolist()):>16}"
                  f"{'ok' if safe else 'BIAS':>6}{exact:>9.5f}"
                  f"{got['full']['p']:>9.5f}{got['delta-only']['p']:>9.5f}"
                  f"{so['p']:>9.5f}{so['flipped']:>9}{so['closest']:>11.4f}"
                  f"{so['unfinished']:>7}")
            rows.append({"rule": rule, "beta": b, "start": n0.tolist(),
                         "saddle": float(basin_boundary(g, b)),
                         "ode_safe": bool(safe), "ode_end_xy": float(xy_end),
                         "p_cme": exact, "cme_valid": bool(res["valid"]),
                         **{m: got[m] for m in ARMS}})

    print(f"\n=== P4 admissibility: cells whose noiseless flow already fails")
    bad = [r for r in rows if not r["ode_safe"]]
    print(f"  {len(bad)} of {len(rows)} inadmissible"
          + (f" -> {[(r['rule'], r['beta']) for r in bad]}" if bad else " (none)"))
    adm = [r for r in rows if r["ode_safe"]]

    print(f"\n=== P2/P3: is the categorical zero broken?")
    print(f"{'rule':>9}{'beta':>7}{'s-only P':>11}{'flips':>8}{'full P':>9}"
          f"{'s/full':>9}")
    for r in adm:
        f = r["full"]["p"]
        print(f"{r['rule']:>9}{r['beta']:>7.2f}{r['s-only']['p']:>11.6f}"
              f"{r['s-only']['flipped']:>8}{f:>9.5f}"
              f"{(r['s-only']['p']/f if f else float('nan')):>9.4f}")
    for rule in args.rules:
        rs = [r for r in adm if r["rule"] == rule]
        z = [r for r in rs if r["beta"] == 0.0]
        nz = [r for r in rs if r["beta"] > 0.0]
        if z:
            ok0 = z[0]["s-only"]["p"] == 0.0 and z[0]["s-only"]["flipped"] == 0
            print(f"  [{rule}] beta=0 exactly zero with no flips? "
                  f"{'YES' if ok0 else 'NO -- P2 FAILS'}")
        live = [r for r in nz if r["s-only"]["flipped"] > 0]
        print(f"  [{rule}] beta>0 cells with s-only sign flips: "
              f"{len(live)}/{len(nz)}")
    both = all(
        any(r["s-only"]["flipped"] > 0
            for r in adm if r["rule"] == rule and r["beta"] > 0)
        for rule in args.rules if any(r["rule"] == rule for r in adm))
    print(f"\n  -> P3 {'HOLDS under every start rule' if both else 'does NOT hold under every start rule'}"
          f": §29's prediction that an additive term breaks the categorical zero is "
          f"{'CONFIRMED' if both else 'NOT confirmed'}.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"gamma": g, "omega": om, "thr": int(thr),
                                    "identity": ident, "rows": rows},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()

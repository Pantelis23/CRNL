"""T15-i: is the restoration boundary a CLOSED FORM? — making §56's trichotomy constructive

§56 proved the non-restoring rate vectors form a convex cone and that **capability** is
combinatorial: some d_r > 0 suffices, because one can load c onto that reaction and pick x
where its bracket dominates. It said nothing about **realisation**, and left the trichotomy
half-finished -- "tuned" is not a criterion, it is the absence of one.

**§55 already contains the missing half without saying so.** It measured P at the SYMMETRIC
FIXED POINT and got (1 - 2 gamma)/3 for AM, vanishing exactly at gamma_c = 1/2, and 1/(2n-1)
vanishing at gamma_c(n). Those are the published critical points of the network. So the
state that decides realisation is not an arbitrary accessible x -- **it is the symmetric
steady state, the decision point the dynamics actually occupies**, and §54 gives P there in
closed form:

    RESTORES  <=>  sum_r c_r d_r B_r(x*) > 0,    d_r = S_X(r) - S_Y(r),  B_r(x*) >= 0

with x* the symmetric steady state. **That is a single linear inequality in the rate
constants at fixed x*** -- computable without any simulation, for any exchange-symmetric
mass-action network, because §43 makes P exist and §54 makes it explicit.

**THE CATCH, AND IT IS WHY THIS IS NOT A TAUTOLOGY.** x* itself depends on c. §56's cone
theorem needed P linear in c at FIXED x; here the evaluation point moves with c, so the
argument does not transport. Whether the realising set is still convex is an open structural
question with a one-counterexample kill test, and it is asked below.

The test is against the DYNAMICS, not against another algebraic rule: integrate the full
mass-action ODE from the symmetric steady state with an antisymmetric kick and ask whether
the normalised spread grows or dies.

PREDICTIONS, written before running; per rule 19 each VERDICT rule is stated with the data
that would make it print the other answer.

  P1  GATE (a). §54's closed form <c, v(x*)> equals §53's P_at(x*) to machine precision.
      Both are published; disagreement means one is wrong at the symmetric point and
      nothing below counts.
  P1  GATE (b). For AM the criterion reproduces §55's PUBLISHED (1 - 2 gamma)/3 and hence
      the boundary gamma_c = 1/2 to 1e-10. This is an absolute check against a stored
      number, not a fit (rule 16).
  P2  **THE TEST. sign(P at x*) predicts what the ODE does.** Restoration is the normalised
      spread growing by >= 10x, non-restoration is it dying to <= 1/10. **This is a
      UNIVERSAL claim, so ONE disagreement refutes it** (rule 19, §53's lesson) -- no
      fraction is thresholded, every disagreement is printed with its network and rates.
      Excluded and COUNTED, not dropped: no symmetric fixed point in (0, 1/2); more than
      one; |P| < 1e-8 (genuinely marginal, where a disagreement is not a counterexample);
      spread ending between 1/10 and 10 (ambiguous); any species going negative.
      **AND THE VERDICT IS GATED ON BOTH BRANCHES FIRING.** A first pass drew 20 networks
      of which the criterion predicted 20 decaying and 0 restoring: "P2 HOLDS" printed off
      a branch that never ran, which is the exact failure rule 19 names. So the sampling is
      stratified to fill both branches, and **P2 may not print HOLDS unless at least
      `--minbranch` cases of EACH kind were decided.** Otherwise it prints UNDER-TESTED.
  P3  **THE STRUCTURAL QUESTION. Is the REALISING set convex?** Scaling all rates by t
      leaves x* fixed and scales P by t, so the set is a cone either way; the question is
      addition. **Predicted NO**, because x*(c1 + c2) is not between x*(c1) and x*(c2) in
      any way P respects. Kill test: two restoring c whose sum does not restore. One
      example settles it. **If none is found, that is reported as a failed search with its
      size, NOT as convexity** -- absence of a counterexample is not a proof.
  P4  **RULE 9, an axis I did not choose.** P is a linearisation, so it can only speak
      locally. Sweep the kick size d0 over decades. If some network is locally
      non-restoring yet restores from a large kick, the criterion is LOCAL-ONLY -- that is
      a subcritical/bistable network and a real finding about the criterion's scope, not a
      failure of it. **This one is a measurement, not a pass/fail**, and is reported as the
      fraction of networks whose verdict moves with d0.

SECOND PASS, written before re-running. The first pass's P4 row at d0 = 0.1 read "20
ambiguous, 0 flip", and reading the cells (rule 18) showed **all 20 were restoring networks
whose spread SATURATED at s1 = 0.70-1.00** -- they restored completely -- but whose ratio
fell short of 10 because s0 was already 0.11-0.43. The normalised spread cannot exceed 1, so
the largest attainable ratio is 1/s0, which was **2.3-9.5 in every one of the 20 cases: the
criterion was unreachable by construction**. That is rule 19 inside P4 itself, and a
conservative geometric estimate caught only 5 of the 20 -- the cells had to be read. The
verdict rule is replaced by a CEILING-AWARE one,

    restore  <=>  s1 >= min(10 * s0, 0.5),      decay  <=>  s1 <= s0 / 10

which is satisfiable at every d0 and is unchanged at the small d0 where P2 runs. The
original rule and its failure stay here per rule 3.
  P5  If P2 holds, the criterion is constructive: report gamma_c in closed form for AM (a
      known number, so a gate) and for several other class-sets (numbers this project has
      never had).
  P6  **RULE 10 GUARD: could the harness be doing the restoring?** The integrator must not
      clip species at zero -- a clip at 0 is an absorbing boundary the chemistry does not
      have and would manufacture a winner. Negative excursions are detected and the cell is
      discarded, never clamped. Conservation is checked along every trajectory.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from crnl.networks.am_reversible import am_reversible
from experiments.amplification_sign import P_at
from experiments.optimal_element import symmetric_classes
from experiments.free_rate_optimum import build_free
from experiments.restoration_cone import v_of


def rhs(net, x):
    return net.stoichiometry_matrix() @ net.fluxes(np.asarray(x, dtype=float))


def _sym_state(u):
    """The symmetric line of the conservative simplex: x = y = u, b = 1 - 2u."""
    return np.array([u, u, 1.0 - 2.0 * u])


def symmetric_fixed_points(net, n=4001, lo=1e-9, hi=0.5 - 1e-9):
    """All roots of dX/dt along the symmetric line. Conservation makes dB/dt follow."""
    us = np.linspace(lo, hi, n)
    f = np.array([rhs(net, _sym_state(u))[0] for u in us])
    roots = []
    for k in range(n - 1):
        if not np.isfinite(f[k]) or not np.isfinite(f[k + 1]):
            continue
        if f[k] == 0.0:
            roots.append(us[k])
        elif f[k] * f[k + 1] < 0:
            roots.append(brentq(lambda u: rhs(net, _sym_state(u))[0],
                                us[k], us[k + 1], xtol=1e-14, rtol=1e-14))
    out = []
    for r in roots:
        if all(abs(r - o) > 1e-9 for o in out):
            out.append(float(r))
    return out


def P_closed(net, u):
    """§54's decomposition evaluated on the symmetric line: sum_r c_r d_r B_r(x*)."""
    v, c = v_of(net, _sym_state(u), "X", "Y")
    if v is None or v.size == 0:
        return None
    return float(np.dot(c, v))


def dynamics_verdict(net, u, d0, T=4000.0, grow=10.0):
    """Integrate the FULL ODE from x* + antisymmetric kick. No clipping (P6)."""
    x0 = _sym_state(u) + np.array([d0 / 2.0, -d0 / 2.0, 0.0])
    if x0.min() < 0:
        return None
    tot0 = x0.sum()
    s0 = abs(x0[0] - x0[1]) / (x0[0] + x0[1])
    neg = {"hit": False}

    def f(t, x):
        return rhs(net, x)

    sol = solve_ivp(f, (0.0, T), x0, method="LSODA", rtol=1e-11, atol=1e-14,
                    dense_output=False, t_eval=np.linspace(0, T, 400))
    if not sol.success:
        return None
    Y = sol.y
    if Y.min() < -1e-9:
        neg["hit"] = True                       # P6: discard, never clamp
    if abs(Y.sum(axis=0) - tot0).max() > 1e-7:
        return None                             # conservation broken -> instrument fault
    den = Y[0] + Y[1]
    if den.min() <= 1e-12:
        return None
    s = np.abs(Y[0] - Y[1]) / den
    ratio = float(s[-1] / s0)
    # ceiling-aware: s <= 1 always, so `ratio >= grow` is unreachable once s0 > 1/grow
    up = min(grow * s0, 0.5)
    return {"s0": float(s0), "s1": float(s[-1]), "ratio": ratio, "neg": neg["hit"],
            "ceiling": float(1.0 / s0),
            "verdict": ("restore" if s[-1] >= up
                        else "decay" if ratio <= 1.0 / grow else "ambiguous")}


def am_from_gamma(gamma):
    return am_reversible(gamma)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nets", type=int, default=140)
    ap.add_argument("--seed", type=int, default=20260812)
    ap.add_argument("--d0", type=float, default=1e-4)
    ap.add_argument("--d0s", type=float, nargs="+",
                    default=[1e-6, 1e-4, 1e-2, 0.1])
    ap.add_argument("--minbranch", type=int, default=20,
                    help="P2 declares no verdict unless BOTH branches have this many")
    ap.add_argument("--pairs", type=int, default=4000)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/restoration_boundary.json"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    t0 = time.time()

    # ---------------- P1 GATE (a): §54's closed form vs §53's P_at ------------------
    print("=== P1 GATE (a): <c, v(x*)> (§54) vs P_at(x*) (§53)")
    worst_a = 0.0
    cls_all = symmetric_classes()
    for _ in range(40):
        k = rng.integers(1, 4)
        ids = rng.choice(len(cls_all), size=int(k), replace=False)
        rates = np.exp(rng.uniform(-2, 2, 2 * len(ids)))
        net = build_free([cls_all[i] for i in ids], rates)
        u = float(rng.uniform(0.05, 0.45))
        a = P_closed(net, u)
        if a is None:
            continue
        b = P_at(net, _sym_state(u), 0, 1, h=1e-6)
        worst_a = max(worst_a, abs(a - b) / max(abs(a), abs(b), 1e-12))
    print(f"  worst relative disagreement over 40 networks: {worst_a:.3e}")
    print(f"  -> P1(a) {'HOLDS' if worst_a < 1e-6 else 'FAILS'}")

    # ---------------- P1 GATE (b): AM must reproduce §55's (1-2g)/3 -----------------
    print("\n=== P1 GATE (b): AM against §55's published P(x*) = (1 - 2 gamma)/3")
    print(f"{'gamma':>8}{'u*':>12}{'P_closed':>14}{'(1-2g)/3':>12}{'rel':>11}")
    worst_b = 0.0
    for g in (0.05, 0.20, 0.35, 0.49):
        net = am_from_gamma(g)
        fps = symmetric_fixed_points(net)
        if len(fps) != 1:
            print(f"{g:>8.2f}   {len(fps)} symmetric fixed points -- SKIPPED")
            continue
        u = fps[0]
        p = P_closed(net, u)
        ref = (1.0 - 2.0 * g) / 3.0
        rel = abs(p - ref) / max(abs(ref), 1e-14)
        worst_b = max(worst_b, rel)
        print(f"{g:>8.2f}{u:>12.8f}{p:>14.9f}{ref:>12.8f}{rel:>11.2e}")
    print(f"  -> P1(b) {'HOLDS' if worst_b < 1e-10 else 'FAILS'} (worst {worst_b:.2e})")

    # gamma_c for AM by root-finding the criterion -- must land on 1/2 exactly
    def am_P(g):
        net = am_from_gamma(g)
        fps = symmetric_fixed_points(net)
        return P_closed(net, fps[0]) if len(fps) == 1 else np.nan

    gc = brentq(am_P, 0.05, 0.95, xtol=1e-14, rtol=1e-15)
    print(f"  gamma_c from the criterion = {gc:.14f}   published 0.5"
          f"   |diff| {abs(gc - 0.5):.2e}")

    # ---------------- P2: the criterion against the DYNAMICS ------------------------
    print(f"\n=== P2: does sign(P at x*) predict the ODE? (one disagreement refutes)")
    idx0 = {c[0]: i for i, c in enumerate(cls_all)}
    DIS0 = idx0[(("X", "Y"), ("B", "B"))]
    REC0 = idx0[(("B", "X"), ("X", "X"))]

    def draw(force_am):
        """A class-set and rates. `force_am` seeds AM's two classes, which §56 says are
        capable -- without it the restoring branch is never populated."""
        k = int(rng.integers(1, 4))
        if force_am:
            extra = [i for i in rng.choice(len(cls_all), size=k, replace=False).tolist()
                     if i not in (DIS0, REC0)][: max(0, k - 2)]
            ids = sorted(set([DIS0, REC0] + extra))
        else:
            ids = sorted(rng.choice(len(cls_all), size=k, replace=False).tolist())
        return ids, np.exp(rng.uniform(-3, 3, 2 * len(ids)))

    rows, excl = [], {"no_fp": 0, "multi_fp": 0, "marginal": 0, "ambiguous": 0,
                      "neg_species": 0, "solver": 0, "branch_full": 0}
    disagree = []
    tried = 0
    quota = {"restore": args.nets // 2, "decay": args.nets - args.nets // 2}
    have = {"restore": 0, "decay": 0}
    while len(rows) < args.nets and tried < 60 * args.nets:
        tried += 1
        ids, rates = draw(force_am=(tried % 2 == 0))
        net = build_free([cls_all[i] for i in ids], rates)
        try:
            fps = symmetric_fixed_points(net, n=1501)
        except Exception:
            excl["solver"] += 1
            continue
        if not fps:
            excl["no_fp"] += 1
            continue
        if len(fps) > 1:
            excl["multi_fp"] += 1
            continue
        u = fps[0]
        p = P_closed(net, u)
        if p is None:
            excl["solver"] += 1
            continue
        if abs(p) < 1e-8:
            excl["marginal"] += 1
            continue
        branch = "restore" if p > 0 else "decay"
        if have[branch] >= quota[branch]:
            excl["branch_full"] += 1
            continue
        try:
            dyn = dynamics_verdict(net, u, args.d0)
        except Exception:
            excl["solver"] += 1
            continue
        if dyn is None:
            excl["solver"] += 1
            continue
        if dyn["neg"]:
            excl["neg_species"] += 1
            continue
        if dyn["verdict"] == "ambiguous":
            excl["ambiguous"] += 1
            continue
        pred = "restore" if p > 0 else "decay"
        agree = pred == dyn["verdict"]
        r = {"classes": ids, "rates": rates.tolist(), "u": u, "P": p,
             "pred": pred, "dyn": dyn["verdict"], "ratio": dyn["ratio"],
             "agree": bool(agree)}
        rows.append(r)
        have[pred] += 1
        if not agree:
            disagree.append(r)

    n_r = have["restore"]
    print(f"  {len(rows)} decidable networks ({n_r} predicted restoring,"
          f" {len(rows)-n_r} predicted decaying) out of {tried} drawn")
    print(f"  excluded and counted: " + ", ".join(f"{k}={v}" for k, v in excl.items()))
    print(f"  DISAGREEMENTS: {len(disagree)}")
    for d in disagree[:8]:
        print(f"    classes {d['classes']} P={d['P']:+.4e} pred {d['pred']}"
              f" but ODE {d['dyn']} (spread x{d['ratio']:.3e})")
    both = min(n_r, len(rows) - n_r)
    if disagree:
        print(f"  -> P2 REFUTED by the cases above")
    elif both < args.minbranch:
        print(f"  -> P2 UNDER-TESTED: the smaller branch has only {both} cases"
              f" (< {args.minbranch}); no verdict is declared")
    else:
        print(f"  -> P2 HOLDS: the criterion IS the boundary,"
              f" with both branches exercised ({n_r} / {len(rows)-n_r})")

    # ---------------- P3: is the REALISING set convex? ------------------------------
    print(f"\n=== P3: two restoring rate vectors whose SUM does not restore?")

    def restores_c(ids, rates):
        net = build_free([cls_all[i] for i in ids], rates)
        fps = symmetric_fixed_points(net, n=1501)
        if len(fps) != 1:
            return None, None
        p = P_closed(net, fps[0])
        if p is None or abs(p) < 1e-10:
            return None, None
        return p > 0, p

    found, trials, both_restore = [], 0, 0
    while len(found) < 4 and trials < args.pairs:
        trials += 1
        ids, _ = draw(force_am=(trials % 2 == 0))
        if len(ids) < 2:
            continue
        c1 = np.exp(rng.uniform(-3, 3, 2 * len(ids)))
        c2 = np.exp(rng.uniform(-3, 3, 2 * len(ids)))
        r1, p1 = restores_c(ids, c1)
        r2, p2 = restores_c(ids, c2)
        if not (r1 and r2):
            continue
        both_restore += 1
        rs, ps = restores_c(ids, c1 + c2)
        if rs is False:
            found.append({"classes": ids, "c1": c1.tolist(), "c2": c2.tolist(),
                          "P1": p1, "P2": p2, "Psum": ps})
            print(f"  COUNTEREXAMPLE: classes {ids}"
                  f"  P(c1)={p1:+.3e}  P(c2)={p2:+.3e}  P(c1+c2)={ps:+.3e}")
    print(f"  drew {trials} pairs, of which {both_restore} had BOTH restoring"
          f" (only those can test convexity); found {len(found)} counterexamples")
    if found:
        print(f"  -> P3 the realising set is NOT convex -- §56's cone does not transport")
    elif both_restore < 50:
        print(f"  -> P3 INCONCLUSIVE: only {both_restore} testable pairs")
    else:
        print(f"  -> P3 no counterexample in {both_restore} testable pairs."
              f" NOT a proof of convexity -- a failed search, reported with its size")

    # ---------------- P4: rule 9 -- does the verdict move with the kick size? -------
    print(f"\n=== P4 (rule 9): is the criterion LOCAL only? verdict vs kick size d0")
    sub = ([r for r in rows if r["pred"] == "restore"][:20]
           + [r for r in rows if r["pred"] == "decay"][:20])
    print(f"  {sum(1 for r in sub if r['pred']=='restore')} restoring +"
          f" {sum(1 for r in sub if r['pred']=='decay')} decaying networks re-integrated")
    print(f"{'d0':>10}{'agree':>9}{'ambig':>8}{'flip':>7}{'dropped':>9}")
    d0tab = {}
    for d0 in args.d0s:
        ok = amb = flip = drop = 0
        for r in sub:
            net = build_free([cls_all[i] for i in r["classes"]], np.array(r["rates"]))
            try:
                dyn = dynamics_verdict(net, r["u"], d0)
            except Exception:
                drop += 1
                continue
            if dyn is None or dyn["neg"]:
                drop += 1
                continue
            if dyn["verdict"] == "ambiguous":
                amb += 1
                continue
            if dyn["verdict"] == r["pred"]:
                ok += 1
            else:
                flip += 1
        d0tab[d0] = {"agree": ok, "ambiguous": amb, "flip": flip, "dropped": drop}
        print(f"{d0:>10.0e}{ok:>9}{amb:>8}{flip:>7}{drop:>9}")
    fl = [v["flip"] for v in d0tab.values()]
    print(f"  -> {'verdict is d0-INDEPENDENT: the criterion is global on this family' if max(fl) == 0 else 'some verdicts MOVE with d0: the criterion is local, those networks are subcritical'}")

    # ---------------- P5: closed-form critical points -------------------------------
    print(f"\n=== P5: gamma_c in closed form, for AM and for other class-sets")
    idx_of = {c[0]: i for i, c in enumerate(cls_all)}
    DIS = idx_of[(("X", "Y"), ("B", "B"))]
    REC = idx_of[(("B", "X"), ("X", "X"))]
    named = {"AM {dis,rec}": [DIS, REC]}
    for i in range(len(cls_all)):
        if i in (DIS, REC):
            continue
        named[f"AM+cls{i}"] = [DIS, REC, i]
    print(f"{'network':>16}  {'gamma_c':>12}{'u* at g_c':>12}   note")
    gcs = {}
    for name, ids in list(named.items())[:9]:
        def crit(g, ids=ids):
            rates = []
            for _ in ids:
                rates += [1.0, g]
            net = build_free([cls_all[i] for i in ids], np.array(rates))
            fps = symmetric_fixed_points(net, n=1501)
            if len(fps) != 1:
                return np.nan
            p = P_closed(net, fps[0])
            return np.nan if p is None else p
        try:
            lo, hi = 1e-4, 0.999
            a, b = crit(lo), crit(hi)
            if not (np.isfinite(a) and np.isfinite(b)):
                print(f"{name:>16}  {'--':>12}{'--':>12}   no unique symmetric fixed point")
                continue
            if a * b > 0:
                where = "restores at EVERY gamma" if a > 0 else "restores at NO gamma"
                print(f"{name:>16}  {'--':>12}{'--':>12}   {where}"
                      f"  (P: {a:+.3e} .. {b:+.3e})")
                continue
            g = brentq(crit, lo, hi, xtol=1e-13, rtol=1e-14)
            net = build_free([cls_all[i] for i in ids],
                             np.array([v for _ in ids for v in (1.0, g)]))
            u = symmetric_fixed_points(net, n=1501)[0]
            gcs[name] = g
            print(f"{name:>16}  {g:>12.8f}{u:>12.8f}")
        except Exception as e:
            print(f"{name:>16}   failed ({type(e).__name__})")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"gate_a": worst_a, "gate_b": worst_b, "gamma_c_AM": gc,
         "rows": rows, "excluded": excl, "disagreements": disagree,
         "nonconvex": found, "convex_trials": trials,
         "convex_testable_pairs": both_restore,
         "d0_table": {str(k): v for k, v in d0tab.items()},
         "gamma_c": gcs}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()

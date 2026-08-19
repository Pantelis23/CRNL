"""§99.1 — is §98 a re-parameterisation of §91?

NOT PRE-REGISTERED, and this docstring does not pretend otherwise (rule 2). The
check below was run ad hoc against §98's stored results after a literature report
raised the objection; the numbers were known before this file existed. It is kept
so the discrimination is reproducible, not so it can be read as a prediction.

THE OBJECTION. §91 parameterised the composition penalty by the noise margin in
units of the upstream rail width; §98 parameterised it by the barrier depth A*Omega.
Near a saddle-node these are not independent -- the margin in rail-width units goes
as mu^(3/4)/sqrt(D) and the WKB action as (margin/sigma)^2, so A*Omega ~ (margin/sigma)^2
is expected on general grounds. Worse, WITHIN one element margin/sigma ~ sqrt(Omega),
so §98's Omega-sweep moves both variables together and cannot separate them at all.
On the calibration sweep alone §98 is §91 rewritten.

WHAT SEPARATES THEM. The proportionality constant is element-dependent. §97's element
has a different landscape and a different coupling, so (margin/sigma)^2 / (A*Omega)
differs between the two elements; interpolating §98's calibration curve at §97's
element in each candidate variable therefore gives two different predictions, and the
measured position picks one.

WHAT THIS CANNOT SETTLE. One out-of-sample point against two variables that agree to
22%. Neither error is compared to an error bar, because none has been measured for the
position. The pre-registered version needs a THIRD element chosen so the ratio differs
from both of these -- see THEORIES T-CASC-l.
"""

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(HERE, "results", "penalty_interpolation.json")
OUT = os.path.join(HERE, "results", "margin_vs_action.json")

# Both margins are quoted at Omega = 30, the reference size at which §94 and §97
# measured them. A is the escape action of each element's own landscape.
A_CAL, MARGIN_CAL = 0.190241, 3.199   # calibration element, §94
A_NEW, MARGIN_NEW = 0.439977, 5.382   # §97's element, barrier 2.3x deeper
OM_REF = 30.0


def load():
    with open(SRC) as fh:
        return json.load(fh)


def ratios():
    """(margin/sigma)^2 / (A*Omega) for each element. Equal => the two axes are one."""
    return (MARGIN_CAL**2 / (A_CAL * OM_REF), MARGIN_NEW**2 / (A_NEW * OM_REF))


def discriminate(data):
    """Interpolate the calibration curve at §97's element in each candidate variable."""
    rows = data["sweep"]
    aoms = np.array([r["AOm"] for r in rows], float)
    poss = np.array([r["pos"] for r in rows], float)
    order = np.argsort(aoms)
    aoms, poss = aoms[order], poss[order]
    assert np.all(np.diff(aoms) > 0), "A*Omega must increase for np.interp (§98.1)"

    # Within one element margin/sigma ~ sqrt(Omega), so (margin/sigma)^2 ~ Omega.
    m2s = MARGIN_CAL**2 * np.array([r["omega"] for r in rows], float)[order] / OM_REF
    measured = float(data["new"]["pos"])

    out = []
    for name, xs, xnew in (
        ("A*Omega", aoms, A_NEW * OM_REF),
        ("(margin/sigma)^2", m2s, MARGIN_NEW**2),
    ):
        inside = bool(xs.min() <= xnew <= xs.max())
        pred = float(np.interp(xnew, xs, poss)) if inside else float("nan")
        out.append(
            {
                "variable": name,
                "x_new": float(xnew),
                "inside": inside,
                "predicted": pred,
                "measured": measured,
                "rel_err": (pred - measured) / measured if inside else float("nan"),
            }
        )
    return out


def main():
    data = load()
    r_cal, r_new = ratios()
    print(f"(margin/sigma)^2 / (A*Omega):  calibration {r_cal:.4f}   §97 {r_new:.4f}"
          f"   ({100 * abs(r_new / r_cal - 1):.0f}% apart)")
    if abs(r_new / r_cal - 1) < 0.02:
        print("  the two axes are indistinguishable across elements; this test is void")

    res = discriminate(data)
    print(f"\n{'variable':>20}{'x at §97':>11}{'predicted':>11}{'measured':>10}{'error':>9}")
    for r in res:
        note = "" if r["inside"] else "  (OUTSIDE traced range -- no prediction)"
        print(f"{r['variable']:>20}{r['x_new']:>11.2f}{r['predicted']:>11.4f}"
              f"{r['measured']:>10.4f}{100 * r['rel_err']:>8.1f}%{note}")

    best = min((r for r in res if r["inside"]), key=lambda r: abs(r["rel_err"]))
    print(f"\ntransfers better: {best['variable']}  (post-hoc; one point; see docstring)")

    with open(OUT, "w") as fh:
        json.dump({"ratios": {"calibration": r_cal, "new": r_new},
                   "discrimination": res, "better": best["variable"]}, fh, indent=2)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

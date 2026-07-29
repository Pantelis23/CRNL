# Working rules for CRNL

This project's value is epistemic. It has produced several confident wrong results,
and the record of *how* they were wrong is worth more than the ones that survived.
These rules exist because each was bought with a withdrawn claim.

## The documents are part of the deliverable, not a report on it

`FINDINGS.md` holds what was measured. `THEORIES.md` holds what we think might be
true, what is open, and what we were wrong about. They are not written up
afterwards — **a result that is not in them does not exist**, and a run whose
numbers are only in a terminal is a run that will be redone.

1. **Every measurement lands in `FINDINGS.md` in the same session it is run**,
   with its scope stated (which network, which γ, which Ω range, how many trials).
   Scope is not a caveat paragraph, it is part of the claim.
2. **Every experiment writes its predictions into its own docstring BEFORE it
   runs.** A prediction recorded after seeing the output is not a prediction.
3. **Refuted predictions stay, with the mechanism that refuted them.** Do not
   quietly drop one and report only what worked. Four failed predictions in a
   single session taught more than the confirmations.
4. **When a result is corrected or withdrawn, the original stays visible** — as a
   quoted block or a struck heading — with what killed it. `THEORIES.md` §4 is the
   catalogue and it earns its length.
5. **Update `THEORIES.md` in the same commit that changes what is known.** Closing
   a question means striking it and pointing at the section; opening one means
   naming its kill test. An open question without a kill test is a wish.
6. **Run `experiments/verify_base.py` (it is in the suite) before committing
   anything that touches `crnl/`.** The tests prove the code agrees with itself;
   the audit proves it still agrees with what is written down. Those come apart
   whenever a behavioural function changes under an already-published section.
7. **Never edit a published number to match a new run.** Either the old number
   stands with the new one beside it, or it is explicitly withdrawn. §12's table is
   left as first printed with the refit in §15.2 for exactly this reason.

## Measurement discipline

8. **Verify twice before moving on.** Re-check the result and reason out the
   prediction before running the next thing.
9. **Constancy along the axis you happened to sweep is not constancy.** Before
   merging two anomalies into one cause, measure each along an axis you did not
   choose for it. This has failed three times.
10. **Watch for the harness doing work the chemistry cannot** — a free `sign()`, a
    mismatched control rail, a reset blank pool, a species clipped at zero.
    Three results were withdrawn for this. Ask of every helper: could the
    chemistry have done this by itself?
11. **A control must share a clock with its arm.** Different reaction orders mean
    different time units; comparing raw rates across them is meaningless (§10.3,
    §20.1).
12. **A conditional mean is only as good as the thing conditioned on being rare.**
    Do not summarise cells that were mostly censored.
13. **An approximation's own numerical parameter is a second axis.** Check
    convergence *within* a level before comparing *between* levels.
14. **Retractions are claims too.** Verify a withdrawal as carefully as an
    assertion — one withdrawal here rested on a 2σ signal and was itself wrong.
15. **Report every candidate extrapolation, never only the flattering one.** If
    the ansätze disagree, the quantity is unresolved and says so.

## Mechanics

- Shared venv: `/home/pantelis/Desktop/Projects/Work/venv`.
- **No coauthor trailers in commits.** Hard rule.
- Do not commit `docs/superpowers/` or `.superpowers/`.
- Run the full suite before committing.
- **Pushing is the user's decision.** Ask; do not infer it from a previous push.
- Prefer the Write/Edit tools over shell heredocs for content with backticks, and
  never `pkill -f <pattern>` when the pattern appears in your own command line —
  both have corrupted work here.

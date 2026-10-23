# Checks of the derivations

Symbolic and enumeration checks of results derived in the appendices of the
paper, each done twice: in SymPy (the `.py` scripts) and in Mathematica
(`checks.nb`). Each check prints PASS or FAIL; each script includes a
control, a deliberately wrong answer which must be rejected.

| Script | Results checked |
|---|---|
| `partition_prior.py` | The prior on a partition: the integral over p, the conditional of p, the number of partitions N_l(k), the prior on the number of changes K (against enumerating every partition), and the values of Table tab:changes |
| `triggered_report.py` | A death which brings about the report: the probability of not being reported, M; of being new to a report and shown as alive; of being shown as dead given new; and that the deaths of cases shown as alive are unchanged (by summing over every value of the kind, T_C and T_D) |

To run the SymPy checks: `pip install sympy==1.14.0`, then
`python partition_prior.py` and `python triggered_report.py` from this
folder; a script exits with status 1 if any check fails. To run the
Mathematica checks, open `checks.nb` and evaluate the notebook.

"""Checks of "A death which brings about the report" (Appendix, "What a report
shows: one onset day"): the probability that a fatal case has not been
reported, M (Step 1); that it is new to a report and shown as alive (Step 2);
that it is shown as dead given that it is new (Step 3); and that the deaths
of cases already shown as alive are unchanged.

Each fatal case has a delay to its report T_C and a time from onset to its
death being recorded T_D, independent, and is of the triggered kind with
probability tau, independently of both; a triggered case is reported at
min(T_C, T_D), any other at T_C. The probabilities are checked by summing
over every value of (kind, T_C, T_D), with symbolic probabilities for each
value of T_C and of T_D on a few days, for every pair of reports s < u.
"""

from itertools import product

from sympy import simplify, symbols

from common import check, finish

tau = symbols("tau", positive=True)
days = range(0, 6)
# symbolic probabilities of T_C = 0..5 and T_D = 1..6, each summing to one
c = list(symbols("c0:5", positive=True))
c.append(1 - sum(c))
d = list(symbols("d1:6", positive=True))
d.append(1 - sum(d))
TC = dict(zip(range(0, 6), c))
TD = dict(zip(range(1, 7), d))


def S_C(x):
    """Pr(T_C > x)."""
    return sum(pr for t, pr in TC.items() if t > x)


def S_D(x):
    """Pr(T_D > x)."""
    return sum(pr for t, pr in TD.items() if t > x)


def prob(event):
    """The probability of event(kind, t_C, t_D), summing over all values;
    kind 1 is triggered."""
    total = 0
    for kind, (tc, pc), (td, pd) in product((0, 1), TC.items(), TD.items()):
        if event(kind, tc, td):
            total += (tau if kind else 1 - tau) * pc * pd
    return total


def reported_at(kind, tc, td):
    return min(tc, td) if kind else tc


def M(x):
    """The formula of Step 1."""
    return S_C(x) * ((1 - tau) + tau * S_D(x))


ok1 = ok2 = ok3 = ok4 = True
for s, u in product(days, days):
    if s >= u:
        continue
    # Step 1: not reported by s
    ok1 &= simplify(prob(lambda k, tc, td: reported_at(k, tc, td) > s) - M(s)) == 0
    # Step 2: new to the report at u (not reported by s, reported by u), and
    # shown as alive there (death not recorded by u)
    new_alive = prob(lambda k, tc, td: s < reported_at(k, tc, td) <= u and td > u)
    ok2 &= simplify(new_alive - (S_C(s) - S_C(u)) * S_D(u)) == 0
    # Step 3: shown as dead given new, with Pr(new) = M(s) - M(u)
    new = prob(lambda k, tc, td: s < reported_at(k, tc, td) <= u)
    new_dead = prob(lambda k, tc, td: s < reported_at(k, tc, td) <= u and td <= u)
    ok3 &= simplify(new - (M(s) - M(u))) == 0
    ok3 &= simplify(new_dead / new - (1 - (S_C(s) - S_C(u)) * S_D(u) / (M(s) - M(u)))) == 0
    # a case reported by s and shown as alive there has its death recorded by
    # u with probability 1 - S_D(u)/S_D(s), whichever its kind, so the
    # distribution of the deaths newly shown among such cases is unchanged
    for kind in (0, 1):
        shown_alive = lambda k, tc, td: k == kind and reported_at(k, tc, td) <= s and td > s
        recorded = prob(lambda k, tc, td: shown_alive(k, tc, td) and td <= u)
        ok4 &= simplify(recorded / prob(shown_alive) - (1 - S_D(u) / S_D(s))) == 0

check("Step 1: Pr(not reported by t') = S^C [(1 - tau) + tau S^D]", ok1)
check("Step 2: Pr(new to report j and alive) = [S^C(t_{j-1}) - S^C(t_j)] S^D(t_j)", ok2)
check("Step 3: Pr(new) = M(t_{j-1}) - M(t_j), and Pr(dead | new) as given", ok3)
check("Deaths of cases shown as alive: 1 - S^D(t_j)/S^D(t_{j-1}) for either kind", ok4)

# a control: the formula of Step 2 without S^D, which counts cases whose death
# has been recorded, is rejected
s, u = 1, 3
check("Step 2 (control): a wrong formula, without S^D(t_j), is rejected",
      simplify(prob(lambda k, tc, td: s < reported_at(k, tc, td) <= u and td > u)
               - (S_C(s) - S_C(u))) != 0)

finish()

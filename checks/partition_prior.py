"""Checks of the prior on a partition (Appendix, "Quantities which change at
unknown times"): the integral over p of the joint prior, eq. (app-partition-
joint); the conditional of p, eq. (app-p); the marginal prior of a partition,
eq. (app-partition-marginal); the number of partitions N_l(k); the prior on
the number of changes K; and the values of Table tab:changes.

The joint prior is
    pi(p, P)  propto  Beta(p | a, b) p^K (1 - p)^(T - 1 - K) 1{P in A_l},
with K = K(P) the number of changes of P. Below, m = K and n = T - 1 - K.
"""

from itertools import combinations

from sympy import (beta, binomial, gamma, gammasimp, hyperexpand, integrate,
                   powsimp, simplify, symbols)

from common import check, finish

a, b, p = symbols("a b p", positive=True)
m, n = symbols("m n", nonnegative=True)

# the kernel of the beta density times the Bernoulli terms, the powers of p
# and of 1 - p combined; the density's constant 1/B(a, b) stays outside
kernel = powsimp(p**(a - 1) * (1 - p)**(b - 1) * p**m * (1 - p)**n)
check("Step 1: the powers combine to p^(a+K-1) (1-p)^(b+T-1-K-1)",
      simplify(kernel - p**(a + m - 1) * (1 - p)**(b + n - 1)) == 0)

# Step 2: the integral over p of the joint prior is B(a + K, b + T - 1 - K) / B(a, b).
# SymPy returns the integral as a hypergeometric function at 1, which
# hyperexpand evaluates by Gauss's theorem
integral = hyperexpand(integrate(kernel, (p, 0, 1), conds="none")) / beta(a, b)
check("Step 2: integral over p equals B(a+K, b+T-1-K)/B(a,b)",
      gammasimp(integral.rewrite(gamma) - (beta(a + m, b + n) / beta(a, b)).rewrite(gamma)) == 0)
# a control: a wrong answer (one more change) is not accepted
check("Step 2 (control): a wrong answer, B(a+K+1, b+T-1-K)/B(a,b), is rejected",
      gammasimp(integral.rewrite(gamma) - (beta(a + m + 1, b + n) / beta(a, b)).rewrite(gamma)) != 0)

# Step 1: the conditional of p given P, the joint prior over its integral, is
# the Beta(a + K, b + T - 1 - K) density
conditional = (kernel / beta(a, b)) / integral
target = p**(a + m - 1) * (1 - p)**(b + n - 1) / beta(a + m, b + n)
check("Step 1: conditional of p is Beta(a+K, b+T-1-K)",
      gammasimp(conditional.rewrite(gamma) - target.rewrite(gamma)) == 0)


def partitions(T, l, l1):
    """Every partition of days 1..T whose first segment has at least l1 days
    and every other at least l, as the list of the start days after day 1."""
    out = []
    for K in range(T):
        for starts in combinations(range(2, T + 1), K):
            lengths = [e - s for s, e in zip((1,) + starts, starts + (T + 1,))]
            if lengths[0] >= l1 and all(x >= l for x in lengths[1:]):
                out.append(starts)
    return out


def count(T, l, l1, k):
    """N(k): partitions into k + 1 segments, the first of at least l1 days
    and the others of at least l; with l1 = l, this is N_l(k) of Step 3."""
    r = T - l1 - l * k
    return binomial(r + k, k) if r >= 0 else 0


# Step 3: N_l(k) = binom(T - l(k+1) + k, k), and with a longer first segment
for T, l, l1 in [(12, 1, 1), (12, 2, 2), (14, 3, 3), (14, 3, 6), (15, 2, 5)]:
    found = partitions(T, l, l1)
    ok = all(sum(1 for P in found if len(P) == k) == count(T, l, l1, k) for k in range(T))
    check(f"Step 3: number of partitions with k changes (T={T}, l={l}, first >= {l1})", ok)


def prior_K(T, l, l1, a_, b_):
    """Pr(K = k) from the formula of Step 3 (with the first segment of at least l1 days)."""
    w = [count(T, l, l1, k) * beta(a_ + k, b_ + T - 1 - k) for k in range(T)]
    total = sum(w)
    return [x / total for x in w]


# Step 3: Pr(K = k) against enumerating every partition, each weighted by
# its marginal prior B(a + K, b + T - 1 - K)
for T, l, l1 in [(14, 3, 3), (14, 3, 6), (15, 2, 5)]:
    a_, b_ = 1, T
    found = partitions(T, l, l1)
    weights = [beta(a_ + len(P), b_ + T - 1 - len(P)) for P in found]
    total = sum(weights)
    enumerated = [sum(w for P, w in zip(found, weights) if len(P) == k) / total for k in range(T)]
    check(f"Step 3: Pr(K = k) against enumeration (T={T}, l={l}, first >= {l1}, a=1, b=T)",
          all(simplify(x - y) == 0 for x, y in zip(enumerated, prior_K(T, l, l1, a_, b_))))

# Table tab:changes, for the fit to all reports (T = 157, a = 1, b = T)
table = {  # quantity: (minimum segment, first segment, Pr(K=0..3), Pr(K>=4), mean, largest)
    "R_t": (7, 7, ["0.57", "0.26", "0.11", "0.04"], "0.02", "0.70", 21),
    "CFR_t": (14, 14, ["0.62", "0.26", "0.09", "0.02"], "<0.01", "0.52", 10),
    "theta^C": (14, 69, ["0.78", "0.19", "0.03", "<0.01"], "<0.01", "0.26", 6),
    "theta^D": (14, 82, ["0.82", "0.16", "0.02", "<0.01"], "<0.01", "0.21", 5),
}


def shown(x):
    """A probability as the table shows it: two decimal places, or <0.01."""
    x = float(x)
    return "<0.01" if x < 0.005 else f"{x:.2f}"


T = 157
for name, (l, l1, p0to3, p4, mean, largest) in table.items():
    pr = prior_K(T, l, l1, 1, T)
    computed = ([shown(x) for x in pr[:4]], shown(sum(pr[4:])),
                f"{float(sum(k * x for k, x in enumerate(pr))):.2f}",
                max(k for k, x in enumerate(pr) if x != 0))
    check(f"Table tab:changes, row {name}: {computed}", computed == (p0to3, p4, mean, largest))

finish()

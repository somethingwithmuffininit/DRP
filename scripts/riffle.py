"""Exact total variation of repeated Gilbert–Shannon–Reeds riffles from an ordered deck."""
from fractions import Fraction
from math import comb, factorial


def rising_counts(n):
    """Eulerian counts indexed by descents (one fewer than rising sequences)."""
    if n < 1:
        raise ValueError("A deck needs at least one card")
    counts = [1]
    for size in range(2, n + 1):
        counts = [(d + 1) * right + (size - d) * left
                  for d, (left, right) in enumerate(zip([0] + counts, counts + [0]))]
    return counts


def riffle_tv(n, shuffles):
    if shuffles < 0:
        raise ValueError("The shuffle count must be nonnegative")
    counts = rising_counts(n)
    a = 2 ** shuffles
    orders, labels = factorial(n), a ** n
    deviation = 0
    for descents, count in enumerate(counts):
        ways = comb(a + n - 1 - descents, n)
        deviation += count * abs(orders * ways - labels)
    return Fraction(deviation, 2 * labels * orders)

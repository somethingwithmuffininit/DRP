"""Independent checks with explicit inverse shuffles on permutation groups."""
from collections import Counter
from fractions import Fraction
from itertools import product, permutations
from math import factorial
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from riffle import riffle_tv, rising_counts


class RiffleChecks(unittest.TestCase):
    def test_class_sizes(self):
        for n in range(1, 8):
            actual = Counter(sum(a > b for a, b in zip(p, p[1:])) for p in permutations(range(n)))
            self.assertEqual(rising_counts(n), [actual[d] for d in range(n)])

    def test_explicit_group_walk(self):
        for n in range(2, 6):
            deck = tuple(range(n))
            orders = list(permutations(deck))
            law = Counter({deck: 1})
            # Stable sort independent fair binary labels; tie order is inherited.
            for t in range(5):
                den = (2 ** n) ** t
                exact = sum(abs(Fraction(law[p], den) - Fraction(1, factorial(n))) for p in orders) / 2
                self.assertEqual(exact, riffle_tv(n, t))
                next_law = Counter()
                for p, weight in law.items():
                    for bits in product(range(2), repeat=n):
                        q = tuple(p[i] for i in sorted(range(n), key=bits.__getitem__))
                        next_law[q] += weight
                law = next_law

    def test_published_table_and_boundaries(self):
        self.assertEqual([round(float(riffle_tv(52, k)), 3) for k in (5, 6, 7, 8, 9, 10)],
                         [.924, .614, .334, .167, .085, .043])
        self.assertEqual(riffle_tv(1, 100), 0)
        for args in ((0, 1), (6, -1)):
            with self.assertRaises(ValueError):
                riffle_tv(*args)


if __name__ == "__main__":
    unittest.main()

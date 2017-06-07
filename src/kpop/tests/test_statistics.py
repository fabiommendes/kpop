import pytest
from kpop.statistics import biallelic_pairwise_fst, biallelic_fst


class TestBiallelicFst:
    def test_fst_fixed_opposite_populations(self):
        assert biallelic_fst([1], [0], [0], [0], [10, 10]) == 1

    def test_fst_equal_populations(self):
        e = 1e-6
        assert abs(biallelic_fst([.25], [.25], [.5], [.5], [1e6, 1e6])) < e
        assert biallelic_fst([1], [1], [0], [0], [10, 10]) == 0

    def test_fst_simple_cases(self):
        fst = biallelic_fst

        p1 = 0.5
        h1 = 2 * p1 * (1 - p1)
        p2 = 0.25
        h2 = 2 * p2 * (1 - p2)

        # Identical populations
        assert fst([p1, p2], [p1, p2], [h1, h2], [h1, h2], [500, 500]) <= 0.0
        assert fst([p1, p2], [p1, p2], [h1, h2], [h1, h2], [50, 50]) <= 0.0

        # Different populations
        assert fst([p1, p2], [p2, p1], [h1, h2 / 2], [h1 / 2, h2],
                   [500, 500]) > 0

        # Very different populations (confirm: Fst --> 0.666?)
        assert fst([1., 0.], [0., 1.], [0., 0.], [0., 0.], [50, 50]) > 0.6

    def test_fst_matrix(self):
        data = biallelic_pairwise_fst(
            [[0.1], [0.2], [0.9]],
            [[0.5], [0.4], [0.1]],
            [100, 100, 100],
        )
        result = [
            [0., 0.03749029, 0.78011826],
            [0.03749029, 0., 0.66045591],
            [0.78011826, 0.66045591, 0.]]

        assert abs(data - result).sum() <= 1e-6


def llist(L):
    return list(map(list, L))

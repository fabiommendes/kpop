import pytest
from numpy.testing import assert_almost_equal

from kpop.utils.frequencies import *


@pytest.fixture
def f_flat():
    return [0.1, 0.9, 0.5]


@pytest.fixture
def f_filled():
    return [[0.1, 0.9], [0.9, 0.1], [0.5, 0.5]]


def test_fill_frequencies(f_flat, f_filled):
    assert_almost_equal(fill_freqs_vector(f_filled), f_filled)
    assert_almost_equal(fill_freqs_vector(f_flat), f_filled)


def test_flatten_frequencies(f_flat, f_filled):
    assert (flatten_frequencies(f_flat) == f_flat).all()
    assert (flatten_frequencies(f_filled) == f_flat).all()


def test_freqs_binomial():
    pop = [
        [[1, 2], [1, 1], [2, 2]],
        [[2, 1], [1, 1], [2, 2]],
    ]

    freqs = freqs_binomial(pop, flat=True)
    assert (freqs == [0.5, 1.0, 0.0]).all()

    freqs = freqs_binomial(pop, flat=True, alpha=1.0)
    assert (freqs == [0.5, 5 / 6, 1 / 6]).all()

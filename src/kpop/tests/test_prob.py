import pytest
from kpop.prob import Prob


def test_mode():
    prob = Prob({'a': 0.2, 'b': 0.3, 'c': 0.5})
    assert prob.max() == 0.5
    assert prob.mode() == 'c'


def test_mode_set():
    prob = Prob({'a': 0.2, 'b': 0.3, 'c': 0.3, 'd': 0.2})
    assert prob.max() == 0.3
    assert prob.mode_set() == {'b', 'c'}

import numpy as np
import pytest

from kpop.prob import Prob


class TestProb:
    @pytest.fixture
    def p1(self):
        return Prob({'a': 0.5, 'b': 0.5})

    @pytest.fixture
    def p2(self):
        return Prob({'a': 1, 'b': 0})

    @pytest.fixture
    def p3(self):
        return Prob({'a': 0.75, 'b': 0.25})

    def test_prob_from_list(self):
        p = Prob([0.5, 0.25, 0.25])
        assert p[0] == 0.5
        assert p[1] == 0.25
        assert p[2] == 0.25

    def test_probability_can_be_non_normalized(self):
        p = Prob([1, 1], normalize=False)
        assert sum(p.values()) == 2

    def test_probability_with_support(self, p1):
        p = Prob(p1, support='abc')
        assert p['c'] == 0.0

    def test_prob_dict_interface(self, p1):
        dic = {'a': 0.5, 'b': 0.5}
        assert len(p1) == 2
        assert dict(p1) == dic
        assert set(p1) == set(dic)
        assert set(p1.items()) == set(dic.items())
        assert set(p1.keys()) == set(dic.keys())
        assert set(p1.values()) == set(dic.values())

    def test_update_probability_support(self, p1):
        p1.update_support('cde')
        assert len(p1) == 5
        assert p1['c'] == 0

    def test_set_probability_support(self, p1):
        p1.set_support('abc')
        assert len(p1) == 3

        with pytest.raises(ValueError):
            p1.set_support('c')

    def test_entropy(self, p1, p2):
        assert p1.entropy() == np.log(2)
        assert p1.entropy() > p2.entropy()

    def test_random_elements(self, p1, p2):
        assert all(p2.random_sequence(2) == ['a', 'a'] for _ in range(5))
        assert set(p1.random_sequence(15)) == {'a', 'b'}

    def test_mode(self):
        prob = Prob({'a': 0.2, 'b': 0.3, 'c': 0.5})
        assert prob.max() == 0.5
        assert prob.mode() == 'c'

    def test_mode_set(self, p1, p2):
        assert p1.max() == 0.5
        assert p1.mode_set() == {'a', 'b'}

        assert p2.max() == 1.0
        assert p2.mode_set() == {'a'}

    def test_sharp_prob(self, p1, p2, p3):
        assert p3.sharp() == p2
        assert p1.sharp() == p1

    def test_kl_divergence(self, p1, p2):
        assert p1.kl_divergence(p1) == 0
        assert p1.kl_divergence(p2) == float('inf')
        assert p2.kl_divergence(p2) == 0
        assert abs(p2.kl_divergence(p1) - 0.69) < 5e-2
        assert p1.kl_divergence(Prob([1])) == float('inf')

    def test_encode(self, p3):
        assert list(p3.encode('abc')) == [0.75, 0.25, 0.0]
        assert list(p3.encode()) == [0.75, 0.25]
        assert list(Prob([3, 1]).encode()) == [0.75, 0.25]

    def test_mixture(self, p1, p2):
        mix1 = Prob.mixture([1, 0], [p1, p2])
        mix2 = Prob.mixture([0, 1], [p1, p2])
        mix3 = Prob.mixture([1, 1], [p1, p2])
        assert mix1 == p1
        assert mix2 == p2
        assert mix3 == Prob({'a': 0.75, 'b': 0.25})

    def test_mixture_with_bad_coeffs(self, p1, p2):
        with pytest.raises(ValueError):
            Prob.mixture([1, 1, 1], [p1, p2])

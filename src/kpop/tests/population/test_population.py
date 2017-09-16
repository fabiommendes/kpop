import numpy as np
from numpy.testing import assert_almost_equal

from kpop import Population, Prob
from kpop.population.population_base import sorted_allele_mapping, \
    biallelic_mapping


class TestAuxiliaryFunctions:
    "Private functions used in base population class"

    def test_sorted_allele_mapping(self):
        prob = Prob({1: 0.6, 2: 0.4})
        assert sorted_allele_mapping(prob) == {}

        prob = Prob({2: 0.6, 1: 0.4})
        assert sorted_allele_mapping(prob) == {2: 1, 1: 2}

        prob = Prob({2: 0.5, 3: 0.4, 1: 0.1})
        assert sorted_allele_mapping(prob) == {2: 1, 3: 2, 1: 3}

    def test_biallelic_mapping(self):
        prob = Prob({1: 0.6, 2: 0.4})
        assert biallelic_mapping(prob) == {}

        prob = Prob({2: 0.6, 1: 0.4})
        assert biallelic_mapping(prob) == {}

        prob = Prob({2: 0.6, 1: 0.2, 3: 0.2})
        assert biallelic_mapping(prob) == {2: 1, 1: 2, 3: 2}


class TestProperties:
    "Test basic attributes that expose kpop population data"

    def test_population_basic_attributes(self, popA):
        assert popA.is_biallelic
        assert popA.num_alleles == 2
        assert popA.ploidy == 2
        assert popA.num_loci == 5
        assert popA.size == 8
        assert popA.id == 'A'
        assert popA.populations == [popA]

    def test_population_frequency_attributes(self, popA):
        assert_almost_equal(popA.freqs_vector, [1.0, 0.0, 0.75, 0.25, 0.5])
        assert_almost_equal(popA.freqs_matrix, [
            (1.0, 0), (0.0, 1), (0.75, 0.25), (0.25, 0.75), (0.5, 0.5)
        ])
        assert popA.freqs == [
            {1: 1.00, 2: 0.00},
            {1: 0.00, 2: 1.00},
            {1: 0.75, 2: 0.25},
            {1: 0.25, 2: 0.75},
            {1: 0.50, 2: 0.50},
        ]
        assert list(popA.freqs_vector) == [1.0, 0.0, 0.75, 0.25, 0.5]
        assert list(popA.freqs_matrix[:, 0]) == [1.0, 0.0, 0.75, 0.25, 0.5]
        assert list(popA.hfreqs_vector) == [0, 0, 0.5, 0.5, 1]


class TestAsArray:
    "Test basic slicing and data transformation interfaces"

    def test_conversion_to_array(self, popB):
        arr = popB.as_array('raw')
        assert arr.shape == (4, 5, 2)
        assert isinstance(arr, np.ndarray)
        assert (arr == [
            [[2, 2], [2, 2], [2, 2], [1, 2], [2, 2]],
            [[2, 2], [1, 2], [1, 1], [2, 2], [2, 1]],
            [[2, 2], [1, 2], [1, 1], [1, 2], [2, 1]],
            [[2, 2], [2, 2], [2, 2], [2, 1], [1, 2]],
        ]).all()

    def test_conversion_to_raw_array(self, popB):
        assert (popB.as_array() == popB.as_array('raw')).all()
        assert popB.as_array('raw-unity').shape == (4, 5, 2)

    def test_flat_conversion(self, popB):
        assert popB.as_array('flat').shape == (4, 10)
        assert popB.as_array('rflat').shape == (4, 10)

        flat = popB.as_array('flat-unity')
        assert_almost_equal(flat.mean(0), 0)

        flat = popB.as_array('rflat-unity')
        assert_almost_equal(flat.mean(0), 0)

    def test_conversion_to_count_array(self, popB):
        arr = popB.as_array('count')
        assert arr.shape == (4, 5)
        assert (arr == [
            [0, 0, 0, 1, 0],
            [0, 1, 2, 0, 1],
            [0, 1, 2, 1, 1],
            [0, 0, 0, 1, 1],
        ]).all()


class _TestTransformations:
    def test_drop_non_biallelic(self):
        pop = Population([
            [[1, 3], [1, 2], [0, 1]],
            [[1, 1], [1, 1], [1, 1]],
        ])
        new = pop.drop_non_biallelic()
        assert list(new.as_array()) == [
            [[1, 2], [0, 1]],
            [[1, 1], [1, 1]],
        ]

    def test_force_biallelic(self):
        pop = Population([
            [[1, 3], [1, 2], [0, 1]],
            [[2, 2], [1, 1], [1, 1]],
            [[1, 2], [1, 1], [1, 1]],
        ])
        new = pop.force_biallelic()
        assert list(new.as_array()) == [
            [[2, 2], [1, 2], [0, 1]],
            [[1, 1], [1, 1], [1, 1]],
            [[2, 1], [1, 1], [1, 1]],
        ]

    def test_sort_by_allele_freqs(self):
        pop = Population([
            [[1, 3], [1, 2], [0, 1]],
            [[2, 2], [1, 1], [0, 0]],
            [[1, 2], [1, 1], [1, 0]],
        ])
        new = pop.sort_by_allele_freq()
        assert list(new.as_array()) == [
            [[2, 3], [1, 2], [0, 1]],
            [[1, 1], [1, 1], [0, 0]],
            [[2, 1], [1, 1], [1, 0]],
        ]


class TestRender:
    "Test functions that render a population to string"

    def test_render_population(self, popA):
        assert popA.io.render() == '''
A1: 11 22 12 12 12
A2: 11 22 11 22 21
A3: 11 22 11 22 21
A4: 11 22 11 21 12
A5: 11 22 21 21 21
A6: 11 22 21 22 21
A7: 11 22 11 12 12
A8: 11 22 21 22 12
    '''.strip()

    def test_render_population_biallelic_freqs(self, popA, popB):
        freqs_render = (popA + popB).stats.render_biallelic_freqs()
        assert freqs_render == '''
locus         A         B
   L1  1.000000  0.000000
   L2  0.000000  0.250000
   L3  0.750000  0.500000
   L4  0.250000  0.375000
   L5  0.500000  0.375000
'''[1:]

        assert (popA + popB).stats.render_biallelic_freqs(sep=', ') == '''
locus,        A,        B
   L1, 1.000000, 0.000000
   L2, 0.000000, 0.250000
   L3, 0.750000, 0.500000
   L4, 0.250000, 0.375000
   L5, 0.500000, 0.375000
'''[1:]

    def _test_render_population_frequencies(self, popA, popB):
        assert (popA + popB).render_biallelic_frequencies(sep=', ') == '''

'''.strip()


class TestRandomPopulations:
    "Test the creation of random populations"

    def test_make_random_population(self):
        pop = Population.random(30, 20)
        assert pop.size == 30
        assert pop.num_loci == 20
        delta = pop.freqs_vector - \
                pop.stats.empirical_freqs(as_matrix=True)[:, 0]
        assert abs(delta).mean() < 0.15

    def test_make_random_population_with_seed(self):
        popA = Population.random(5, 10, seed=0)
        popB = Population.random(5, 10, seed=0)
        assert (popA.freqs_vector == popB.freqs_vector).all()
        assert popA == popB

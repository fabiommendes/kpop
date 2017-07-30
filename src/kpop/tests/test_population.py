import pytest
import numpy as np
from numpy.testing import assert_almost_equal

from kpop import Population, Individual


#
# Test rendering and attributes
#
class TestProperties:
    "Test basic attributes that expose kpop population data"

    def test_population_basic_attributes(self, popA):
        assert popA.is_biallelic
        assert popA.num_alleles == 2
        assert popA.ploidy == 2
        assert popA.num_loci == 5
        assert popA.size == 8
        assert popA.label == 'A'
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

    def test_fst_statistics(self, popA, popB):
        popC = Population(popA, label='C', copy=True)
        pop = popA + popB + popC
        assert pop.render_fst() == (
            '           A         B\n'
            'B:  0.409558\n'
            'C: -0.028571  0.409558'
        )


class TestRender:
    "Test functions that render a population to string"

    def test_render_population(self, popA):
        assert popA.render() == '''
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
        freqs_render = (popA + popB).render_biallelic_freqs()
        assert freqs_render == '''
locus         A         B
   L1  1.000000  0.000000
   L2  0.000000  0.250000
   L3  0.750000  0.500000
   L4  0.250000  0.375000
   L5  0.500000  0.375000
'''[1:]

        assert (popA + popB).render_biallelic_freqs(sep=', ') == '''
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


#
# Test evolution and creation of new offspring
#
def test_create_individual_labels(popA):
    assert popA[0].label == 'A1'

def test_create_offspring(popA):
    ind3 = popA.random_individual()
    ind2 = popA.new_offspring(label='Achild')
    ind1 = popA.new(label='Arandom')

    for ind in ind1, ind2, ind3:
        print(ind)
        assert ind.population is popA
        assert ind.label is not None
        assert ind[0, 0] == 1
        assert ind[1, 0] == 2


def test_fill_population_uses_the_correct_frequencies(popA):
    popA.fill(10)

    # These are fixed alleles
    assert popA[9][0, 0] == 1
    assert popA[9][0, 1] == 1
    assert popA[9][1, 0] == 2
    assert popA[9][1, 1] == 2


def test_population_genetic_drift(popA):
    size = 1e6  # we need a very large size to avoid going back to the same
    # frequence by pure chance

    # Frequencies do not change with evolution of zero generations
    pop2 = popA.genetic_drift(0, sample_size=0, population_size=size)
    assert popA.freqs == pop2.freqs

    # Genetic drift
    pop2 = popA.genetic_drift(10, sample_size=0, population_size=size)
    assert popA.freqs != pop2.freqs

    # Fixed alleles do not change
    assert popA.freqs[0] == pop2.freqs[0]
    assert popA.freqs[1] == pop2.freqs[1]

    # Non-fixed alleles should always change
    for i in range(2, 5):
        assert popA.freqs[i] != pop2.freqs[i]


def test_make_random_population():
    pop = Population.make_random(30, 20)
    assert pop.size == 30
    assert pop.num_loci == 20
    delta = pop.freqs_vector - pop.empirical_freqs(as_matrix=True)[:, 0]
    assert abs(delta).mean() < 0.15


def test_add_remove_individual(popA):
    popA.add(popA.new_offspring())
    popA.remove()

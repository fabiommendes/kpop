import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from kpop import Population, Prob, MultiPopulation
from kpop.population.population_base import sorted_allele_mapping, \
    biallelic_mapping


def deeplist(array):
    """
    Convert array in a list of lists.
    """

    def conv(arr, depth):
        if depth <= 1:
            return list(arr)
        return [conv(line, depth - 1) for line in arr]

    return conv(array, array.ndim)


class TestPopulation:
    """
    Test PopulationBase interface.
    """

    #
    # Test basic attributes that expose kpop population data
    #
    def test_population_basic_attributes(self, popA):
        assert popA.is_biallelic
        assert popA.num_alleles == 2
        assert popA.ploidy == 2
        assert popA.num_loci == 5
        assert popA.size == 8
        assert popA.id == 'A'
        assert popA.populations == [popA]
        assert popA.num_populations == 1
        assert popA.has_missing_data is False
        assert popA.missing_data_ratio == 0
        assert popA.missing_data_total == 0.0

    def test_population_special_attributes_aliases(self, popA):
        assert popA.admixture is popA.admix
        assert popA.classification is popA.cls
        assert popA.clusterization is popA.cluster
        assert popA.projection is popA.proj
        assert popA.simulation is popA.sim
        assert popA.statistics is popA.stats

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

    #
    # Basic indexing, slicing and data extraction interfaces
    #
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

    def test_accept_fetch_individual_by_label(self, popB):
        assert popB['B1'] == popB[0]

    def test_conversion_to_raw_array(self, popB):
        assert (popB.as_array() == popB.as_array('raw')).all()
        assert popB.as_array('raw-norm' ).shape == (4, 5, 2)

    def test_flat_conversion(self, popB):
        assert popB.as_array('flat').shape == (4, 10)
        assert popB.as_array('rflat').shape == (4, 10)

        flat = popB.as_array('flat-norm' )
        assert_almost_equal(flat.mean(0), 0)

        flat = popB.as_array('rflat-norm' )
        assert_almost_equal(flat.mean(0), 0)

    def test_conversion_to_count_array(self, popB):
        arr = popB.as_array('count')
        assert arr.shape == (4, 5)
        assert deeplist(arr) == [[0, 0, 0, 1, 0],
                                 [0, 1, 2, 0, 1],
                                 [0, 1, 2, 1, 1],
                                 [0, 0, 0, 1, 1]]

    def test_conversion_to_count_snp(self, popB):
        arr = popB.as_array('count-snp')
        assert arr.shape == (4, 5)
        assert arr.dtype == np.dtype(float)

    def test_conversion_to_count_center(self, popB):
        arr = popB.as_array('count-center')
        assert arr.shape == (4, 5)
        assert deeplist(arr) == [[-1., -1, -1, 0, -1],
                                 [-1., 0, 1, -1, 0],
                                 [-1., 0, 1, 0, 0],
                                 [-1., -1, -1, 0, 0]]

    def test_invalid_conversion_method(self, popA):
        with pytest.raises(ValueError):
            popA.as_array('bad-method')

    def test_individual_id_labels(self, popB):
        assert popB.individual_ids == ['B1', 'B2', 'B3', 'B4']
        assert popB[0].id == 'B1'

    #
    # Private functions used in base population class
    #
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

    #
    # Transformations
    #
    def test_drop_individuals(self, popA):
        new = popA.drop_individuals(range(popA.size - 2))
        assert new.size == 2
        assert new[1] == popA[-1]

    def test_drop_non_biallelic(self):
        pop = Population([
            [[1, 3], [1, 2], [0, 1]],
            [[1, 1], [1, 1], [1, 1]],
        ])
        new, indexes = pop.drop_non_biallelic()
        assert deeplist(new.as_array()) == [
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
        assert deeplist(new.as_array()) == [
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
        assert deeplist(new.as_array()) == [
            [[2, 3], [1, 2], [0, 1]],
            [[1, 1], [1, 1], [0, 0]],
            [[2, 1], [1, 1], [1, 0]],
        ]

    def test_keep_loci(self, popA):
        pop = popA.keep_loci([0])
        assert pop.as_array().shape == (popA.size, 1, popA.ploidy)

    def test_map_alleles(self, popA):
        maps = [{} for _ in range(popA.num_loci)]
        assert popA == popA.map_alleles(maps)


class TestMultiPopulationInteface(TestPopulation):
    @pytest.fixture
    def popA(self):
        from kpop.tests.conftest import popA, popA_data

        pop = popA(popA_data())
        return MultiPopulation([pop], id='A')


class TestEmptyPopulations:
    def test_cannot_initialize_with_empty_frequency_list(self):
        with pytest.raises(ValueError):
            Population(freqs=[])
        with pytest.raises(ValueError):
            Population(freqs=[[[]]])

    def test_initialize_an_empty_population(self):
        pop1 = Population(freqs=[{1: 0.5, 2: 0.5}, {1: 1.0, 2: 0.0}])
        pop2 = Population(freqs=[0.5, 1.0])
        pop3 = Population(freqs=[[0.5, 0.5], [1.0, 0.0]])
        assert pop1.freqs == pop2.freqs
        assert pop2.freqs == pop3.freqs

    def test_set_freqs_attribute(self, popB):
        with pytest.raises(ValueError):
            popB.freqs = [0.5, 0.5]

        popB.freqs = [Prob([0.5, 0.5]) for _ in range(popB.num_loci)]


class TestSinglePopulations:
    """
    Test specific funtions for the PopulationSingle class"
    """

    def test_random_invalid_args(self):
        with pytest.raises(ValueError):
            Population.random(10, num_loci=-1)

    def test_can_init_from_other_population(self, popA):
        pop = Population(popA)
        assert pop == popA
        assert pop.individual_ids == popA.individual_ids

    def test_population_equality_with_multi_population(self, popA):
        assert popA == MultiPopulation([popA])

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

    def test_copy_with_bad_data(self, popA):
        with pytest.raises(ValueError):
            popA.copy(data=[0, 1])

    def test_functions_with_bad_data_argument(self, popA):
        data = popA.as_array()

        with pytest.raises(TypeError):
            maps = [{} for _ in range(popA.num_loci)]
            popA.map_alleles(maps, data=data)

        with pytest.raises(TypeError):
            popA.keep_loci([0, 1], data=data)

        with pytest.raises(TypeError):
            popA.keep_individuals([0, 1], data=data)


class TestMultiPopulation:
    """
    Test specific interfaces in MultiPopulation.
    """

    def test_multi_population(self, popA, popB):
        pop = popA + popB
        assert pop[0] == popA[0]
        assert pop[len(popA)] == popB[0]
        assert len(pop) == len(popA) + len(popB)
        assert len(list(pop)) == len(popA) + len(popB)

    def test_population_attribute_behaves_as_a_list(self, popA, popB):
        pop = popA + popB
        del pop.populations[1]
        assert pop == popA

    def test_multipopulation_with_3_populations(self, popA, popB):
        popC = Population(popB[:-1], id='C')
        pop = popA + popB
        assert type(pop) == MultiPopulation
        assert len(pop.populations) == 2
        assert pop.populations == [popA, popB]

        pop = popA + popB + popC
        assert type(pop) == MultiPopulation
        assert len(pop.populations) == 3
        assert pop.populations == [popA, popB, popC]

        pop = popA + (popB + popC)
        assert type(pop) == MultiPopulation
        assert len(pop.populations) == 3
        assert pop.populations == [popA, popB, popC]

    def test_prevent_repeated_populations(self, popA):
        with pytest.raises(ValueError):
            MultiPopulation([popA, popA])

    def test_filter_individuals(self, popA, popB):
        pop = popA + popB
        idx = pop.slice_indexes([0, len(popA)])
        assert deeplist(idx[0]) == [0]
        assert deeplist(idx[1]) == [0]

        new = pop.keep_individuals([0, len(popA)])
        assert new.size == 2
        assert new[0] == popA[0]
        assert new[1] == popB[0]

    def test_random_creates_multipopulation_when_given_a_list_of_sizes(self):
        pops = Population.random([1, 2, 3], 5)
        assert isinstance(pops, MultiPopulation)
        assert len(pops.populations) == 3

        p1, p2, p3 = pops.populations
        assert (p1.freqs_vector != p2.freqs_vector).all()

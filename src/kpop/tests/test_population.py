import pytest
from numpy.testing import assert_almost_equal

from kpop import Population, Individual


@pytest.fixture
def individuals_data():
    return [
        '11 22 12 12 12',
        '11 22 11 22 21',
        '11 22 11 22 21',
        '11 22 11 21 12',
        '11 22 21 21 21',
        '11 22 21 22 21',
        '11 22 11 12 12',
        '11 22 21 22 12',
    ]


@pytest.fixture
def pop(individuals_data):
    return Population([Individual(x) for x in individuals_data], label='A')


def test_population_attributes(pop):
    assert pop.is_biallelic
    assert pop.num_alleles == 2
    assert pop.ploidy == 2
    assert pop.num_loci == 5
    assert pop.size == 8
    assert pop.label == 'A'


def test_population_frequencies(pop):
    assert_almost_equal(pop.freqs_vector, [1.0, 0.0, 0.75, 0.25, 0.5])
    assert_almost_equal(pop.freqs_matrix, [
        (1.0, 0), (0.0, 1), (0.75, 0.25), (0.25, 0.75), (0.5, 0.5)
    ])
    assert pop.freqs == [
        {1: 1.00, 2: 0.00},
        {1: 0.00, 2: 1.00},
        {1: 0.75, 2: 0.25},
        {1: 0.25, 2: 0.75},
        {1: 0.50, 2: 0.50}
    ]


def test_create_individual_labels(pop):
    assert pop[0].label == 'A1'


def test_render_population(pop):
    assert pop.render() == '''
A1: 11 22 12 12 12
A2: 11 22 11 22 21
A3: 11 22 11 22 21
A4: 11 22 11 21 12
A5: 11 22 21 21 21
A6: 11 22 21 22 21
A7: 11 22 11 12 12
A8: 11 22 21 22 12
'''.strip()


def test_create_offspring(pop):
    ind3 = pop.random()
    ind2 = pop.new_offspring(label='Achild')
    ind1 = pop.new(label='Arandom')

    for ind in ind1, ind2, ind3:
        print(ind)
        assert ind.population is pop
        assert ind.label is not None
        assert ind[0, 0] == 1
        assert ind[1, 0] == 2


def test_fill_population(pop):
    pop.fill(10)
    assert pop[9][0, 0] == 1
    assert pop[9][0, 1] == 1
    assert pop[9][1, 0] == 2
    assert pop[9][1, 1] == 2


def test_population_genetic_drift(pop):
    # Frequencies do not change with evolution of zero generations
    pop2 = pop.genetic_drift(0, sample_size=0)
    assert pop.freqs == pop2.freqs

    # Fixed alleles do not change
    pop2 = pop.genetic_drift(10, sample_size=0)
    assert pop.freqs[0] == pop2.freqs[0]
    assert pop.freqs[1] == pop2.freqs[1]

    # Non-fixed alleles should always change
    assert pop.freqs[2] != pop2.freqs[2]
    assert pop.freqs[3] != pop2.freqs[3]
    assert pop.freqs[4] != pop2.freqs[4]


def test_can_analyze_pca(pop):
    data = pop.pca(k=3, method='count')
    assert data.shape == (8, 3)

    data = pop.pca(k=3, method='flatten')
    assert data.shape == (8, 3)


def test_make_random_population():
    pop = Population.make_random(30, 20)
    assert pop.size == 30
    assert pop.num_loci == 20
    delta = pop.freqs_vector - pop.empirical_freqs(as_matrix=True)[:, 0]
    assert abs(delta).mean() < 0.15


def test_add_remove_individual(pop):
    pop.add(pop.new_offspring())
    pop.remove()


def test_populations_list(pop):
    assert list(pop.populations) == [pop]

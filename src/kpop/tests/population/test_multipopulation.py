import pytest

from kpop import Population, MultiPopulation


@pytest.fixture
def numloci():
    return 50  # that should make the classification more efficient


class TestMultiPopulation:
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

    def test_multipopulation_with_3_populatoins(self, popA, popB):
        popC = Population(popB[:-1], id='C', copy=True)
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


#
# Clusterization and classification
#
def test_sharp_classifier(popA, popB):
    pop = popA + popB
    for ind in popA:
        assert pop.classify(ind) == 'A'
    for ind in popB:
        assert pop.classify(ind) == 'B'


def test_probabilistic_classifier(popA, popB):
    pop = popA + popB
    pAs = []

    for ind in popA:
        prob = pop.prob_classify(ind)
        pA = prob['A']
        pAs.append(pA)

    print(pAs)
    assert all([p > 0.75 for p in pAs])

import pytest


@pytest.fixture
def numloci():
    return 50  # that should make the classification more efficient


def test_multi_population(popA, popB):
    pop = popA + popB
    assert pop[0] == popA[0]
    assert pop[5] == popB[0]
    assert len(pop) == len(popA) + len(popB)
    assert len(list(pop)) == len(popA) + len(popB)


def test_classification(popA, popB, pop):
    for ind in popA:
        assert pop.classify(ind) == 'A'
    for ind in popB:
        assert pop.classify(ind) == 'B'


def test_admixture(popA, pop):
    pAs = []
    for ind in popA:
        admx = pop.prob_classify(ind)
        pA = admx['A']
        pAs.append(pA)

    print(pAs)
    assert all([p > 0.75 for p in pAs])

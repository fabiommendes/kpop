import pytest

from kpop import Population


@pytest.fixture
def numloci():
    return 20


@pytest.fixture
def popA(numloci):
    return Population.make_random(5, numloci, label='A', min_prob=0.1)


@pytest.fixture
def popB(numloci):
    return Population.make_random(5, numloci, label='B', min_prob=0.9)


@pytest.fixture
def pop(popA, popB):
    return popA + popB

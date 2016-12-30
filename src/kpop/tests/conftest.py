import pytest

from kpop.individual import Individual


@pytest.fixture
def popA():
    return [Individual.random_biallelic([0.1] * 20, label='A%s' % (i + 1))
            for i in range(10)]


@pytest.fixture
def popB():
    return [Individual.random_biallelic([0.9] * 20, label='A%s' % (i + 1))
            for i in range(10)]


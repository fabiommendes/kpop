import matplotlib as mpl
import pytest

mpl.use('Agg')

from kpop import Population, Individual
from tests import load_data, data_path, temporary_location


#
# Example population
#
@pytest.fixture
def popA_data():
    return popA_data_()


def popA_data_():
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
def popB_data():
    return [
        '22 22 22 12 22',
        '22 12 11 22 21',
        '22 12 11 12 21',
        '22 22 22 21 12',
    ]


@pytest.fixture
def popA(popA_data):
    return popA_(popA_data)


def popA_(data):
    return Population([Individual(x) for x in data], id='A')


@pytest.fixture
def popB(popB_data):
    return Population([Individual(x) for x in popB_data], id='B')


#
# Random populations
#
@pytest.fixture
def num_loci():
    return 20


@pytest.fixture
def popA_random(num_loci):
    return Population.random(5, num_loci, id='A', min_prob=0.1)


@pytest.fixture
def popB_random(num_loci):
    return Population.random(5, num_loci, id='B', min_prob=0.9)


@pytest.fixture
def popAB(popA, popB):
    return popA + popB


@pytest.fixture(scope='session')
def loader():
    return load_data


@pytest.fixture(scope='session')
def path():
    return data_path


@pytest.fixture(scope='session')
def temp_path():
    return temporary_location

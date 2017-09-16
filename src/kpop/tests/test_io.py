import os
import pytest

from kpop import Population, load_csv


def load(f):
    return os.path.join(os.path.dirname(__file__), 'data', f)


class TestIoLoaders:
    @pytest.fixture
    def pop(self):
        return Population('11 11 12\n11 22 11')

    def test_load_csv_file(self, pop):
        pop_csv = load_csv(open(load('pop-small.csv')))
        assert pop_csv == pop
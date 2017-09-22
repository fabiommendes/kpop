import pytest

from kpop import Population, load_csv


class TestIoLoaders:
    @pytest.fixture
    def pop(self):
        return Population('11 11 12\n11 22 11')

    def test_load_csv_file(self, pop, loader):
        pop_csv = load_csv(loader('pop-small.csv'))
        assert pop_csv == pop

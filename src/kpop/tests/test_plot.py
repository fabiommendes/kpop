import pytest
from kpop import Population

class TestPlots:
    @pytest.fixture
    def pop(self):
        return Population('112211\n121212')

    def test_pca_plot(self, pop):
        pop.plot.pca()
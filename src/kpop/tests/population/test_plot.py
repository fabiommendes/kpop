import pytest
from kpop import Population

class TestPlots:
    @pytest.fixture
    def pop(self):
        return Population('11 22 11\n12 12 12', ploidy=2)

    def test_pca_plot(self, pop):
        assert pop.shape == (2, 3, 2)
        pop.plot.pca()
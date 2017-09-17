import pytest

from kpop import Population


class TestProjection:
    def test_can_analyze_pca(self, popA):
        data = popA.projection.pca(k=3, data='count')
        assert data.shape == (8, 3)

        data = popA.projection.pca(k=3, data='flat')
        assert data.shape == (8, 3)

    def test_pca_do_not_accept_pca_arg(self, popB):
        with pytest.raises(TypeError):
            popB.projection.pca(pca=5)

    def test_transformation_with_reduction(self, popB):
        data = popB.projection.ica(pca=3)
        assert data.shape == (4, 2)

    def test_kernel_pca(self, popB):
        data = popB.projection.kernel_pca()
        assert data.shape == (4, 2)

    def test_ica(self, popB):
        data = popB.projection.ica()
        assert data.shape == (4, 2)

    def test_nmf(self, popB):
        data = popB.projection.nmf()
        assert data.shape == (4, 2)

    def test_mds(self, popB):
        data = popB.projection.mds()
        assert data.shape == (4, 2)

    def test_spectral(self, popB):
        data = popB.projection.spectral()
        assert data.shape == (4, 2)

    def test_tsne(self, popB):
        data = popB.projection.tsne()
        assert data.shape == (4, 2)

    def test_lle(self, popB):
        data = popB.projection.lle()
        assert data.shape == (4, 2)

    def test_isomap(self, popB):
        data = popB.projection.isomap()
        assert data.shape == (4, 2)

    def test_factor(self, popB):
        data = popB.projection.factor()
        assert data.shape == (4, 2)

    def test_lda_projection(self, popB):
        data = popB.projection.lda_projection()
        assert data.shape == (4, 2)

    def test_transformer(self, popB):
        from sklearn.decomposition import PCA

        data = popB.projection(PCA)
        assert data.shape == (4, 2)

    def test_project_alias(self, popB):
        data = popB.projection.project('pca')
        assert data.shape == (4, 2)

    def test_called_with_invalid_method(self, popB):
        with pytest.raises(ValueError):
            popB.projection('wrong-value')

        with pytest.raises(ValueError):
            popB.projection(42)

    def test_invalid_arguments(self, popB):
        with pytest.raises(ValueError):
            popB.projection.lda_projection(k=2, n_populations=2)

        with pytest.raises(TypeError):
            popB.projection.nmf(pca=3)


class TestPlots:
    @pytest.fixture
    def pop(self):
        return Population('11 22 11\n12 12 12', ploidy=2)

    def test_pca_plot(self, pop):
        assert pop.shape == (2, 3, 2)
        pop.plot.pca()

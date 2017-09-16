class TestProjectors:
    def test_can_analyze_pca(self, popA):
        data = popA.projection.pca(k=3, data='count')
        assert data.shape == (8, 3)

        data = popA.projection.pca(k=3, data='flat')
        assert data.shape == (8, 3)

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
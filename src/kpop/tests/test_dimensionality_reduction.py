def test_can_analyze_pca(popA):
    data = popA.projection.pca(k=3, method='count')
    assert data.shape == (8, 3)

    data = popA.projection.pca(k=3, method='flatten')
    assert data.shape == (8, 3)

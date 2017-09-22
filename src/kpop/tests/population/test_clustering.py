class TestClustering:
    """
    Tests the Population.clusterization attribute.
    """


class TestClassification:
    """
    Tests the Population.classification attribute.
    """

    def test_simple_classification(self, popA, popB):
        cls = (popA + popB).classification(['A'] * 8 + ['B'] * 4)
        assert cls.classify(popB) == ['B', 'B', 'B', 'B']

    def test_ancestry_classification(self, popA, popB):
        cls = (popA + popB).classification('ancestry')
        assert cls.classify(popB) == ['B', 'B', 'B', 'B']

    def test_sharp_classifier(self, popA, popB):
        pop = popA + popB
        for ind in popA:
            assert pop.cls.parent(ind) == 'A'
        for ind in popB:
            assert pop.cls.parent(ind) == 'B'

    def test_probabilistic_classifier(self, popA, popB):
        pop = popA + popB
        pAs = []

        for ind in popA:
            prob = pop.cls.parent_prob(ind)
            pA = prob['A']
            pAs.append(pA)

        print(pAs)
        assert all([p > 0.75 for p in pAs])

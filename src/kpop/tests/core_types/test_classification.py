import pytest
from sidekick import L
from sklearn.naive_bayes import MultinomialNB

from kpop.classifiers import SklearnClassifier


class TestClassification:
    """
    Tests the Population.classification attribute.
    """

    @pytest.fixture
    def cls(self, popA, popB):
        return (popA + popB).classification('ancestry')

    def test_simple_classification(self, popA, popB):
        cls = (popA + popB).classification(['A'] * 8 + ['B'] * 4)
        assert list(cls.classify(popB)) == ['B', 'B', 'B', 'B']

    def test_ancestry_classification(self, popA, popB):
        cls = (popA + popB).classification('ancestry')
        assert list(cls.classify(popB)) == ['B', 'B', 'B', 'B']

    def test_random_classification(self, cls, popB):
        # Prob matrix
        prob, = cls.prob_matrix(popB[0])
        assert prob[0] < prob[1]

        # Prob list
        prob, = cls.prob_list(popB[0])
        assert prob['A'] < prob['B']

        # Probability dataframe
        prob = cls.prob_table(popB[0])
        assert prob.loc[0, 'A'] < prob.loc[0, 'B']

    def test_can_be_created_from_classifier(self, popB):
        cls = popB.classification(list('aabb'), MultinomialNB, data='count')
        assert isinstance(cls, SklearnClassifier)

    def test_can_be_created_from_name(self, popB):
        cls = popB.classification(list('aabb'), 'naive-bayes', data='count')
        assert isinstance(cls, SklearnClassifier)

        with pytest.raises(ValueError):
            popB.classification([1, 1, 2, 2], 'bad-method')
        with pytest.raises(ValueError):
            popB.classification([1, 1, 2, 2], ['bad-type'])

    def test_normalize_labels(self, popB):
        norm = popB.classification._normalize_labels
        labels = [1, 1, 2, 2]
        popB.meta['foo'] = labels
        assert norm(labels) == labels
        assert norm('ancestry') == ['B', 'B', 'B', 'B']
        assert (norm('foo') == labels).all()

        with pytest.raises(ValueError):
            norm()
        popB.meta['labels'] = labels
        assert (norm() == labels).all()

        # Bad cases
        with pytest.raises(ValueError):
            norm('bad')
        with pytest.raises(ValueError):
            norm([1, 2])

    def test_naive_bayes_with_bad_method(self, popB):
        labels = list('aabb')
        for method in ['count-unity', 'flat-unity', 'count-snp']:
            with pytest.raises(ValueError):
                popB.classification.naive_bayes(labels, data=method)

    def test_bernoulli_naive_bayes(self, popB):
        # cls = popB.classification.naive_bayes(['a', 'a', 'a', 'b'], data='flat')
        # assert cls.prob_matrix(popB).shape == (4, 2)
        pass

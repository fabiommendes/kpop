from lazyutils import lazy

from .libs import lazy_module
from .libs import np
from .prob import Prob

population = lazy_module('kpop.population')


class Classifier:
    """
    Abstract classifier interface.
    """

    _label_set = lazy(lambda x: sorted(set(x.labels)))
    _labels_map = lazy(lambda x: {x: i for i, x in enumerate(x._label_set)})

    @lazy
    def _labels_encoder(self):
        sorted_items = sorted(self._labels_map.items(), key=lambda _: _[1])
        return [k for k, i in sorted_items]

    @lazy
    def _integer_labels(self):
        map = self._labels_map
        return np.array([map[label] for label in self.labels])

    def __init__(self, training_data, transformer, labels):
        self.training_data = training_data
        self.labels = labels
        self._transformer = transformer

    def __call__(self, pop):
        """
        self(pop) <==> self.classify(pop)
        """
        return self.classify(pop)

    def transform(self, pop):
        """
        Transform an individual or population to a raw np.ndarray data set.
        """
        if getattr(pop, 'is_individual', False):
            pop = population.Population([pop.data])
        return self._transformer(pop)

    def classify(self, pop):
        """
        Classify a population object.
        """
        is_individual = getattr(pop, 'is_individual', False)
        data = self.transform(pop)
        result = self.classify_data(data)
        return result[0] if is_individual else result

    def classify_data(self, data):
        """
        Classify raw data returning a list of labels.
        """
        raise NotImplementedError

    def prob_matrix(self, pop):
        """
        Return a matrix with the probability that each individual is classified
        with each label. Individuals are represented in the rows and labels
        in the columns.

        Label indexes are assigned by ordering, e.g., If the original labels
        contain 'foo', 'bar' and 'baz', 'bar' will be assigned a column index
        of 0  (because it is the first in alphabetical order), 'baz' will be the
        second column and 'foo' the third.
        """

        logp = self.log_prob_matrix(pop)
        logp -= logp.max()
        probs = np.exp(logp)
        probs /= probs.sum()
        return probs

    def prob_list(self, pop):
        """
        Return a list of :cls:`kpop.Prob` objects with the probabilities
        assigned to each label classification.
        """

        values = self._labels_encoder
        return [Prob(zip(values, row)) for row in self.prob_matrix(pop)]

    def prob_table(self, pop):
        """
        Return a pandas dataframe with the probabilities that each individual
        belongs to each label.
        """
        from pandas import DataFrame

        data = self.prob_matrix(pop)
        return DataFrame(data, columns=self._labels_encoder)

    def log_prob_matrix(self, pop):
        """
        Like :meth:`prob_matrix`, but returns the log probabilities.
        """
        if type(self).prob_matrix is not Classifier.prob_matrix:
            return np.log(self.prob_matrix(pop))

        raise NotImplementedError(
            "either 'log_prob_matrix' or 'prob_matrix' must be defined"
        )

    def log_prob_list(self, pop):
        """
        Like :meth:`prob_list`, but returns the log probabilities.
        """
        values = self._labels_encoder
        return [Prob(zip(values, row), normalize=False)
                for row in self.log_prob_matrix(pop)]

    def log_prob_table(self, pop):
        """
        Like :meth:`prob_table`, but returns the log probabilities.
        """
        from pandas import DataFrame

        data = self.log_prob_matrix(pop)
        return DataFrame(data, columns=self._labels_encoder)


class SklearnClassifier(Classifier):
    """
    A Kpop classifier façade for Scikit learn classifiers.
    """

    def __init__(self, classifier, training_data, transformer, labels,
                 **kwargs):
        super().__init__(training_data, transformer, labels)

        # Train the classifier
        self.classifier = classifier(**kwargs)
        self.classifier.fit(self.training_data, self._integer_labels)

    def classify_data(self, data):
        int_labels = self.classifier.predict(data)
        map = self._labels_encoder
        return [map[i] for i in int_labels]

    def log_prob_matrix(self, pop):
        data = self.transform(pop)
        return self.classifier.predict_log_proba(data)

    def prob_matrix(self, pop):
        data = self.transform(pop)
        return self.classifier.predict_proba(data)

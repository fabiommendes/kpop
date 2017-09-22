import numpy as np

from .utils import lazy_module

population = lazy_module('kpop.population')


class Classifier:
    """
    Abstract classifier interface.
    """

    def __init__(self, training_data, transformer, labels):
        self.training_data = training_data
        self.labels = labels
        self._transformer = transformer

    def __call__(self, pop):
        """
        self(pop) <==> self.classify(pop)
        """
        return self.classify(pop)

    def classify(self, pop):
        """
        Classify a population object.
        """
        is_individual = isinstance(pop, population.Individual)
        if is_individual:
            pop = population.Population([pop])
        data = self._transformer(pop)
        result = self.classify_data(data)
        return result[0] if is_individual else result

    def classify_data(self, data):
        """
        Classify raw data returning a list of labels.
        """
        raise NotImplementedError


class SklearnClassifier(Classifier):
    """
    A Kpop classifier façade for Scikit learn classifiers.
    """

    def __init__(self, classifier, training_data, transformer, labels,
                 **kwargs):
        super().__init__(training_data, transformer, labels)

        # Prepare labels
        self._labels_map = map = {x: i for i, x in enumerate(set(self.labels))}
        self._labels_encoder = [
            k for k, i in
            sorted(map.items(), key=lambda x: x[1])
            ]
        self._integer_labels = np.array([map[label] for label in self.labels])

        # Train the classifier
        self.classifier = classifier(**kwargs)
        self.classifier.fit(self.training_data, self._integer_labels)

    def classify_data(self, data):
        int_labels = self.classifier.predict(data)
        map = self._labels_encoder
        return [map[i] for i in int_labels]
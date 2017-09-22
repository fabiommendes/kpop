import numpy as np

from kpop.classifiers import SklearnClassifier
from .attr import Attr
from ..admixture import likelihood
from ..prob import Prob
from ..utils.checks import is_sklearn_classifier


class Classification(Attr):
    """
    Implements the population.classification attribute.
    """

    _methods = ('naive_bayes',)

    def __call__(self, labels=None, which='naive_bayes', **kwargs):
        if is_sklearn_classifier(which):
            return self.sklearn(which, labels, **kwargs)
        elif callable(which):
            raise NotImplementedError('do not accept function classifiers')
        elif isinstance(which, str):
            which_ = which.lower().replace('-', '_')
            if which_ in self._methods:
                method = getattr(self, which_)
                return method(labels, **kwargs)

    def _normalize_labels(self, labels):
        "Normalizes the labels attribute and return a sequence of labels."

        pop = self._population
        if labels is None or labels == '':
            try:
                return pop.meta['labels']
            except KeyError:
                raise ValueError('could not fetch labels from metadata')
        elif isinstance(labels, str):
            if labels in pop.meta:
                return pop.meta[labels]
            if labels == 'ancestry':
                return ancestry_labels(pop)
            raise ValueError('could not find %r metadata' % labels)
        else:
            labels = list(labels)
            if len(labels) != pop.size:
                raise ValueError(
                    'list of labels must have the same size as the population'
                )
            return labels

    def naive_bayes(self, labels=None, data='count', **kwargs):
        """
        Classify objects using the naive_bayes classifier.
        """
        from sklearn import naive_bayes

        if data == 'count':
            classifier = naive_bayes.MultinomialNB
        elif data == 'flat':
            classifier = naive_bayes.BernoulliNB
        else:
            raise ValueError(
                'naive bayes only accets "count" and "flat" for the data '
                'argument'
            )
        return self.sklearn(classifier, labels, data=data, **kwargs)

    def svm(self, labels=None, data='count', **kwargs):
        """
        Classify objects using the Support Vector Machine (SVM) classifier.
        """
        from sklearn.svm import SVC as classifier
        return self.sklearn(classifier, labels, data=data, **kwargs)

    def sklearn(self, classifier, labels=None, data='count', **kwargs):
        """
        Uses a scikit learn classifier to classify population.

        Args:
            classifier:
                A scikit learn classifier class (e.g.,
                sklearn.naive_bayes.BernoulliNB)
            labels:
                A sequence of labels used to train the classifier.
            data (str):
                The method used to convert the population to a usable data set.
                It uses the same options as in the :meth:`Population.as_array`
                method.
        """
        func = lambda pop: pop.as_array(data)
        raw_data = func(self._population)
        labels = self._normalize_labels(labels)
        return SklearnClassifier(classifier, raw_data, func, labels, **kwargs)

    def populations(self, populations=None):
        """
        Normalize the populations argument and return a list of populations.

        If populations is not give or if it is None, return the .populations
        attribute of the current population.
        """

        if populations is None:
            return self._population.populations
        try:
            return populations.populations
        except AttributeError:
            return list(populations)

    def prior(self, populations=None):
        """
        Prior probability for each sub-population.
        """

        size = len(populations or self._population.populations)
        if size == 0:
            raise ValueError('empty list of populations')
        return np.ones(size, dtype=float) * (1.0 / size)

    def parent(self, ind, **kwargs):
        """
        Classify individual in one of the parental populations.

        Returns:
            The id or index (if no id is defined) of the selected
            population.

        See Also:
            :meth:`parent_prob` for a better description of the algorithm. This
            method accept the same arguments.
        """

        return self.parent_prob(ind, **kwargs).mode()

    def parent_prob(self, ind, *, prior=None, populations=None):
        """
        Classify individual in one of parental populations and return a
        probability distribution over population ids.

        Args:
            ind:
                a :class:`Individual` instance or an array with genotype data.
            prior:
                optional list of prior probabilities for each population.
            populations:
                optional list of parental populations or a MultiPopulation
                object. If no parental is given, it uses pop.populations.

        Returns:
            A :class:`Prob` mapping from population ids to probabilities.
        """

        pop = self._population
        prior = np.asarray(self.prior(populations))
        populations = self.populations(populations)

        if pop.is_biallelic:
            freqs = np.array([pop.freqs_vector for pop in populations])
            probs = likelihood.bayes(ind.data, freqs, prior=prior)
            data = {}
            for i, pop in enumerate(populations):
                data[pop.id or i] = probs[i]
            return Prob(data)
        else:
            raise NotImplementedError('only biallelic data is supported')


def ancestry_labels(pop):
    """
    Return a list of labels from a population object using subpopulation ids as
    data.
    """

    labels = []
    for i, subpop in enumerate(pop.populations):
        label = subpop.id or i
        labels.extend([label] * subpop.size)
    return labels

from .attr import Attr
from ..classifiers import SklearnClassifier
from ..libs import np
from ..libs import sk_naive_bayes, sk_svm
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
        raise ValueError('invalid method: %r' % which)

    def _normalize_labels(self, labels=None):
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

    def naive_bayes(self, labels=None, data='count', prior='uniform',
                    alpha=0.5):
        """
        Classify objects using the naive_bayes classifier.


        Args:
            labels:
                List of labels or a string with the metadata column used as
                label. Optionally, the 'ancestry' string classify using the
                sub-populations as labels.
            alpha:
                Additive (Laplace/Lidstone) smoothing parameter (0 for no
                smoothing).
            prior:
                The prior probability for each label. Must be either a Prob()
                object, the string 'uniform' or None. The default value is
                'uniform' that assigns a fixed uniform prior. If prior is None,
                it learns priors from data. Finally, it can also be specified
                as a Prob() object or a mapping from labels to probabilities.
        """

        kwargs = {'alpha': alpha}

        # Prepare prior data
        if prior == 'uniform':
            kwargs['fit_prior'] = False
            kwargs['class_prior'] = None
        elif prior is None:
            kwargs['fit_prior'] = True
            kwargs['class_prior'] = None
        else:
            labels = self._normalize_labels(labels)
            prob = Prob(prior)
            prob_vector = np.zeros(len(prob))
            for i, label in enumerate(sorted(set(labels))):
                prob_vector[i] = prob[label]
            kwargs['fit_prior'] = True
            kwargs['class_prior'] = prob_vector

        if data == 'count':
            classifier = sk_naive_bayes.MultinomialNB
        else:
            raise ValueError(
                'naive bayes only accepts "count" and "flat" for the data '
                'argument'
            )
        return self.sklearn(classifier, labels, data=data, **kwargs)

    def svm(self, labels=None, data='count', **kwargs):
        """
        Classify objects using the Support Vector Machine (SVM) classifier.
        """

        classifier = sk_svm.SVC
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

import numpy as np

from .attr import Attr
from ..result import transform_result


class Clusterization(Attr):
    """
    Implements all clustering methods for a population in a separate namespace.
    """

    _data = property(lambda _: np.array(_.population))
    _methods = {'kmeans', }

    def __call__(self, which, k=5, **kwargs):
        if callable(which):
            method = kwargs.pop('method', 'count')
            data = self._as_array(method)
            return transform_result(which, data, k, **kwargs)

        elif isinstance(which, str):
            which_ = which.lower().replace('-', '')
            method = getattr(self, which_)
            return method(k, **kwargs)

        raise ValueError('invalid method: %r' % which)

    def _as_array(self, data='count'):
        return self._population.as_array(data)

    def clusterize(self, which, k=5, **kwargs):
        """
        An alias to Population.cluster(...)
        """
        return cluster(which, k, **kwargs)

    def labels(self, which, k=5, **kwargs):
        """
        Clusterize using the given method and return a list of labels.
        """

    def sklearn(self, classifier, k=5, use_similarity=False, **kwargs):
        """
        Runs a scikit learn unsupervised classifier.
        """

    def lda(self, k, **kwargs):
        data = self._population.as_array('count')

    def missclassifications(self, which, labels='populations'):
        """
        Runs a clustering algorithm and compare with the expected labels.

        Args:
            which:
                A string describing the clusterization method.
            labels:
                An optional list of labels. If no label is given, it uses
                sub-population names or indexes.

        Returns:
            An object with classifications.
        """

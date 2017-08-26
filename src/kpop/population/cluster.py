class cluster:
    """
    Implements all clustering methods for a population in a separate namespace.
    """

    _data = property(lambda _: np.array(_.population))
    _methods = {'kmeans',}

    def __init__(self, population):
        self._population = population

    def __call__(self, which, k=5, **kwargs):
        pass

    def _as_array(self, data='count'):
        return self._population.as_array(data)

    def cluster(self, which, k=5, **kwargs):
        """
        An alias to Population.cluster(...)
        """
        return cluster(which, k, **kwargs)

    def labels(self, which, k=5, **kwargs):
        """
        Clusterize using the given method
        Args:
            which:
            k:
            **kwargs:

        Returns:

        """
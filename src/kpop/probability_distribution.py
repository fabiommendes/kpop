import collections
from math import log
from random import random

import numpy as np


class Prob(collections.Mapping):
    """
    A dictionary-like object that behaves as a mapping between categories to
    their respective probabilities.
    """

    __slots__ = ('_data',)

    def __init__(self, data, normalize=True, support=None):
        self._data = dict(data)
        if normalize:
            norm = sum(self._data.values())
            if norm != 1:
                for k, v in self._data.items():
                    self._data[k] = v / norm
        if support:
            self.update_support(support)

    def __getitem__(self, key):
        return self._data[key]

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self._data)

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def update_support(self, support):
        """
        Force all elements in support to be explicitly present in distribution
        (possibly with null probability).

        Args:
            support:
                a list of elements in the support set for probability
                distribution.
        """
        setdefault = self._data.setdefault
        for elem in support:
            setdefault(elem, 0.0)

    def set_support(self, support):
        """
        Defines the support set of distribution.

        If elements exist in support, they are forced to exist in distribution,
        possibly with null probability. If element exists in the distribution
        but is not present in support, raises a ValueError.
        """

        support = set(support)
        for x in self:
            if x not in support:
                raise ValueError('%r not in support' % x)
        self.update_support(support)

    def entropy(self):
        """
        Return the Shannon entropy for the probability distribution.
        """

        return sum(-x * log(x) for x in self.values() if x)

    def random(self):
        """
        Returns a random element.
        """

        r = random()
        cum_prob = 0
        for elem, p in self.items():
            cum_prob += p
            if cum_prob >= r:
                return elem

    def random_sequence(self, size):
        """
        Returns a sequence of random elements.
        """

        r = self.random
        return [r() for _ in range(size)]

    def mode(self):
        """
        Return the element with the maximum probability.

        If more than one element shares the maximum probability, return an
        arbitrary value within this set.
        """

        p_mode = 0.0
        mode = None
        for elem, p in self.items():
            if p >= p_mode:
                mode = elem
        return mode

    def mode_set(self):
        """
        Return a set of elements that share the maximum probability.
        """

        p_mode = 0.0
        mode_set = set()
        for elem, p in self.items():
            if p > p_mode:
                mode_set = {elem}
            elif p == p_mode:
                mode_set.add(elem)
        return mode_set

    def kl_divergence(self, distr: collections.Mapping):
        """
        Return the Kullback-Leibler divergence with probability distribution.
        """

        prob = self._data.get
        divergence = 0.0
        visited = 0

        for k, q in distr.items():
            visited += 1
            p = prob(k, 0.0)
            divergence += p and p * log(p / q)

        if len(self._data) != visited:
            return float('inf')

        return divergence

    def encode(self, coding):
        """
        Encode probability distribution as a vector.

        Args:
            coding: a sequence of ordered categories.

        Example:
            >>> prob = Prob({'a': 0.75, 'b': 0.25})
            >>> prob.encode(['b', 'a'])
            [0.25, 0.75]
        """

        prob = self._data.get
        return np.array([prob(x, 0.0) for x in coding], dtype=float)

import collections
from collections import defaultdict
from math import log
from random import random

import numpy as np


class Prob(collections.Mapping):
    """
    A dictionary-like object that behaves as a mapping between categories to
    their respective probabilities.
    """

    __slots__ = ('_data',)

    @classmethod
    def mixture(cls, coeffs, probs):
        """
        Create a mixture probability from the given coeffs and list of Probs
        objects.

        Args:
            coeffs:
                Mixture coefficients. These coefficients do not have to be
                normalized.
            probs:
                List of Prob objects.

        Returns:
            A Prob object representing the mixture.
        """
        if len(coeffs) != len(probs):
            raise ValueError('coeffs and probs must be aligned')
        data = defaultdict(float)

        for q, prob in zip(coeffs, probs):
            for k, p in prob.items():
                data[k] += q * p
        return Prob(data)

    def __init__(self, data, normalize=True, support=None):
        try:
            self._data = dict(data)
        except TypeError:
            self._data = dict(enumerate(data))

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

    def __eq__(self, other):
        if isinstance(other, collections.Mapping):
            nonzero_a = {k: v for k, v in self.items() if v}
            nonzero_b = {k: v for k, v in other.items() if v}
            return nonzero_a == nonzero_b
        return NotImplemented

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
        possibly with zero probability. If element exists in the distribution
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
        raise ValueError('probability distribution do not sum to 1.0')

    def random_sequence(self, size):
        """
        Returns a sequence of random elements.
        """

        r = self.random
        return [r() for _ in range(size)]

    def max(self):
        """
        Return the value of maximum probability.
        """

        return max(self.values())

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
                p_mode = p
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
                p_mode = p
            elif p == p_mode:
                mode_set.add(elem)
        return mode_set

    def sharp(self, mode_set=True):
        """
        Return a sharp version of the probability distribution.

        All elements receive probability zero, except the mode which receives
        probability one.
        """

        data = {k: 0.0 for k in self}
        if mode_set:
            mode_set = self.mode_set()
            p_mode = 1. / len(mode_set)
            for k in mode_set:
                data[k] = p_mode
        else:
            data[self.mode()] = 1.0
        return Prob(data)

    def kl_divergence(self, q: collections.Mapping):
        """
        Return the Kullback-Leibler divergence with probability distribution.

        This is given by the formula:

            $KL = \sum_i p_i \ln \frac {p_i} {q_i},$

        in which p_i comes from the probability object and q_i comes from the
        argument.
        """

        prob = self._data.get
        divergence = 0.0
        visited = 0

        for k, q in q.items():
            visited += 1
            p = prob(k, 0.0)
            if p:
                try:
                    divergence += p and p * log(p / q)
                except ZeroDivisionError:
                    return float('inf')

        if len(self._data) != visited:
            return float('inf')

        return divergence

    def encode(self, coding=None):
        """
        Encode probability distribution as a vector.

        Args:
            coding: a sequence of ordered categories.

        Example:
            >>> prob = Prob({'a': 0.75, 'b': 0.25})
            >>> prob.encode(['b', 'a'])
            [0.25, 0.75]
        """

        if coding is None:
            types = {type(x) for x in self._data}
            if types == {int}:
                coding = range(max(self) + 1)
            else:
                coding = sorted(self._data)

        prob = self._data.get
        return np.array([prob(x, 0.0) for x in coding], dtype=float)

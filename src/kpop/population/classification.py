import numpy as np

from .attr import Attr
from ..admixture import likelihood
from ..prob import Prob


class Classification(Attr):
    """
    Implements the population.classification attribute.
    """

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

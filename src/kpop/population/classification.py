from kpop.admixture import likelihood

from .attr import Attr
from ..prob import Prob


class Classification(Attr):
    """
    Implements the population.classification attribute.
    """

    def classify(self, ind, prior=None):
        """
        Classify individual in one of the parent populations.

        Args:
            ind: a :class:`Individual` instance.

        Returns:
            The id or index (if no id is defined) of the selected
            population.
        """

        return self.prob_classify(ind, prior=prior).mode()

    def prob_classify(self, ind, prior=None):
        """
        Classify individual in one of parental populations and return a
        probability distribution over population ids.

        Args:
            ind: a :class:`Individual` instance.

        Returns:
            A :class:`Prob` mapping from population ids to probabilities.
        """

        if prior is None:
            prior = self.prior()

        if self.is_biallelic:
            probs = likelihood.bayes(ind.data, self.freqs_vector_pop,
                                     prior=prior)
            data = {}
            for i, pop in enumerate(self.populations):
                data[pop.id or i] = probs[i]
            return Prob(data)
        else:
            raise NotImplementedError

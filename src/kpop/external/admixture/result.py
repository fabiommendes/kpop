from lazyutils import delegate_to, lazy

from kpop import MultiPopulation, Population
from kpop.prob import Prob


class AdmixtureResult:
    """
    Store ADMIXTURE _results.
    """

    num_iter = delegate_to('_parsed')
    log_like = delegate_to('_parsed')
    fst = delegate_to('_parsed')

    @lazy
    def num_loci(self):
        return self.freqs_vector_pop.shape[1]

    @lazy
    def num_individuals(self):
        return len(self.q)

    @lazy
    def _parsed(self):
        # parsed = AttrDict()
        raw = self._raw
        del self._raw
        print(raw)
        raise NotImplementedError

    @property
    def freqs(self):
        return self.freqs_vector_pop

    def __init__(self, out, freqs, q, ploidy):
        self._raw = out
        self.freqs_vector_pop = freqs
        self.q = q
        self.ploidy = ploidy
        self.num_alleles = 2

    def make_parental(self, labels=None):
        """
        Return a MultiPopulation() object that defines all parental
        frequencies.
        """

        parental = MultiPopulation(ploidy=self.ploidy,
                                   num_alleles=2,
                                   num_loci=self.num_loci)
        for i, freqs in enumerate(self.freqs_vector_pop):
            if labels is None:
                label = 'pop-{0!s}'.format((i + 1))
            else:
                label = labels[i]
            pop = Population(freqs=freqs, ploidy=self.ploidy, label=label)
            parental.add_population(pop)
        return parental

    def make_admixture_labels(self, individuals):
        """
        Return a copy of the list of individuals the with the admixture
        coefficients saved as Prob() object in the ind.admixture_q attribute.
        """

        result = []
        for ind, freqs in zip(individuals, self.q):
            prob = Prob(enumerate(freqs))
            ind = ind.copy(admixture_q=prob)
            result.append(ind)
        return result

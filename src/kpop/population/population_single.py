import numpy as np
from lazyutils import lazy

from kpop.population.utils import parse_population_data
from .individual import IndividualProxy, random_data_from_freqs_matrix
from .population_base import PopulationBase
from .populations import ImmutablePopulationList
from ..utils.frequencies import random_frequencies


class Population(PopulationBase):
    """
    A Population is a collection of individuals.
    """

    _multi_population_class = None

    @lazy
    def populations(self):
        return ImmutablePopulationList([self])

    @classmethod
    def random(cls, size=0, num_loci=0, alleles=2, ploidy=2, id=None,
               min_prob=0.0, seed=None):
        """
        Creates a new random population.

        Args:
            size:
                Number of individuals.
            num_loci:
                Number of loci in the genotype.
            alleles:
                Number of alleles for all loci.
            ploidy:
                Ploidy of genotype.
            min_prob:
                Minimum value for a frequency probability.

        Returns:
            A new population object.
        """
        if num_loci <= 0:
            raise ValueError('num_loci must be at least one!')
        if seed is not None:
            np.random.seed(seed)

        # Frequencies
        freqs = random_frequencies(num_loci, alleles, clip=min_prob, seed=seed)

        # Create data
        data = []
        for _ in range(size):
            ind = random_data_from_freqs_matrix(freqs, ploidy=ploidy)
            data.append(ind)
        data = np.array(data)

        # Return population
        return Population(
            data, freqs=freqs, id=id,
            num_loci=num_loci, num_alleles=alleles, ploidy=ploidy
        )

    def __init__(self, data=(), id=None, individual_ids=None, **kwargs):
        if isinstance(data, str):
            data, _labels = parse_population_data(data)
            individual_ids = individual_ids or _labels
        elif isinstance(data, PopulationBase):
            individual_ids = individual_ids or data.individual_ids
            data = data.as_array()
        else:
            data = np.array(data, dtype='uint8')

        self._data = np.asarray(data)
        kwargs.update(id=id, individual_ids=individual_ids)
        super(Population, self).__init__(**kwargs)

    def __len__(self):
        return len(self._data)

    def __getitem__(self, idx):
        return IndividualProxy(self, idx)

    def __iter__(self):
        for idx in range(self.size):
            yield IndividualProxy(self, idx)

    def __str__(self):
        return '\n'.join(str(x) for x in self)

    def __add__(self, other):
        if isinstance(other, Population):
            return self._multi_population_class([self, other])
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Population):
            x = self._data
            y = other._data
            return x.shape == y.shape and (x == y).all()
        else:
            return super().__eq__(other)

    def _as_array(self):
        return np.array(self._data)

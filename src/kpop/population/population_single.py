import numpy as np
from lazyutils import lazy

from kpop.population.individual import Individual
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

        freqs = random_frequencies(num_loci, alleles=alleles, clip=min_prob,
                                   seed=seed)
        pop = Population(freqs=freqs, id=id, num_loci=num_loci,
                         num_alleles=alleles, ploidy=ploidy)
        pop.fill(size, seed=seed)
        pop.freqs_matrix = freqs
        return pop

    def __init__(self, data=(), copy=False, **kwargs):
        self._data = []

        # Initialize data
        if isinstance(data, str):
            for ind in data.splitlines():
                self.add(ind, copy=copy)
        elif isinstance(data, PopulationBase):
            for ind in data:
                self.add(ind, copy=True)
        else:
            for ind in data:
                self.add(ind, copy=copy)

        super(Population, self).__init__(**kwargs)

        # Recompute ids, if they are not given
        if self.id:
            for ind in self._data:
                if getattr(ind, 'id', None) is None:
                    ind.id = self._next_id()
                ind.population = self
                ind._container = self

    def __len__(self):
        return len(self._data)

    def __getitem__(self, idx):
        return self._data[idx]

    def __iter__(self):
        return iter(self._data)

    def __str__(self):
        return '\n'.join(str(x) for x in self)

    def __add__(self, other):
        if isinstance(other, Population):
            return self._multi_population_class([self, other])
        return NotImplemented

    def add(self, ind, copy=False):
        """
        Adds individual to population.
        """

        if not isinstance(ind, Individual):
            ind = Individual(ind, population=self)
        else:
            if ind._container not in (self, None):
                if copy:
                    ind = ind.copy()
                else:
                    msg = 'individual already belongs to population'
                    raise ValueError(msg)

                ind.population = self
                ind._container = self

        self._data.append(ind)

    def remove(self, idx=-1):
        """
        Remove and return individual at given index.

        If no index is given, remove the last individual.
        """

        ind = self._data.pop(idx)
        ind._container = None
        return ind

    def fill(self, size, seed=None):
        """
        Fill population with new random individuals until it reaches the
        given size.
        """

        if seed is not None:
            np.random.seed(seed)

        while self.size < size:
            new = self.simulation.new_individual(id=self._next_id())
            self.add(new)

    def fill_missing(self):
        """
        Fills missing data.
        """

        freqs = self.freqs
        for ind in self:
            data = ind.data
            missing = np.where(data == 0)
            for i, j in zip(*missing):
                data[i, j] = freqs[i].random_individual()

import collections
from random import randrange

import numpy as np
from lazyutils import lazy

from kpop import Individual
from kpop.probability_distribution import Prob


class PopulationBase(collections.Sequence):
    """
    Base class for Population and Multipopulation.
    """

    @lazy
    def freqs(self):
        if self._freqs is None:
            self._freqs = self.empirical_freqs()
        return self._freqs

    @property
    def freqs_matrix(self):
        data = np.zeros((self.num_loci, self.num_alleles), dtype=float)
        pass

    @lazy
    def is_biallelic(self):
        pass

    @lazy
    def num_alleles(self):
        return max(len(freq) for freq in self.freqs)

    @lazy
    def ploidy(self):
        return self[0].ploidy

    @lazy
    def num_loci(self):
        return self[0].num_loci

    @property
    def size(self):
        return len(self)

    def __str__(self):
        return self.render(6)

    def render(self, label_align=None):
        """
        Renders population data to string.
        """

        return '\n'.join(ind.render(label_align=label_align) for ind in self)

    def save(self, file, **kwargs):
        """
        Saves population data to file.
        """

        data = self.render(**kwargs)
        with open(file, 'w') as F:
            F.write(data)

    def empirical_freqs(self, alpha=0.0):
        """
        Return an array with an empirical estimate of population frequencies.
        """

        counters = [collections.Counter() for _ in range(self.num_loci)]
        for ind in self:
            for counter, genotype in zip(counters, ind):
                counter.update(genotype)

        if alpha == 0:
            return [Prob(c) for c in counters]
        else:
            raise NotImplementedError

    def genetic_drift(self, n_generations, population_size=None):
        """
        Applies a simple genetic drift model for mutating the frequencies in
        freqs. Returns a new empty population with the new frequencies.

        Args:
            n_generations:
                Number of generations to apply the drift model
            population_size:
                Effective size of population. Defaults to the current population
                size.

        Example:
            >>> pop = Population.from_freqs(np.ones(100) / 2)
            >>> pop2 = pop.genetic_drift(20, population_size=10)
            >>> pop2.fill(10)

            After many generations, genetic drift in a small population tends
            to fixate alleles. The above case has a change greater than 50%
            chance of fixating each allele.

            >>> fixed = fout * (1 - fout) == 0
            >>> fixed.any() # it is random, but we can be pretty sure that at least 1
            ...             # allele has been fixed!
            True

        """

        raise NotImplementedError
        N = self.ploidy * population_size
        fout = unfillfreqs(freqs).copy()
        for _ in range(n_generations):
            for j, f in enumerate(fout):
                fout[j] = np.random.binomial(N, f) / N
        return asfreqshape(fout, freqs.shape)


class Population(PopulationBase):
    """
    A Population is a collection of individuals
    """

    @classmethod
    def random(cls, size, num_loci, alleles=2, ploidy=2, label_prefix='ind'):
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
            label_prefix:
                A string prefix to prepend to all individual labels.

        Returns:
            A new population object.
        """
        freqs = random_frequencies(num_loci, alleles=alleles)
        pop = Population(freqs=freqs)

        data = random_pop_data(size, freqs, ploidy=ploidy)
        for i, ind in enumerate(data, start=1):
            label = '%s%s' % (label_prefix, i)
            ind = Individual(ind, copy=False, population=pop, label=label)
            pop.add(ind)
        return pop

    def __init__(self, data=(), freqs=None):
        self._data = list(data)
        self._freqs = normalize_freqs(freqs)
        self.allele_names = None

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
            return MultiPopulation((self, other))
        return NotImplemented

    def new(self):
        """
        Return a new random individual from population
        """

    def new_offspring(self, i=None, j=None):
        """
        Return a new offspring created by breeding two random elements in the
        population.
        """

        size = len(self)
        if i is None:
            i = randrange(size)
        if j is None:
            j = randrange(size)
        return self[i].breed(self[j])

    def add(self, ind):
        """
        Adds individual to population.
        """

        if not isinstance(ind, Individual):
            ind = Individual(ind, population=self)
        else:
            ind.population = self
        self._data.append(ind)


def normalize_freqs(freqs):
    """
    Normalize frequencies to be a sequence of probability distributions.
    """

    if freqs is None:
        return None

    if isinstance(freqs, np.ndarray):
        result = []
        for prob in freqs:
            data = enumerate(prob)
            result.append(Prob(data))
        return result
    else:
        result = []
        for prob in freqs:
            if not isinstance(prob, collections.Mapping):
                prob = enumerate(prob)
            prob = Prob(prob)
            result.append(prob)
        return result


class MultiPopulation(PopulationBase):
    """
    A population formed by several sub-populations.
    """

    def __init__(self, populations=()):
        self.populations = list(populations)

    def __len__(self):
        return sum(len(x) for x in self.populations)

    def __getitem__(self, idx):
        if isinstance(idx, int):
            i = idx
            for pop in self.populations:
                size = len(pop)
                if i >= size:
                    i -= size
                else:
                    return pop[i]
            raise IndexError(idx)

    def __iter__(self):
        for pop in self.populations:
            yield from pop

    def __add__(self, other):
        if isinstance(other, Population):
            populations = self.populations + [other]
            return MultiPopulation(populations)
        elif isinstance(other, MultiPopulation):
            populations = self.populations + other.populations
            return MultiPopulation(populations)
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Population):
            populations = [other] + self.populations
            return MultiPopulation(populations)
        return NotImplemented


def random_frequencies(num_loci, alleles=2):
    """
    Return a random frequency distribution for SNPs.
    """

    if alleles <= 1:
        raise ValueError('needs at least 2 different alleles')
    if alleles == 2:
        return fill_frequencies(np.random.uniform(0, 1, size=num_loci))
    else:
        data = np.random.uniform(0, 1, size=num_loci * alleles)
        data = data.reshape((num_loci, alleles))
        data /= data.sum(axis=1)[:, None]
        return data


def fill_frequencies(freqs):
    """
    Return a Nx2 array with N tuples of (f, 1-f) for each element f of freqs.
    """

    size = len(freqs)
    out = np.ones((size, 2), dtype=float)
    out[:, 0] = freqs
    out[:, 1] -= freqs
    return out


def random_pop_data(size, freqs, ploidy=2, dtype=np.int8, alleles=2):
    """
    Create a random population data of bi-allele individuals from the given
    frequency.

    freqs can be a list of frequencies f[j] or a list of (f[j], 1-f[j]) tuples.
    In both cases, the frequency f[j] represents the probability that the
    j-th feature has a value of 1.
    """

    num_loci = len(freqs)
    data = np.random.uniform(size=num_loci * ploidy * int(size))
    data = data.reshape((size, num_loci, ploidy))
    alleles = freqs.shape[1]
    if alleles == 2:
        return np.asarray(data < freqs[:, 0, None], dtype=dtype)
    else:
        data = np.zeros((size, num_loci, ploidy), dtype=dtype)
        mask = data >= freqs[:, 0, None]
        for i in range(1, alleles):
            data[mask] = i
            mask &= data >= freqs[:, i, None]
        return data


popA = Population.random(5, 20)
popB = Population.random(5, 20)
print(popA + popB)

from pprint import pprint
pprint(popA.empirical_freqs())

print(popA.freqs)
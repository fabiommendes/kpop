import numpy as np
from lazyutils import lazy

from kpop.population.population_base import PopulationBase
from kpop.population.population_single import Population
from kpop.population.populations import PopulationList


class MultiPopulation(PopulationBase):
    """
    A population formed by several sub-populations.
    """

    is_multi_population = True

    @lazy
    def freqs_vector_pop(self):
        return np.array([pop.freqs_vector for pop in self.populations])

    @lazy
    def freqs_matrix_pop(self):
        return np.array([pop.freqs_matrix for pop in self.populations])

    @lazy
    def freqs_pop(self):
        return [pop.freqs for pop in self.populations]

    def prior(self):
        """
        Prior probability for each sub-population.
        """

        size = len(self.populations)
        if size == 0:
            return np.array([], dtype=float)
        return np.zeros(size, dtype=float) + 1 / size

    def __init__(self, populations=(), freqs=None, **kwargs):
        if freqs is not None:
            raise ValueError('cannot specify frequencies on MultiPopulation '
                             'objects')
        self.populations = PopulationList()
        for population in populations:
            self.add_population(population)
        super().__init__(**kwargs)

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
        if isinstance(other, (Population, MultiPopulation)):
            populations = list(self.populations)
            return MultiPopulation(populations + list(other.populations))
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, (Population, MultiPopulation)):
            populations = list(self.populations)
            return MultiPopulation([other] + populations)
        return NotImplemented

    def add_population(self, population):
        """
        Adds a new sub-population.

        Args:
            population: A :class:`Population` instance.
        """

        if population in self.populations:
            raise ValueError('sub-population already present!')
        self.populations.append(population)

    def slice_indexes(self, indexes):
        """
        Map indexes to a list of indexes for each sub-population.
        """

        indexes = [[] for pop in self.populations]
        breakpoints = np.add.accumulate([len(pop) for pop in self.populations])
        for idx in indexes:
            for k, idx_max in enumerate(breakpoints):
                if idx < idx_max:
                    if k > 1:
                        idx -= breakpoints[k - 1]
                    indexes[k].append(idx)
                    break
            else:
                raise IndexError('index out of bounds: %s' % idx)

        return [np.array(lst) for lst in indexes]

    def copy(self, populations=None, **kwargs):
        new = super().copy(**kwargs)
        if populations is not None:
            new.populations = PopulationList(populations)
        return new

    def keep_individuals(self, indexes, **kwargs):
        sliced_indexes = self.slice_indexes(indexes)
        data = zip(self.populations, sliced_indexes)
        populations = [pop.keep_individuals(indexes) for pop, indexes in data]
        return self.copy(populations=populations, **kwargs)

    def keep_loci(self, indexes, **kwargs):
        sliced_indexes = self.slice_indexes(indexes)
        data = zip(self.populations, sliced_indexes)
        populations = [pop.keep_loci(indexes) for pop, indexes in data]
        return self.copy(populations=populations, **kwargs)

    def map_alleles(self, alleles_mapping, **kwargs):
        map = alleles_mapping
        populations = [pop.map_alleles(map) for pop in self.populations]
        return self.copy(populations=populations, **kwargs)

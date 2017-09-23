from .population_base import PopulationBase
from .population_list import PopulationList
from .population_single import Population
from ..libs import np, pd


class MultiPopulation(PopulationBase):
    """
    A population formed by several sub-populations.
    """

    is_multi_population = True

    @property
    def meta(self):
        # Collect common columns
        cols, *tail = [set(pop.meta.columns) for pop in self.populations]
        for new in tail:
            cols.intersection_update(new)
        cols = sorted(cols)

        # Create dataframe
        metas = [pop.meta[cols] for pop in self.populations]
        df = pd.DataFrame(columns=cols)
        for meta in metas:
            df = df.append(meta)
        return df

    def __init__(self, populations=(), freqs=None, **kwargs):
        if freqs is not None:
            raise ValueError(
                'cannot specify frequencies on MultiPopulation objects'
            )

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

    def _as_array(self):
        stack = [pop._as_array() for pop in self.populations]
        return np.vstack(stack)

    def add_population(self, population):
        """
        Adds a new sub-population.

        Args:
            population: A :class:`Population` instance.
        """

        for sub in self.populations:
            if population.id is not None and population.id == sub.id:
                raise ValueError(
                    'sub-population with id "%s" already present!' % sub.id
                )
        self.populations.append(population)

    def slice_indexes(self, indexes):
        """
        Map indexes to a list of indexes for each sub-population.
        """

        slices = [[] for pop in self.populations]
        breakpoints = np.add.accumulate([len(pop) for pop in self.populations])
        for idx in indexes:
            for k, idx_max in enumerate(breakpoints):
                if idx < idx_max:
                    if k > 0:
                        idx -= breakpoints[k - 1]
                    slices[k].append(idx)
                    break
            else:
                raise IndexError('index out of bounds: %s' % idx)

        return [np.array(lst) for lst in slices]

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
        populations = [pop.keep_loci(indexes) for pop in self.populations]
        return self.copy(populations=populations, **kwargs)

    def map_alleles(self, alleles_mapping, **kwargs):
        map = alleles_mapping
        populations = [pop.map_alleles(map) for pop in self.populations]
        return self.copy(populations=populations, **kwargs)

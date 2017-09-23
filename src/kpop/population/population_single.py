from lazyutils import lazy

from .individual import IndividualProxy
from .population_base import PopulationBase
from .population_list import ImmutablePopulationList
from .utils import parse_population_data, discard_attrs
from ..libs import np, pd


class Population(PopulationBase):
    """
    A Population is a collection of individuals.
    """

    _multi_population_class = None

    @lazy
    def populations(self):
        return ImmutablePopulationList([self])

    def __init__(self, data=(), id=None, individual_ids=None, **kwargs):
        # Prepare data
        if isinstance(data, str):
            data, _labels = parse_population_data(data)
            individual_ids = individual_ids or _labels
        elif isinstance(data, PopulationBase):
            individual_ids = individual_ids or data.individual_ids
            data = data.as_array()
        elif len(data) != 0 and getattr(data[0], 'is_individual', False):
            if individual_ids is None:
                individual_ids = [ind.id for ind in data]
            if set(individual_ids) == {None}:
                individual_ids = None
            data = np.array([ind.data for ind in data])
        else:
            data = np.array(data, dtype='uint8')

        # Save values
        self._data = np.asarray(data)
        kwargs.update(id=id)
        super(Population, self).__init__(**kwargs)

        # Individual ids and meta data
        if individual_ids is None:
            fmt = ('ind' if id is None else id) + '%s'
            individual_ids = [fmt % (i + 1) for i in range(self.size)]
        self.individual_ids = list(individual_ids)
        self.meta = pd.DataFrame(index=individual_ids)

    def __len__(self):
        return len(self._data)

    def __getitem__(self, idx):
        return IndividualProxy(self, idx)

    def __iter__(self):
        for idx in range(self.size):
            yield IndividualProxy(self, idx)

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

    def copy(self, data=None, **kwargs):
        new = super().copy(**kwargs)

        if data is not None:
            data = np.asarray(data, dtype='uint8')
            if data.ndim != 3:
                raise ValueError('invalid number of dimensions: %s' % data.ndim)
            new._data = data

            if data.shape != self.shape:
                discard_attrs(new, new._shape_attrs)

        return new

    def keep_loci(self, indexes, **kwargs):
        check_data(kwargs)
        new_data = self._data[:, indexes, :]
        return self.copy(data=new_data, **kwargs)

    def keep_individuals(self, indexes, **kwargs):
        check_data(kwargs)
        new_data = self._data[indexes]
        return self.copy(data=new_data, **kwargs)

    def map_alleles(self, alleles_mapping, **kwargs):
        check_data(kwargs)
        data = self._data.copy()
        for i, ind in enumerate(data):
            for j, locus in enumerate(ind):
                map = alleles_mapping[j]
                for k, value in enumerate(locus):
                    data[i, j, k] = map.get(value, value)

        return self.copy(data=data, **kwargs)


def check_data(kwargs):
    if 'data' in kwargs:
        raise TypeError('invalid argument: "data"')

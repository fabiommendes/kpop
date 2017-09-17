import collections


class PopulationList(collections.MutableSequence):
    """
    Class for the .population attribute of a MultiPopulation instance.
    """

    def __init__(self, data=()):
        self._data = list(data)

    def __len__(self):
        return len(self._data)

    def __setitem__(self, index, value):
        index = self._conv_index(index)
        self._data[index] = value

    def __delitem__(self, index):
        index = self._conv_index(index)
        del self._data[index]

    def __getitem__(self, index):
        index = self._conv_index(index)
        return self._data[index]

    def __repr__(self):
        return repr(self._data)

    def __eq__(self, other):
        if isinstance(other, str):
            return NotImplemented

        if isinstance(other, collections.Sequence):
            zipped = zip(self, other)
            return len(self) == len(other) and all(x == y for x, y in zipped)

        return NotImplemented

    def insert(self, index, value):
        index = self._conv_index(index)
        self._data.insert(index, value)

    def _conv_index(self, idx):
        if isinstance(idx, (int, slice)):
            return idx
        else:
            for i, pop in enumerate(self):
                if idx == pop.id:
                    return i
            else:
                raise KeyError('invalid id: %r' % idx)

    def population_ids(self, prefix='pop'):
        """
        Return a list of labels or indexes for sub-populations.
        """

        def gen():
            for i, pop in enumerate(self):
                yield (prefix + str(i) if pop.id is None else pop.id)

        return list(gen())


class ImmutablePopulationList(PopulationList):
    """
    Immutable version of PopulationList.
    """

    def __setitem__(self, idx, value):
        raise self._immutable()

    def __delitem__(self, idx):
        raise self._immutable()

    def insert(self, index, value):
        raise self._immutable()

    def _immutable(self):
        return TypeError('population list is immutable')

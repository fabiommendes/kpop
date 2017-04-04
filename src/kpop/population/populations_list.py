import collections


class PopulationsList(collections.MutableSequence):
    """
    Class for the .population attribute of a MultiPopulation instance.
    """

    def __init__(self, data=()):
        self._data = list(data)

    def __len__(self):
        return len(self._data)

    def __setitem__(self, index, value):
        index = self._conv_index(index)
        self._check_new_pop(value)
        self._data[index] = value

    def __delitem__(self, index):
        index = self._conv_index(index)
        del self._data[index]

    def __getitem__(self, index):
        index = self._conv_index(index)
        return self._data[index]

    def __repr__(self):
        return repr(self._data)

    def insert(self, index, value):
        index = self._conv_index(index)
        self._check_new_pop(value)
        self._data.insert(index, value)

    def _check_new_pop(self, value):
        if self._data:
            pass

    def _conv_index(self, idx):
        if isinstance(idx, (int, slice)):
            return idx
        else:
            for i, pop in enumerate(self):
                if idx == pop.label:
                    return i
            else:
                raise KeyError('invalid label: {0!r}'.format(idx))

    def labels(self, prefix='pop'):
        """
        Return a list of labels or indexes for sub-populations.
        """

        def gen():
            for i, pop in enumerate(self):
                yield (prefix + str(i) if pop.label is None else pop.label)
        return list(gen())


class ImmutablePopulationList(PopulationsList):
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

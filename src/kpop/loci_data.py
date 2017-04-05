import collections
from collections import defaultdict
from maputil import InvMap

null_getter = defaultdict(lambda: None)


class Locus:
    """
    Represents a single locus of genetic data.
    """

    def __init__(self, name, num_alleles=2, chromosome=None, position=None,
                 genetic_distance=None, names=None):
        self.name = name
        self.num_alleles = num_alleles
        self.chromosome = chromosome
        self.position = position
        self.genetic_distance = genetic_distance
        self.allele_names = InvMap(names or {})
        self._has_complex_names = any(len(x) > 1 for x in self.allele_names)

    def render_locus(self, data):
        """
        Renders data in given locus.
        """

        components = self.render_components(data)
        if self._has_complex_names:
            return ','.join(components)
        return ''.join(components)

    def render_components(self, data):
        """
        Return a list of strings associated with data.
        """

        map = self.allele_names
        return [map[x] for x in data]

    def register_allele_name(self, name, value=None):
        """
        Register a new allele name.

        Args:
            value:
                Value associated with the name. If no value is given, tries to
                find a suitable value.
            name:
                Name associated with value.
        """

        if len(name) > 1:
            self._has_complex_names = True

        if value is None:
            raise NotImplementedError
        self.allele_names[value] = name

    def __repr__(self):
        return '{0!s}({1!r}, num_alleles={2!r})'.format(self.name, self.num_alleles)

    def render_map(self):
        """
        Render locus data as a line in a .MAP file.
        """

        return '{0!s} {1!s} {2!s} {3!s}'.format(
            self.chromosome or 0,
            self.name,
            self.position or 0,
            self.genetic_distance or -1
        )


class LociSequence(collections.Sequence):
    """
    A sequence of loci.
    """

    def __init__(self, names, num_alleles=None, chromosomes=None,
                 positions=None, genetic_distances=None):

        if isinstance(names, int):
            names = ['marker{0!s}'.format(i) for i in range(names)]

        self._data = data = []
        self._names_map = names_map = {}
        chromosomes = chromosomes or null_getter
        positions = positions or null_getter
        genetic_distances = genetic_distances or null_getter
        num_alleles = num_alleles or defaultdict(lambda: 2)

        for i, name in enumerate(names):
            locus = Locus(name,
                          num_alleles=num_alleles[i],
                          chromosome=chromosomes[i],
                          position=positions[i],
                          genetic_distance=genetic_distances[i])
            data.append(locus)
            if name in names_map:
                raise ValueError('repeated marker name: {0!r}'.format(name))
            names_map[name] = i

    def __getitem__(self, index):
        if isinstance(index, str):
            idx = self._names_map[index]
            return self._data[idx]
        return self._data[index]

    def __len__(self):
        return len(self._data)

    def render_map(self):
        """
        Render loci data as a .MAP file.
        """

        return '\n'.join(locus.render_map() for locus in self)

    def render_data(self, data, raw=False):
        """
        Renders some data

        Args:
            data:
            raw:
                If True, render raw numeric values instead of converting them
                to their corresponding character representations.

        Returns:

        """

        return ' '.join(locus.render_locus(x) for locus, x in zip(self, data))

import os
import pickle

from kpop.libs import lazy_module
from .attr import Attr

io = lazy_module('kpop.io')


class Io(Attr):
    """
    Implements the population.io attribute.
    """

    @staticmethod
    def load(file, format='auto', **kwargs):
        """
        Loads population from file.

        It accepts a few formats:

        * 'pickle': Python serialization protocol
        * 'csv': comma separated values

        See Also:
            :func:`kpop.load_csv` for more options for the csv loader.
        """

        if format == 'auto':
            format = format_from_path(file)

        if format == 'pickle':
            if isinstance(file, str):
                with open(file, 'r+b') as F:
                    return pickle.load(F)
            else:
                return pickle.load(file)
        elif format == 'csv':
            return io.load_csv(file, **kwargs)
        else:
            raise NotImplementedError

    def csv(self, **kwargs):
        """
        Return population as a CSV string.
        """

        pop = self._population
        names = pop.allele_names
        if names is None:
            names = ['loci%s' % i for i in range(1, pop.size + 1)]
        header = 'id,' + ','.join(names)
        data = '\n'.join(x.render_csv(**kwargs) for x in self._population)
        return '%s\n%s' % (header, data)

    def plink_ped(self):
        """
        Renders population as a PLINK's .ped string.
        """

        lines = []
        memo = {}
        for i, ind in enumerate(self._population, start=1):
            line = ind.render_ped(individual_id=i, memo=memo)
            lines.append(line)
        return '\n'.join(lines)

    def plink_map(self):
        """
        Renders population .map file for use with PLINK.
        """

        data = []
        for j in range(1, self._population.num_loci + 1):
            data.append('1 snp%s 0 %s' % (j, j))
        return '\n'.join(data)

    def save(self, file, format='auto', **kwargs):
        """
        Saves population data to file.
        """

        if format == 'auto':
            format = format_from_path(file)

        if format == 'pickle':
            with open(file, 'w+b') as F:
                pickle.dump(self._population, F)
        elif format in ['csv', 'ped', 'map', 'plink_ped', 'plink_map']:
            attr = {'ped': 'plink_ped', 'map': 'plink_map'}
            render = getattr(self, attr.get(format, format))
            data = render(**kwargs)
            with open(file, 'w') as F:
                F.write(data)
        else:
            raise ValueError('invalid file format: %r' % format)

    def render(self, id_align=None, max_loci=None, max_ind=None):
        """
        Renders population data to string.

        Args:
            id_align:
                Size of column used to print id values. Ids are aligned to
                the right in this column. If not given, it tries to infer
                an optimal value.
            max_loci:
                Maximum number of loci printed before collapsing.
            max_ind:
                Maximum number of individuals printed before collapsing.
        """

        pop = self._population
        size = len(pop)
        if max_ind and size > max_ind:
            good_idx = set(range(max_ind // 2))
            good_idx.update(range(size - max_ind // 2, size))
        else:
            good_idx = set(range(size))

        # Find best align if no id align is set
        if id_align is None:
            id_align = max(len(x.id) for x in pop)

        # Render individuals
        data = [x.render(id_align=id_align, max_loci=max_loci)
                for i, x in enumerate(pop) if i in good_idx]

        # Add ellipsis for large data sets
        if max_ind and size > max_ind:
            data.insert(max_ind // 2 + 1, '...')

        return '\n'.join(data)


def format_from_path(path):
    if path.endswith('.pickle'):
        return 'pickle'
    elif path.endswith('.csv'):
        return 'csv'
    elif path.endswith('.ped'):
        return 'ped'
    else:
        raise ValueError('invalid type: %r' % os.path.splitext(path)[-1])

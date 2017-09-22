import pickle

from .attr import Attr


class Io(Attr):
    """
    Implements the population.io attribute.
    """

    @staticmethod
    def load(file, format='pickle'):
        """
        Loads population from file.
        """

        if format == 'pickle':
            if isinstance(file, str):
                with open(file, 'r+b') as F:
                    return pickle.load(F)
            else:
                return pickle.load(file)
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
        header = 'label,' + ','.join(names)
        data = '\n'.join(x.render_csv(**kwargs) for x in self)
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
            if file.endswith('.pickle'):
                return self.save(file, 'pickle', **kwargs)
            elif file.endswith('.csv'):
                return self.save(file, 'csv', **kwargs)
            else:
                raise ValueError('could not determine file format')

        if format == 'pickle':
            with open(file, 'w+b') as F:
                pickle.dump(self, F)
        elif format in ['csv', 'ped', 'map', 'plink_ped', 'plink_map']:
            attr = {'ped': 'plink_ped', 'map': 'plink_map'}
            render = getattr(self, attr.get(format, format))
            data = render(**kwargs)
            with open(file, 'w') as F:
                F.write(data)
        else:
            raise ValueError('invalid file format: %r' % format)

    def render(self, id_align=None, limit=None, ind_limit=None):
        """
        Renders population data to string.
        """

        pop = self._population
        size = len(pop)
        if ind_limit and size > ind_limit:
            good_idx = set(range(limit // 2))
            good_idx.update(range(size - limit // 2, size))
        else:
            good_idx = set(range(size))

        # Find best align if no label align is set
        if id_align == 'best':
            id_align = max(len(x.label) + 1 for x in pop)

        # Render individuals
        data = [x.render(id_align=id_align, limit=limit)
                for i, x in enumerate(pop) if i in good_idx]

        # Add ellipsis for large data sets
        if ind_limit and size > ind_limit:
            data.insert(ind_limit // 2 + 1, '...')

        return '\n'.join(data)

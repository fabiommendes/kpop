import abc
import pickle


class RenderablePopulationMixin(abc.ABC):
    """
    Rendering and serialization functions for population classes
    """

    # Default properties
    size = num_loci = ploidy = 0
    allele_names = {}

    @abc.abstractmethod
    def __len__(self):
        pass

    @abc.abstractmethod
    def __iter__(self):
        pass

    def render(self, label_align=None, limit=None, ind_limit=None):
        """
        Renders population data to string.
        """

        size = len(self)
        if ind_limit and size > ind_limit:
            good_idx = set(range(limit // 2))
            good_idx.update(range(size - limit // 2, size))
        else:
            good_idx = set(range(size))

        # Find best align if no label align is set
        if label_align == 'best':
            label_align = max(len(x.label) + 1 for x in self)

        # Render individuals
        data = [x.render(label_align=label_align, limit=limit)
                for i, x in enumerate(self) if i in good_idx]

        # Add ellipsis for large data sets
        if ind_limit and size > ind_limit:
            data.insert(ind_limit // 2 + 1, '...')

        return '\n'.join(data)

    def render_csv(self, **kwargs):
        """
        Return population as CSV data.
        """

        names = self.allele_names
        if names is None:
            names = ['loci{0!s}'.format(i) for i in range(1, self.size + 1)]
        header = 'label,' + ','.join(names)
        data = '\n'.join(x.render_csv(**kwargs) for x in self)
        return '{0!s}\n{1!s}'.format(header, data)

    def render_ped(self):
        """
        Renders population as a plink's .ped file.
        """

        lines = []
        memo = {}
        for i, ind in enumerate(self, start=1):
            line = ind.render_ped(individual_id=i, memo=memo)
            lines.append(line)
        return '\n'.join(lines)

    def render_map(self):
        """
        Renders population .map file for use with plink.
        """

        data = []
        for j in range(1, self.num_loci + 1):
            data.append('1 snp{0!s} 0 {1!s}'.format(j, j))
        return '\n'.join(data)

    def save(self, file, format='pickle', **kwargs):
        """
        Saves population data to file.
        """

        if format == 'auto':
            if file.endswith('.pickle'):
                return self.save(file, 'pickle', **kwargs)
            elif file.endswith('.csv'):
                return self.save(file, 'csv', **kwargs)

        if format == 'pickle':
            with open(file, 'w+b') as F:
                pickle.dump(self, F)
        elif format in ['csv', 'ped', 'map']:
            render = getattr(self, 'render_' + format)
            data = render(**kwargs)
            with open(file, 'w') as F:
                F.write(data)
        else:
            raise ValueError('invalid file format: {0!r}'.format(format))

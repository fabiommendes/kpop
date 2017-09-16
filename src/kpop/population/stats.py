import collections

import numpy as np

from .attr import Attr
from ..prob import Prob
from ..statistics import biallelic_pairwise_fst
from ..utils import freqs_to_matrix


class Stats(Attr):
    """
    Implements the population.stats attribute.
    """

    def allele_count(self, allele=1):
        """
        Return an array of (size, num_loci) with the counts for the number of
        times the given allele appears in each individual at each locus.
        """

        data = self._population.as_array('raw')
        return (data == allele).sum(axis=2)

    def fst_matrix(self, sizes=None):
        """
        Return the Fst divergence matrix for all sub-populations.

        See Also:
            :func:`kpop.statistics.biallelic_pairwise_fst`
        """

        pops = self._population.populations
        sizes = list(map(len, pops)) if sizes is None else sizes
        freqs = [pop.freqs_vector for pop in pops]
        hfreqs = [pop.hfreqs_vector for pop in pops]
        return biallelic_pairwise_fst(freqs, hfreqs, sizes)

    def render_fst(self, sizes=None):
        """
        Returns a pretty string with the Fst divergence for population.

        See Also:
            :func:`.fst_matrix`
        """

        data = self.fst_matrix(sizes)
        populations = self._population.populations
        names = [
            pop.id or str(i + 1) for i, pop in enumerate(populations)
            ]
        colsize = max(map(len, names))
        rownames = [name.ljust(colsize) for name in names[1:]]
        colwidth = max(colsize, 9)
        colnames = [name.rjust(colwidth) for name in names[:-1]]
        lines = [' ' * (colsize + 2) + ' '.join(colnames)]

        for i, rowname in enumerate(rownames):
            row = ['%.6f' % x for x in data[i + 1][:i + 1]]
            row = [x.rjust(colwidth) for x in row]
            row = '%s: %s' % (rowname, ' '.join(row))
            lines.append(row)

        return '\n'.join(lines)

    def empirical_freqs(self, alpha=0.0, as_matrix=False):
        """
        Return an array with an empirical estimate of population frequencies.

        Args:
            alpha:
                Not used
            as_matrix (bool):
                If True, return a frequency matrix instead of a list of Prob()
                instances.

        Returns:
            A list of Prob() instances with the distribution of each locus.
            If ``as_matrix=True`` return a frequency matrix.
        """
        pop = self._population
        counters = [collections.Counter() for _ in range(pop.num_loci)]
        for ind in pop:
            for counter, genotype in zip(counters, ind):
                counter.update([x for x in genotype if x])
        freqs = [Prob({k: v + alpha for k, v in c.items()}) for c in counters]

        # Normalize return value
        if as_matrix:
            return freqs_to_matrix(freqs, pop.num_alleles)
        else:
            return freqs

    def non_biallelic(self, *, _data=None):
        """
        Return a list of all non biallelic loci.
        """

        pop = self._population
        data = np.array([ind.data for ind in pop]) if _data is None else _data
        num_alleles = data.max(axis=(0, 2))
        return np.where(num_alleles > 2)[0]

    def render_biallelic_freqs(self, sep='  ', align=True,
                               decimal_places=6, **kwargs):
        """
        Return a string with a list of frequencies of the first allele in a
        CSV-compatible format.

        Each line corresponds to a different loci and each column is a
        different population.

        Args:
            sep (str):
                String used to separate each item in a line.
            align (bool):
                If True (default), align columns to the right.
            decimal_places (int):
                Force all numeric data to be rendered with the given number
                of decimal places.
        """

        from kpop.io import csv_lines

        pop = self._population
        loci_names = ['L%s' % (i + 1) for i in range(pop.num_loci)]
        columns = [pop.id for pop in pop.populations]
        columns.insert(0, 'locus')
        freqs = [pop.freqs_vector for pop in pop.populations]
        data = list(zip(loci_names, *freqs))

        # Call csv_lines with the prepared data
        kwargs['align'] = align
        kwargs['sep'] = sep
        kwargs['decimal_places'] = decimal_places
        lines = csv_lines(data, columns=columns, **kwargs)
        return ''.join(lines)

import collections
import copy

import numpy as np
from lazyutils import lazy
from sidekick import _
from sklearn import preprocessing

from .admixture import Admixture
from .classification import Classification
from .clustering import Clustering
from .io import Io
from .plot import Plot
from .projection import Projection
from .simulation import Simulation
from .stats import Stats
from .utils import normalize_freqs_arg, freqs_fastest, get_freqs, set_freqs, \
    hfreqs_vector
from ..utils.frequencies import fill_freqs_vector, freqs_to_matrix
from ..utils.functional import fn_lazy, fn_property


class PopulationBase(collections.Sequence):
    """
    Base class for Population and MultiPopulation.

    Attrs:
        freqs:
            A list of :cls:`kpop.Prob` objects representing the probabilities
            of each loci.
        freqs_matrix:
            A full matrix with the shape (num individuals, max num of alleles)
            with the probability for each allele.
        freqs_vector:
            Frequencies for allele 1. This is more useful for biallelic data,
            since the frequency of the second allele is simply the complement.
        hfreqs_vector:
            Vector of frequencies of heterozygotes.
    """

    # General shape
    size = property(len)
    num_loci = lazy(lambda _: _[0].num_loci)
    ploidy = lazy(lambda _: _[0].ploidy)
    shape = property(lambda _: (_.size, _.num_loci, _.ploidy))
    data_size = fn_property(_.size * _.num_loci * _.ploidy)
    dtype = np.dtype('uint8')

    # Frequencies
    freqs = property(get_freqs, set_freqs)
    freqs_matrix = lazy(lambda _: freqs_to_matrix(_.freqs))
    freqs_vector = lazy(lambda _: np.ascontiguousarray(_.freqs_matrix[:, 0]))
    hfreqs_vector = lazy(hfreqs_vector)
    _freqs_fastest = property(freqs_fastest)

    # Allele statistics
    allele_names = None
    is_biallelic = fn_lazy(_.num_alleles == 2)
    num_alleles = lazy(lambda _: max(max(freq) for freq in _.freqs))

    # Multi population
    is_multi_population = False
    num_populations = fn_property(lambda _: len(_.populations))

    # Missing data
    has_missing = property(lambda _: any(ind.has_missing for ind in _))
    missing_total = property(lambda _: sum(ind.missing_total for ind in _))
    missing_ratio = fn_property(_.missing_total / _.data_size)

    # Special attributes. These will be inserted later via monkey patching
    populations = ()
    admixture = Admixture()
    clustering = Clustering()
    classification = Classification()
    io = Io()
    plot = Plot()
    projection = Projection()
    simulation = Simulation()
    stats = Stats()

    def __init__(self, freqs=None, allele_names=None, id=None, parent=None,
                 ploidy=None, num_loci=None, num_alleles=None):

        # Normalize frequencies
        self._freqs = normalize_freqs_arg(freqs)
        if isinstance(freqs, np.ndarray):
            if freqs.ndim == 2:
                self.freqs_matrix = np.array(freqs)
                self.num_alleles = freqs.shape[1]
            elif freqs.ndim == 1:
                self.freqs_vector = np.array(freqs)
                self.num_alleles = 2
                self.freqs_matrix = fill_freqs_vector(self.freqs_vector)
            else:
                raise ValueError('invalid frequency data')
        if self._freqs:
            num_loci = len(self._freqs)

        # Save required attributes
        self.allele_names = allele_names
        self.id = id
        self.parent = parent
        self._last_id_index = 0

        # Save optional given lazy attributes
        if ploidy is not None:
            self.ploidy = ploidy
        if num_loci is not None:
            self.num_loci = num_loci
        if num_alleles is not None:
            self.num_alleles = num_alleles

    def __repr__(self):
        return self.io.render(id_align='best', limit=20, ind_limit=10)

    def __str__(self):
        return self.io.render(id_align='best')

    def __eq__(self, other):
        # Treat other as a sequence and compare each item.
        try:
            N = len(other)
        except TypeError:
            return NotImplemented

        if len(self) != N:
            return False
        return all(x == y for (x, y) in zip(self, other))

    def _population(self, *args, **kwargs):
        from kpop import Population
        return Population(*args, **kwargs)

    def as_array(self, which='raw'):
        """
        Convert to a numpy data array using the requested conversion method.
        This is a basic pre-processing step in many dimensionality reduction
        algorithms.

        Genotypes are categorical data and usually it doesn't make sense to
        treat the integer encoding used in kpop as ordinal data (there is
        no ordering implied when treating say, allele 1 vs allele 2 vs allele
        3).

        Conversion methods:
            * raw:
                An 3 dimensional array of (size, num_loci, ploidy) for raw
                genotype data. Each component represents the value of a single
                allele.
            * flat:
                Like raw, but flatten the last dimension into a (size,
                num_loci * ploidy) array. This creates a new feature per
                loci for each degree of ploidy in the data.
            * rflat:
                Flatten data, but first shuffle the positions of alleles at
                each loci. This is recommended if data does not carry reliable
                haplotype information.
            * raw-unity, flat-unity, rflat-unity:
                Normalized versions of "raw", "flat", and "rflat" methods. All
                components are rescaled with zero mean and unity variance.
            * count:
                Force conversion to biallelic data and counts the number of
                occurrences of the first allele. Most methdds will require
                normalization, so you probably should consider an specific
                method such as count-unity, count-snp, etc
            * count-unity:
                Normalized version of count scaled to zero mean and unity
                variance.
            * count-snp:
                Normalizes each feature using the standard deviation expected
                under the assumption of Hardy-Weinberg equilibrium. This
                procedure is described at Patterson et. al., "Population
                Structure and Eigenanalysis" and is recommended for SNPs
                subject to genetic drift.
            * count-center:
                Instead of normalizing, simply center data by subtracting half
                the ploidy to place it into a symmetric range. This
                normalization puts data into a cube with a predictable
                origin and range. For diploid data, the components will be
                either -1, 0, or 1.

        Returns:
            An ndarray with transformed data.
        """
        data = np.array([ind.data for ind in self])

        # Raw conversion
        if which == 'raw':
            return data
        elif which == 'raw-unity':
            data = data - data.mean(axis=0)
            std = data.std(axis=0)
            data /= np.where(std, std, 1)
            return data

        # Flattened representations
        elif which in {'flat', 'flat-unity'}:
            data = data.reshape(self.size, self.num_loci * self.ploidy)
            if which == 'flat-unity':
                return preprocessing.scale(data.astype(float))
            return data
        elif which in {'rflat', 'rflat-unity'}:
            return self.shuffled_loci().as_array(which[1:])

        # Counters
        elif which in {'count', 'count-unity', 'count-snp', 'count-center'}:
            count = (np.array(data) == 1).sum(axis=2)
            if which == 'count-unity':
                return preprocessing.scale(count.astype(float))
            elif which == 'count-snp':
                mu = count.mean(axis=0)
                p = mu / self.ploidy
                norm = np.sqrt(p * (1 - p))
                norm = np.where(norm, norm, 1)
                return (count - mu) / norm
            elif which == 'count-center':
                return count - self.ploidy / 2
            else:
                return count

        raise ValueError('invalid conversion method: %r' % which)

    def drop_markers(self, indexes, *, _data=None):
        """
        Create a new population with all loci in the given index removed.

        Args:
            indexes: a list of loci indexes

        Returns:
            A new population.
        """

        data = np.array([ind.data for ind in self]) if _data is None else _data
        indexes = set(indexes)
        good_loci = np.array([idx for idx in range(self.num_loci)
                              if idx not in indexes])
        data = data[:, good_loci]

        # Create new individuals
        new_inds = []
        for ind, gene in zip(self, data):
            new_inds.append(ind.copy(gene, num_alleles=2))

        return self.transformed_copy(new_inds)

    def drop_non_biallelic(self, *, _data=None):
        """
        Creates a new population that remove all non-biallelic loci.

        Returns:
            A (population, removed) tuple with the new population and a list of
            of all dropped locus indexes.
        """

        data = np.array([ind.data for ind in self]) if _data is None else _data
        bad_loci = self.stats.non_biallelic(_data=data)
        return self.drop_markers(bad_loci, _data=data), bad_loci

    def transformed_copy(self, individuals, force_multi=False, **kwargs):
        """
        Return a transformed copy of population with the given list of
        individuals. The number of individuals must match the current
        population.

        Args:
            individuals:
                A list of new, transformed individuals.

        Returns:
            The transformed Population or MultiPopulation.
        """

        # Return population or multi-population
        kwargs.setdefault('parent', self)
        kwargs.setdefault('id', self.id)
        if not self.is_multi_population:
            pop = self.population_class(individuals, **kwargs)
            if force_multi:
                return self.multi_population_class([pop], **kwargs)
            return pop
        else:
            populations = []
            idx = 0
            for pop in self.populations:
                data = individuals[idx: idx + len(pop)]
                idx += len(pop)
                pop = self._population(data, id=pop.id, parent=pop.parent)
                populations.append(pop)
            return self.multi_population_class(populations, **kwargs)

    def copy(self):
        """
        Return a copy of population.
        """
        return copy.deepcopy(self)

    def shuffled_loci(self):
        """
        Return a copy with shuffled contents of each locus.
        """

        pop = self.copy()
        for ind in pop:
            for loci in ind.data:
                np.random.shuffle(loci)
        return pop

    def _next_id(self):
        self._last_id_index += 1
        return '%s%s' % (self.id or 'ind', self._last_id_index)

import collections
import copy

import numpy as np
import pandas as pd
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
from .utils import freqs_fastest, get_freqs, set_freqs, hfreqs_vector
from ..prob import Prob
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

    # Meta information
    individual_ids = property(lambda _: _.meta['ids'])

    # Special attributes. These will be inserted later via monkey patching
    populations = ()
    admixture = Admixture()
    clustering = Clustering()
    classification = Classification()
    io = Io()
    plot = Plot()
    projection = Projection()
    simulation = Simulation()
    statistics = Stats()

    # Aliases
    admix = lazy(lambda self: self.admixture)
    cls = lazy(lambda self: self.classification)
    sim = lazy(lambda self: self.simulation)
    stats = lazy(lambda self: self.statistics)

    def __init__(self, freqs=None, allele_names=None, id=None, ploidy=None,
                 num_loci=None, num_alleles=None, individual_ids=None):

        # Normalize frequencies
        if freqs is None:
            self._freqs = None
        elif len(freqs) == 0:
            raise ValueError('cannot initialize from empty frequencies')
        elif isinstance(freqs[0], collections.Mapping):
            self._freqs = [Prob(p) for p in freqs]
        else:
            freqs = np.asarray(freqs)

            if freqs.ndim == 2:
                self._freqs = [Prob(dict(enumerate(p, 1))) for p in freqs]
                self.freqs_matrix = np.array(freqs)
                self.num_alleles = freqs.shape[1]
            elif freqs.ndim == 1:
                self._freqs = [Prob({1: p, 2: 1 - p}) for p in freqs]
                self.freqs_vector = np.array(freqs)
                self.freqs_matrix = fill_freqs_vector(self.freqs_vector)
                self.num_alleles = 2
            else:
                raise ValueError('invalid frequency data')

        # Fix num_loci from data
        if self._freqs is not None:
            self.num_loci = len(self._freqs)
            if num_loci is not None and num_loci != self.num_loci:
                raise ValueError('invalid value for num_loci')
        elif num_loci is not None:
            self.num_loci = num_loci

        # Individual ids
        if individual_ids is None:
            fmt = 'ind%s' if id is None else '%s%%s' % id
            individual_ids = [fmt % i for i in range(1, self.size + 1)]

        # Save required attributes
        self.allele_names = allele_names
        self.id = id
        self._last_id_index = 0
        self.meta = pd.DataFrame({'ids': individual_ids})

        # Save optional given lazy attributes
        if ploidy is not None:
            self.ploidy = ploidy
        if num_alleles is not None:
            self.num_alleles = num_alleles

    def __repr__(self):
        return self.io.render(id_align='best', limit=20, ind_limit=10)

    def __str__(self):
        return self.io.render(id_align='best')

    def __eq__(self, other):
        if not isinstance(other, PopulationBase):
            return NotImplemented
        if self.shape != other.shape:
            return False
        return all(x == y for x, y in zip(self, other))

    def _population(self, *args, **kwargs):
        from kpop import Population
        return Population(*args, **kwargs)

    def _next_id(self):
        self._last_id_index += 1
        return '%s%s' % (self.id or 'ind', self._last_id_index)

    def _as_array(self):
        return NotImplementedError('must be implemented on subclasses')

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
        data = self._as_array()

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
            return self.shuffle_loci().as_array(which[1:])

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

    def drop_non_biallelic(self):
        """
        Creates a new population that remove all non-biallelic loci.

        Returns:
            A (population, removed) tuple with the new population and a list of
            of all dropped locus indexes.
        """
        bad_loci = self.statistics.non_biallelic()
        return self.drop_loci(bad_loci), bad_loci

    def force_biallelic(self):
        """
        Return a new population with forced biallelic data.

        If a locus has more than 2 alleles, the most common allele is picked
        as allele 1 and the alternate allele 2 comprises all the other alleles.
        """
        alleles_mapping = [biallelic_mapping(prob) for prob in self.freqs]
        return self.map_alleles(alleles_mapping)

    def sort_by_allele_freq(self):
        """
        Return a new population in which the index attributed to each allele
        in each locus is sorted by the frequency in the population. After that,
        allele 1 will be the most common, allele 2 is the second most common
        and so on.
        """
        alleles_mapping = [sorted_allele_mapping(prob) for prob in self.freqs]
        return self.map_alleles(alleles_mapping)

    def map_alleles(self, alleles_mapping):
        """
        Create new population reorganizing all allele values by the given
        list of allele values mappings.

        Args:
            alleles_mapping:
                A list with num_loci elements. Each element must be a mapping
                from the old allele values to the new ones. If an element is
                an empty dictionary, no remapping is done.
        """
        raise NotImplementedError('must be implemented in subclasses')

    def drop_loci(self, indexes):
        """
        Create a new population with all loci in the given indexes removed.
        """
        indexes = set(indexes)
        keep = np.array([i for i in range(self.num_loci) if i not in indexes])
        return self.keep_loci(keep)

    def drop_individuals(self, indexes):
        """
        Creates new population removing the individuals in the given indexes.
        """
        indexes = set(indexes)
        keep = np.array([i for i in range(self.size) if i not in indexes])
        return self.keep_individuals(keep)

    def keep_loci(self, indexes):
        """
        Creates a new population keeping only the loci in the given indexes.
        """
        raise NotImplementedError('must be implemented in subclasses')

    def keep_individuals(self, indexes):
        """
        Creates new population removing the individuals in the given indexes.
        """
        raise NotImplementedError('must be implemented in subclasses')

    def shuffle_loci(self):
        """
        Return a copy with shuffled contents of each locus.
        """

        pop = self.copy()
        for ind in pop:
            for loci in ind.data:
                np.random.shuffle(loci)
        return pop

    def copy(self, id=None):
        """
        Return a copy of population.
        """

        new = copy.deepcopy(self)
        if id is not None:
            new.id = id
        return new


def sorted_allele_mapping(prob):
    mapping = sorted(prob.items(), key=lambda x: x[1], reverse=True)
    return {x: i for i, (x, y) in enumerate(mapping, 1) if i != x}


def biallelic_mapping(prob):
    if len(prob) <= 2:
        return {}
    else:
        idx = prob.mode()
        mapping = {i: 2 for i in prob}
        mapping[idx] = 1
        return mapping

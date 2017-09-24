import abc
import collections
import copy

from lazyutils import lazy
from sidekick import _

from .admixture import Admixture
from .classification import Classification
from .clusterization import Clusterization
from .io import Io
from .plot import Plot
from .projection import Projection
from .simulation import Simulation
from .statistics import Statistics
from .utils import discard_attrs
from .utils import get_freqs, set_freqs, hfreqs_vector, random_individual_data
from ..data_converter import DataConverter
from ..libs import kpop
from ..libs import np
from ..prob import Prob
from ..utils import fill_freqs_vector, freqs_to_matrix, fn_lazy, fn_property
from ..utils import random_frequencies


class PopulationBase(collections.Sequence, metaclass=abc.ABCMeta):
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
    dtype = lazy(lambda _: np.dtype('uint8'))
    _shape_attrs = (
        'size', 'num_loci', 'ploidy', 'shape', 'data_size',
    )

    # Frequencies
    freqs = property(get_freqs, set_freqs)
    freqs_matrix = lazy(lambda _: freqs_to_matrix(_.freqs))
    freqs_vector = lazy(lambda _: np.ascontiguousarray(_.freqs_matrix[:, 0]))
    hfreqs_vector = lazy(hfreqs_vector)

    # Allele statistics
    allele_names = None
    is_biallelic = fn_lazy(_.num_alleles == 2)
    num_alleles = lazy(lambda _: max(max(freq) for freq in _.freqs))

    # Multi population
    is_multi_population = False
    num_populations = fn_property(lambda _: len(_.populations))

    # Missing data
    has_missing_data = property(lambda _: any(ind.has_missing for ind in _))
    missing_data_total = property(
        lambda _: sum(ind.missing_data_total for ind in _))
    missing_data_ratio = fn_property(_.missing_data_total / _.data_size)

    # Meta information
    individual_ids = lazy(lambda _: list(_.meta.index))

    # Special attributes. These will be inserted later via monkey patching
    populations = ()
    admixture = Admixture()
    clusterization = Clusterization()
    classification = Classification()
    io = Io()
    plot = Plot()
    projection = Projection()
    simulation = Simulation()
    statistics = Statistics()

    # Aliases
    admix = property(lambda self: self.admixture)
    cls = property(lambda self: self.classification)
    cluster = property(lambda self: self.clusterization)
    proj = property(lambda self: self.projection)
    sim = property(lambda self: self.simulation)
    stats = property(lambda self: self.statistics)

    # List of cacheable attributes
    _cacheable_attributes = (
        'has_missing', 'missing_total', 'missing_ratio',
        'is_biallelic', 'num_alleles',
        'admixture', 'clustering', 'classification', 'io', 'plot', 'projection',
        'simulation', 'statistics',
    )

    @classmethod
    def random(cls, size=0, num_loci=0, alleles=2, ploidy=2, id=None,
               seed=None):
        """
        Creates a new random population.

        Args:
            size:
                Number of individuals. If a list of numbers is given, creates
                a Multipopulation object with sub-populations of the assigned
                sizes.
            num_loci:
                Number of loci in the genotype.
            alleles:
                Number of alleles for all loci.
            ploidy:
                Ploidy of genotype.
            min_prob:
                Minimum value for a frequency probability.

        Returns:
            A new population object.
        """
        if num_loci <= 0:
            raise ValueError('num_loci must be at least one!')

        is_multipopulation = isinstance(size, collections.Sequence)
        sizes = [size] if not is_multipopulation else size
        seeds = get_seeds(len(sizes), seed)

        # Create frequencies and data
        all_data = []
        all_freqs = [random_frequencies(num_loci, alleles, seed=k)
                     for k in seeds]
        for pre_seed, freqs, size in zip(seeds, all_freqs, sizes):
            data = []
            ind_seeds = get_seeds(size, pre_seed)
            for seed in ind_seeds:
                ind = random_individual_data(freqs, ploidy=ploidy, seed=seed)
                data.append(ind)
            all_data.append(np.array(data))

        # Return population
        if is_multipopulation:
            sub_populations = []
            for i in range(len(sizes)):
                id_i = None if id is None else '%s%s' % (id, i + 1)
                pop = kpop.Population(
                    all_data[i], freqs=all_freqs[i], id=id_i,
                    num_loci=num_loci, num_alleles=alleles, ploidy=ploidy
                )
                sub_populations.append(pop)
            return kpop.MultiPopulation(sub_populations, id=id)
        else:
            return kpop.Population(
                all_data[0], freqs=all_freqs[0], id=id,
                num_loci=num_loci, num_alleles=alleles, ploidy=ploidy
            )

    def __init__(self, freqs=None, allele_names=None, id=None, ploidy=None,
                 num_loci=None, num_alleles=None):

        # Normalize frequencies
        self._init_freqs(freqs)

        # Fix num_loci from data
        if self._freqs is not None:
            self.num_loci = len(self._freqs)
            if num_loci is not None and num_loci != self.num_loci:
                raise ValueError('invalid value for num_loci')
        elif num_loci is not None:
            self.num_loci = num_loci

        # Save required attributes
        self.allele_names = allele_names
        self.id = id

        # Save optional given lazy attributes
        if ploidy is not None:
            self.ploidy = ploidy
        if num_alleles is not None:
            self.num_alleles = num_alleles

    def _init_freqs(self, freqs):
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

    def __repr__(self):
        return self.io.render(max_loci=20, max_ind=10)

    def __str__(self):
        return self.io.render()

    def __eq__(self, other):
        if not isinstance(other, PopulationBase):
            return NotImplemented
        if self.shape != other.shape:
            return False
        return all(x == y for x, y in zip(self, other))

    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self._getitem_by_index(idx)
        elif isinstance(idx, str):
            return self._getitem_by_label(idx)
        elif isinstance(idx, slice):
            return self._getslice(idx)
        else:
            typename = idx.__class__.__name__
            raise TypeError('invalid index type: %s' % typename)

    def _getitem_by_label(self, key):
        idx = self.meta.index.get_loc(key)
        return self._getitem_by_index(idx)

    def _getitem_by_index(self, idx):
        raise NotImplementedError

    def _getslice(self, slice):
        item = self._getitem_by_index
        data = [item(i) for i in range(*slice.indices(self.size))]
        return kpop.Population(data, id=self.id)

    def _population(self, *args, **kwargs):
        from kpop import Population
        return Population(*args, **kwargs)

    def _clear_caches(self):
        discard_attrs(self, self._cacheable_attributes)

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
            * raw-norm, flat-norm, rflat-norm:
                Normalized versions of "raw", "flat", and "rflat" methods. All
                components are rescaled with zero mean and unity variance.
            * count:
                Force conversion to biallelic data and counts the number of
                occurrences of the first allele. Most methdds will require
                normalization, so you probably should consider an specific
                method such as count-unity, count-snp, etc
            * count-norm:
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
        data_converter = DataConverter(self._as_array())
        return data_converter(which)

    def drop_non_biallelic(self, **kwargs):
        """
        Creates a new population that remove all non-biallelic loci.

        Returns:
            A (population, removed) tuple with the new population and a list of
            of all dropped locus indexes.
        """
        bad_loci = self.statistics.non_biallelic()
        return self.drop_loci(bad_loci, **kwargs), bad_loci

    def force_biallelic(self, **kwargs):
        """
        Return a new population with forced biallelic data.

        If a locus has more than 2 alleles, the most common allele is picked
        as allele 1 and the alternate allele 2 comprises all the other alleles.
        """
        alleles_mapping = [biallelic_mapping(prob) for prob in self.freqs]
        return self.map_alleles(alleles_mapping, **kwargs)

    def sort_by_allele_freq(self, **kwargs):
        """
        Return a new population in which the index attributed to each allele
        in each locus is sorted by the frequency in the population. After that,
        allele 1 will be the most common, allele 2 is the second most common
        and so on.
        """
        alleles_mapping = [sorted_allele_mapping(prob) for prob in self.freqs]
        return self.map_alleles(alleles_mapping, **kwargs)

    @abc.abstractmethod
    def map_alleles(self, alleles_mapping, **kwargs):
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

    def drop_loci(self, indexes, **kwargs):
        """
        Create a new population with all loci in the given indexes removed.
        """
        indexes = set(indexes)
        keep = np.array([i for i in range(self.num_loci) if i not in indexes])
        return self.keep_loci(keep, **kwargs)

    def drop_individuals(self, indexes, **kwargs):
        """
        Creates new population removing the individuals in the given indexes.
        """
        indexes = set(indexes)
        keep = np.array([i for i in range(self.size) if i not in indexes])
        return self.keep_individuals(keep, **kwargs)

    @abc.abstractmethod
    def keep_loci(self, indexes, **kwargs):
        """
        Creates a new population keeping only the loci in the given indexes.
        """
        raise NotImplementedError('must be implemented in subclasses')

    @abc.abstractmethod
    def keep_individuals(self, indexes, **kwargs):
        """
        Creates new population removing the individuals in the given indexes.
        """
        raise NotImplementedError('must be implemented in subclasses')

    def shuffle_loci(self, **kwargs):
        """
        Return a copy with shuffled contents of each locus.
        """

        pop = self.copy(**kwargs)
        for ind in pop:
            for loci in ind.data:
                np.random.shuffle(loci)
        return pop

    def copy(self, id=None):
        """
        Return a copy of population.
        """

        new = copy.copy(self)
        new.populations = copy.copy(self.populations)
        new._clear_caches()
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


def get_seeds(n, seed, salt=0):
    """
    Return a list of seeds from the given initial seed.
    """
    if seed is None:
        return [None] * n
    else:
        np.random.seed(seed + salt)
        return list(np.random.randint(0, 2 ** 31 - 1, size=n))

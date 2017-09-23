import collections

from lazyutils import delegate_to, lazy
from sidekick import _

from kpop.population.utils import id_from_parents, random_individual_data
from .utils import parse_population_data
from ..libs import np
from ..utils import fn_lazy, fn_property

NOT_GIVEN = object()


class IndividualBase:
    """
    Base class for regular Individual and IndividualProxy objects.
    """

    # Shape
    num_loci = fn_lazy(_.data.shape[0])
    ploidy = fn_lazy(_.data.shape[1])
    dtype = delegate_to('data')
    flatten = delegate_to('data')

    # Biallelic data
    is_biallelic = fn_lazy(_.num_alleles == 2)

    # Missing data
    has_missing = fn_property(lambda _: 0 in _.data)
    missing_data_total = fn_property(lambda _: (_.data == 0).sum())

    @property
    def missing_data_ratio(self):
        return self.missing_data_total / (self.num_loci * self.ploidy)

    # Other
    _allele_names = None
    population = None
    data = None
    admixture_q = None

    # Simple queries
    is_individual = True
    is_population = False

    @lazy
    def num_alleles(self):
        if self.population:
            return self.population.num_alleles
        else:
            return self.data.max()

    @property
    def allele_names(self):
        if self._allele_names is None:
            if self.population:
                return self.population.allele_names
        return self._allele_names

    @lazy
    def admixture_vector(self):
        if self.admixture_q is None:
            return None
        values = sorted(self.admixture_q.items())
        return np.array([y for x, y in values])

    def __getitem__(self, idx):
        return self.data[idx]

    def __iter__(self):
        return iter(self.data)

    def __repr__(self):
        return 'Individual(%r)' % self.render(max_loci=20)

    def __str__(self):
        return self.render()

    def __len__(self):
        return len(self.data)

    def __eq__(self, other):
        if isinstance(other, IndividualBase):
            return (self.data == other.data).all()
        elif isinstance(other, (np.ndarray, list, tuple)):
            return (self.data == other).all()
        else:
            return NotImplemented

    def haplotypes(self):
        """
        Return a sequence of ploidy arrays with each haplotype.

        This operation is a simple transpose of genotype data.
        """

        return self.data.T

    def copy(self, data=None, *, meta=NOT_GIVEN, **kwargs):
        """
        Creates a copy of individual.
        """

        kwargs.setdefault('id', self.id)
        kwargs.setdefault('population', self.population)
        kwargs.setdefault('allele_names', self.allele_names)
        dtype = kwargs.setdefault('dtype', self.dtype)
        kwargs['meta'] = self.meta if meta is NOT_GIVEN else meta
        if data is None:
            data = np.array(self.data, dtype=dtype)
        return Individual(data, **kwargs)

    def render(self, id_align=None, max_loci=None):
        """
        Renders individual genotype.
        """

        # Choose locus rendering function
        if self.allele_names is None:
            def render_locus(idx):
                locus = self[idx]
                return ''.join(map(str, locus))
        else:
            def render_locus(idx):
                locus = self[idx]
                try:
                    mapping = self.allele_names[idx]
                except IndexError:
                    mapping = {}

                return ''.join(str(mapping.get(x, x)) for x in locus)

        size = len(self)
        id_align = -1 if id_align is None else int(id_align)
        data = [('%s:' % (self.id or 'ind')).rjust(id_align + 1)]

        # Select items
        if max_loci and size > max_loci:
            good_idx = set(range(max_loci // 2))
            good_idx.update(range(size - max_loci // 2, size))
        else:
            good_idx = set(range(size))

        # Render locus
        for i in range(len(self)):
            if i in good_idx:
                data.append(render_locus(i))

        # Add ellipsis for large data
        if max_loci and size > max_loci:
            data.insert(max_loci // 2 + 1, '...')

        return ' '.join(data)

    def render_ped(self,
                   family_id='FAM001', individual_id=0,
                   paternal_id=0, maternal_id=0,
                   sex=0, phenotype=0, memo=None):
        """
        Render individual as a line in a plink's .ped file.

        Args:
            family_id:
                A string or number representing the individual's family.
            individual_id:
                A number representing the individual's id.
            paternal_id, maternal_id:
                A number representing the individuals father/mother's id.
            sex:
                The sex (1=male, 2=female, other=unknown).
            phenotype:
                A number representing the optional phenotype.
        """

        data = '  '.join(' '.join(map(str, locus)) for locus in self.data)
        return '%s  %s  %s %s  %s  %s  %s' % (
            family_id, individual_id, paternal_id, maternal_id, sex, phenotype,
            data
        )

    def render_csv(self, sep=','):
        """
        Render individual in CSV.
        """

        data = [self.id]
        data.extend(''.join(map(str, x)) for x in self)
        return sep.join(data)

    def breed(self, other, id=None, **kwargs):
        """
        Breeds with other individual.

        Creates a new genotype in which features are selected from both
        parents.
        """

        # Haploid individuals mix parent's genome freely
        if self.ploidy == 1:
            which = np.random.randint(0, 2, size=self.num_loci)
            data1 = self.data[:, 0]
            data2 = other.data[:, 0]
            data = np.where(which, data2, data1)
            data = data[:, None]

        # Diploid create 2 segments for each parent and fuse the results
        elif self.ploidy == 2:
            which = np.random.randint(0, 2, size=self.num_loci)
            data = np.where(which, self.data[:, 0], self.data[:, 1])

            which = np.random.randint(0, 2, size=other.num_loci)
            data_other = np.where(which, other.data[:, 0], other.data[:, 1])

            data = np.stack([data, data_other], axis=1)
        else:
            raise NotImplementedError

        kwargs['id'] = id or id_from_parents(self.id, other.id)
        return self.copy(data, **kwargs)


class Individual(IndividualBase, collections.Sequence):
    """
    Represents a single individual genotype.

    A genotype data must be an integer array of shape (num_loci, ploidy).

    Args:
        data:
            Can be either a string of values or a list of raw genotype values
            represented as integers.
        population:
            Population to which individual belongs to.

    Attributes:
        num_loci:
            Number of loci in the raw genotype data
        ploidy:
            Genotype's ploidy
        data:
            A numpy array of integers with genotype data. Allele types are
            represented sequentially by 1, 2, 3, etc. Missing data is
            represented by zero. By default, data is stored in uint8 form. This
            supports up to 255 different allele types plus zero.
        allele_names:
            A list of mappings between allele integer values to a character
            representation. If not given, it inherits from parent population.
    """

    @classmethod
    def from_freqs(cls, freqs, ploidy=2, **kwargs):
        """
        Returns a random individual from the given frequency distribution.

        Args:
            freqs:
                A frequency distribution. Can be a sequence of Prob() elements
                or an square array of frequencies.
            ploidy:
                Individuals ploidy.
            **kwargs:
                Additional keyword arguments passed to the constructor.

        Returns:
            A new :class:`Individual` instance.
        """

        if isinstance(freqs[0], collections.Mapping):
            raise NotImplementedError
        else:
            data = random_individual_data(freqs, ploidy=ploidy)
        return cls(data, **kwargs)

    def __init__(self, data, id=None, population=None,
                 allele_names=None, dtype=None, meta=None, admixture_q=None,
                 num_alleles=None):

        # Convert string initial data
        if isinstance(data, str):
            data, _id = parse_population_data(data)
            if _id:
                id = id or _id[0]
            data = data[0]

        # Initialize data
        dtype = np.uint8 if dtype is None else dtype
        data = np.asarray(data, dtype=dtype)
        if len(data.shape) == 1:
            data = data[:, None]
        elif len(data.shape) > 2:
            raise ValueError('invalid data: expect an 2D array of loci')

        # Normalize allele_names
        if isinstance(allele_names, collections.Mapping):
            allele_names = dict(allele_names)
            allele_names = [allele_names for _ in data]

        # Save attributes
        self.id = id
        self.data = data
        self.population = population
        self.admixture_q = admixture_q
        self.meta = dict(meta or {})
        self._allele_names = allele_names
        self._container = None

        # Optional attributes
        if num_alleles is not None:
            self.num_alleles = num_alleles


class IndividualProxy(IndividualBase):
    """
    Individual-like instance attatched to a Population object. This class simply
    hold a reference to the population and the index of the individual in that
    population.
    """
    data = property(lambda _: _.population._data[_._idx])

    @lazy
    def id(self):
        return self.population.individual_ids[self._idx]

    @lazy
    def meta(self):
        return dict(self.population.meta.iloc[self._idx])

    def __init__(self, population, idx):
        self.population = population
        self._idx = idx

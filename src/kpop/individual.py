import collections

import numpy as np
from lazyutils import delegate_to, lazy
from sidekick import _

from kpop.kpop_parser import str_to_data
from kpop.utils import fn_lazy, fn_property

NOT_GIVEN = object()


class Individual(collections.Sequence):
    """
    Represents a single individual genotype.

    A genotype data must be an integer array of shape (num_loci, ploidy).

    Args:
        data:
            Can be either a string of values

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

    num_loci = fn_lazy(_.data.shape[0])
    ploidy = fn_lazy(_.data.shape[1])
    is_biallelic = fn_lazy(_.num_alleles == 2)
    has_missing = fn_property(lambda _: 0 in _.data)
    missing_total = fn_property(lambda _: (_.data == 0).sum())
    missing_ratio = fn_property(_.missing_total / (_.num_loci * _.ploidy))
    flatten = delegate_to('data')
    dtype = delegate_to('data')

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
            data = random_data_from_freqs_matrix(freqs, ploidy=ploidy)
        return cls(data, **kwargs)

    def __init__(self, data, copy=True, label=None, population=None,
                 allele_names=None, dtype=np.uint8, meta=None, admixture_q=None,
                 num_alleles=None):

        # Convert string initial data
        if isinstance(data, str):
            if ':' in data:
                prefix, data = data.split(':')
                if label is None:
                    label = prefix

            data = data.strip()
            data, allele_names_ = str_to_data(data, allele_names)
            copy = False
            if population is None and allele_names is None:
                allele_names = allele_names_
            if population and population.allele_names is None:
                population.allele_names = allele_names_

        # Initialize data
        data = (np.array(data, copy=True, dtype=dtype) if copy else
                np.asarray(data, dtype=dtype))
        if len(data.shape) == 1:
            data = data[:, None]
        elif len(data.shape) > 2:
            raise ValueError('invalid data: expect an 2D array of loci')

        # Normalize allele_names
        if isinstance(allele_names, collections.Mapping):
            allele_names = dict(allele_names)
            allele_names = [allele_names for _ in data]

        # Save attributes
        self.data = data
        self._allele_names = allele_names
        self.population = population
        self.admixture_q = admixture_q
        self._container = None
        self.label = label
        self.meta = dict(meta or {})

        # Save lazy overrides
        if num_alleles is not None:
            self.num_alleles = num_alleles

    def __getitem__(self, idx):
        return self.data[idx]

    def __iter__(self):
        return iter(self.data)

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.render(limit=20))

    def __str__(self):
        return self.render()

    def __len__(self):
        return len(self.data)

    def haplotypes(self):
        """
        Return a sequence of ploidy arrays with each haplotype.

        This operation is a simple transpose of genotype data.
        """

        return self.data.T

    def shuffled_loci(self):
        """
        Randomize the position of values in same locus *inplace*.
        """

        new = self.copy()
        for loci in new.data:
            np.random.shuffle(loci)
        return new
    
    #
    # Coping, saving and serialization
    #
    def copy(self, data=None, *, copy=True, meta=NOT_GIVEN, **kwargs):
        """
        Creates a copy of individual.
        """

        kwargs.setdefault('label', self.label)
        kwargs.setdefault('population', self.population)
        kwargs.setdefault('allele_names', self.allele_names)
        dtype = kwargs.setdefault('dtype', self.dtype)
        kwargs['meta'] = self.meta if meta is NOT_GIVEN else meta
        kwargs['copy'] = copy
        if data is None:
            data = np.array(self.data, copy=copy, dtype=dtype)
        return Individual(data, **kwargs)

    def render(self, label_align=None, limit=None):
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
        label_align = -1 if label_align is None else label_align
        data = [((self.label or 'ind') + ':').rjust(label_align + 1)]

        # Select items
        if limit and size > limit:
            good_idx = set(range(limit // 2))
            good_idx.update(range(size - limit // 2, size))
        else:
            good_idx = set(range(size))

        # Render locus
        for i in range(len(self)):
            if i in good_idx:
                data.append(render_locus(i))

        # Add ellipsis for large data
        if limit and size > limit:
            data.insert(limit // 2 + 1, '...')

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

        data = [self.label]
        data.extend(''.join(map(str, x)) for x in self)
        return sep.join(data)

    def breed(self, other, label=None, **kwargs):
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

        # Diploid create 2 segments for each parent and fuse the _results
        elif self.ploidy == 2:
            which = np.random.randint(0, 2, size=self.num_loci)
            data = np.where(which, self.data[:, 0], self.data[:, 1])

            which = np.random.randint(0, 2, size=other.num_loci)
            data_other = np.where(which, other.data[:, 0], other.data[:, 1])

            data = np.stack([data, data_other], axis=1)
        else:
            raise NotImplementedError

        kwargs['label'] = label or label_from_parents(self.label, other.label)
        return self.copy(data, copy=False, **kwargs)


def label_from_parents(l1, l2):
    """
    Creates a new label from parents labels.
    """

    if l1 == l2:
        return l1
    elif l1 is None:
        return l2 + '_'
    elif l2 is None:
        return l1 + '_'

    common = []
    for c1, c2 in zip(l1, l2):
        if c1 == c2:
            common.append(c1)
            continue
        break
    common = ''.join(common)
    l1 = l1[len(common):] or '_'
    l2 = l2[len(common):] or '_'
    return '%s-%s,%s' % (common, l1, l2)


def random_data_from_freqs_matrix(freqs, ploidy=2):
    """
    Creates a random biallelic individual with given ploidy.

    Freqs can be a list of (p, 1 - p) pairs or a list of p's.
    """

    freq_array = np.asarray(freqs)
    num_loci = len(freq_array)
    values = np.random.uniform(size=num_loci * ploidy)
    values = values.reshape((num_loci, ploidy))
    if freq_array.ndim == 2:
        if freq_array.shape[1] != 2:
            raise NotImplementedError('must be biallelic!')
        return np.asarray(values >= freq_array[:, 0, None], dtype=np.uint8) + 1
    else:
        return np.asarray(values >= freq_array[:, None], dtype=np.uint8) + 1

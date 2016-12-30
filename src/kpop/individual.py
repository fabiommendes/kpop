import collections
from random import choice

import numpy as np
from lazyutils import lazy, delegate_to


class Individual(collections.Sequence):
    """
    Represents a single individual genotype.

    A genotype data must be an integer array of shape (num_loci, ploidy).

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

    @property
    def num_loci(self):
        return self._data.shape[0]

    @property
    def ploidy(self):
        return self._data.shape[1]

    @property
    def data(self):
        return self._data

    @property
    def allele_names(self):
        if self._allele_names is None:
            if self.population:
                return self.population.allele_names
        return self._allele_names

    @property
    def has_missing(self):
        return 0 in self._data

    @property
    def missing(self):
        return self._data == 0

    @property
    def missing_pc(self):
        return sum(self.missing) / (self.num_loci * self.ploidy)

    # Inherited attributes
    flatten = delegate_to('_data')
    dtype = delegate_to('_data')

    @classmethod
    def random(cls, data, ploidy=2, freq=None, alleles=None, **kwargs):
        """
        Create a new random individual.

        Args:
            data:
                The first argument can be interpreted in several ways depending
                on the type.

                int:
                    Number of loci in population. It accepts optional keyword
                    arguments ``freqs`` that maps alleles to frequencies.
                sequence:
                    Interpreted as a sequence of distributions for each allele.
        """

        raise NotImplementedError

    @classmethod
    def random_biallelic(cls, data, ploidy=2, **kwargs):
        """
        Creates random biallelic individual. Alleles values are 0 and 1.

        Args:
            data:
                If integer, describes the number of loci. Alleles are chosen
                with equal probability.

                Can also provides a list of probabilities for the first allele
                of value 0.
            ploidy:
                Individual ploidy.
        """

        if isinstance(data, int):
            data = [0.5] * data
        data = random_ind_from_freqs(data, ploidy=ploidy)
        return cls(data, copy=False, **kwargs)

    def __init__(self, data, copy=True, label=None, population=None,
                allele_names=None, dtype=np.uint8):

        # Convert string initial data
        if isinstance(data, str):
            data, allele_names_ = str_to_data(data, allele_names)
            copy = False
            if population is None and allele_names is None:
                allele_names = allele_names_
            if population and population.allele_names is None:
                population.allele_names = allele_names_

        # Initialize data
        if copy:
            data = np.array(data, copy=True, dtype=dtype)
        else:
            data = np.asarray(data, dtype=dtype)
        if len(data.shape) != 2:
            raise ValueError('invalid data: expect an 2D array of loci')

        # Normalize allele_names
        if isinstance(allele_names, collections.Mapping):
            allele_names = dict(allele_names)
            allele_names = [allele_names for _ in data]

        # Save attributes
        self._data = data
        self._allele_names = allele_names
        self.population = population
        self.label = label

    def __getitem__(self, idx):
        return self._data[idx]

    def __iter__(self):
        return iter(self._data)

    def __str__(self):
        return self.render(label_align=6)

    def __len__(self):
        return len(self._data)

    def render(self, label_align=None):
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
        if size > 20:
            good_idx = set(range(10))
            good_idx.update(range(size - 10, size))
        else:
            good_idx = set(range(size))

        # Render locus
        for i in range(len(self)):
            if i in good_idx:
                data.append(render_locus(i))

        # Add ellipsis for large data
        if size > 20:
            data.insert(11, '...')

        return ' '.join(data)

    def breed(self, other, label=None):
        """
        Breeds with other individual.

        Creates a new genotype in which features selected from both parents.
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

        return Individual(data, copy=False,
                          label=label,
                          population=self.population,
                          allele_names=self.allele_names)

    def chromosomes(self):
        """
        Return a sequence of ploidy arrays with each chromosome genotype.

        This operation is a simple transpose of genotype data.
        """

        return self._data.T

    def shuffle(self):
        """
        Randomize the position of alleles at same locus *inplace*.
        """

        for loci in self.data:
            np.random.shuffle(loci)


def check_valid_names(names):
    """
    Raises ValueError if names mapping contains invalid values.
    """

    if 0 in names or '0' in names:
        raise ValueError('missing data should be represented by a dash')


def tokenize_locus(st):
    """
    Split string with genotypes into its components.
    """

    return st.split(',') if ',' in st else list(st)


def max_non_numeric_value(names_map):
    return max(v for k, v in names_map.items() if not k.isdigit())


def update_names_map(names_map, tokens):
    """
    Update names mapping with given tokens.
    """

    names_map['-'] = 0
    if '0' in tokens:
        raise ValueError(
            'missing data must be represented by a dash "-". Do not use zero '
            'directly.')

    # Update mapping with new tokens.
    if not all(tk in names_map for tk in tokens):
        numeric = (tk for tk in tokens if tk.isdigit())
        names_map.update({tk: int(tk) for tk in numeric})
        non_numeric = sorted(tk for tk in tokens
                             if not tk.isdigit() and not tk in names_map)
        value = max_non_numeric_value(names_map) + 1
        for tk in non_numeric:
            names_map[tk] = value
            value += 1


def str_to_data(data, names=None):
    """
    Convert string of the form "aA aa AA" to a tuple of
    (genotype_data, allele_names).

    Args:
        data:
            Input string data.
        names:
            Optional list of mappings between allele name to its numerical
            value. If a single dictionary is given, it assumes that it uses a
            common mapping to all alleles.
    """
    genotypes = data.split()

    # Normalize names
    if names is None:
        names = [{} for _ in genotypes]
    elif isinstance(names, collections.Mapping):
        mapping = dict(names)
        names = [mapping for _ in genotypes]
    else:
        names = list(names)
        if len(names) < len(genotypes):
            names.extend([{} for _ in range(len(genotypes) - len(names))])
    for x in names:
        check_valid_names(x)

    # Save data
    result = []
    for names_map, genotype in zip(names, genotypes):
        tokens = tokenize_locus(genotype)
        update_names_map(names_map, tokens)
        locus = [names_map[tk] for tk in tokens]
        result.append(locus)
    return np.array(result), [{v: k for k, v in x.items()} for x in names]


def random_ind_from_freqs(freqs, ploidy=2):
    """
    Creates a random biallelic individual with given ploidy.
    """

    freq_array = np.asarray(freqs)
    N = len(freq_array)
    R = np.random.uniform(size=N * ploidy).reshape((N, ploidy))
    if freq_array.ndim == 2:
        return np.asarray(R < freq_array[:, 0, None], dtype=np.int8)
    else:
        return np.asarray(R < freq_array[:, None], dtype=np.int8)


data, map = str_to_data('11 ab aA ba')
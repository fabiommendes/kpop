import collections
import pickle
from random import randrange

import numpy as np
from lazyutils import lazy

import kpop
from kpop import Individual
from kpop.freqmatrix import fill_frequencies
from kpop.population.util import normalize_freqs_arg
from kpop.prob import Prob


class PopulationBase(collections.Sequence):
    """
    Base class for Population and Multipopulation.
    """

    @property
    def _freqs_fastest(self):
        try:
            return self.freqs_vector
        except ValueError:
            pass
        if self.num_alleles <= 5:
            return self.freqs_matrix
        return self.freqs

    @property
    def freqs(self):
        if self._freqs is None:
            self._freqs = self.empirical_freqs()
        return self._freqs

    @freqs.setter
    def freqs(self, value):
        self._freqs = normalize_freqs_arg(value)
        del self.freqs_vector
        del self.freqs_matrix

    @lazy
    def freqs_matrix(self):
        return freqs_to_matrix(self.freqs)

    @lazy
    def freqs_vector(self):
        if not self.is_biallelic:
            raise ValueError(
                'only biallelic individuals define a freqs_vector')
        return np.ascontiguousarray(self.freqs_matrix[:, 0])

    @lazy
    def is_biallelic(self):
        return self.num_alleles == 2

    @lazy
    def num_alleles(self):
        return max(max(freq) for freq in self.freqs)

    @lazy
    def ploidy(self):
        return self[0].ploidy

    @lazy
    def num_loci(self):
        return self[0].num_loci

    @property
    def size(self):
        return len(self)

    @property
    def plot(self):
        from kpop.population.plots import PlotAttribute
        return PlotAttribute(self)

    @property
    def has_missing(self):
        return any(ind.has_missing for ind in self)

    @property
    def missing_ratio(self):
        a, b, c = self.shape
        return self.missing_total / (a * b * c)

    @property
    def missing_total(self):
        return sum(ind.missing_total for ind in self)

    @property
    def num_populations(self):
        return len(self.populations)

    @property
    def shape(self):
        return self.size, self.num_loci, self.ploidy

    class population_class:

        def __get__(self, obj, cls=None):
            cls = PopulationBase.population_class = kpop.Population
            return cls

    population_class = population_class()

    class multi_population_class:

        def __get__(self, obj, cls=None):
            cls = PopulationBase.multi_population_class = kpop.MultiPopulation
            return cls

    multi_population_class = multi_population_class()

    is_multi_population = False

    def __init__(self, freqs=None, allele_names=None, label=None, parent=None,
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
                self.freqs_matrix = fill_frequencies(self.freqs_vector)
            else:
                raise ValueError('invalid frequency data')
        if self._freqs:
            num_loci = len(self._freqs)

        # Save required attributes
        self.allele_names = allele_names
        self.label = label
        self.parent = parent
        self._last_label_index = 0

        # Save optional given lazy attributes
        if ploidy is not None:
            self.ploidy = ploidy
        if num_loci is not None:
            self.num_loci = num_loci
        if num_alleles is not None:
            self.num_alleles = num_alleles

    def __repr__(self):
        return self.render(label_align='best', limit=20, ind_limit=10)

    def __str__(self):
        return self.render(label_align='best')

    def compact_data(self, copy=False):
        """
        Return all genotype data as a contiguous array.

        Args:
            copy:
                If True, forces to return a copy, otherwise, it tries to recycle
                the array if data is already compactified.

        Returns:
            A (size, num_loci, ploidy) array.
        """

        return np.array([ind.data for ind in self])

    def allele_count(self, allele=1):
        """
        Return an array of (size, num_loci) with the counts for the number of
        times the given allele appears in each individual at each locus.
        """

        data = self.compact_data()
        return (data == allele).sum(axis=2)

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
            names = ['loci%s' % i for i in range(1, self.size + 1)]
        header = 'label,' + ','.join(names)
        data = '\n'.join(x.render_csv(**kwargs) for x in self)
        return '%s\n%s' % (header, data)

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
            data.append('1 snp%s 0 %s' % (j, j))
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
            raise ValueError('invalid file format: %r' % format)

    @classmethod
    def load(cls, file, format='pickle', **kwargs):
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

        counters = [collections.Counter() for _ in range(self.num_loci)]
        for ind in self:
            for counter, genotype in zip(counters, ind):
                counter.update([x for x in genotype if x])
        freqs = [Prob({k: v + alpha for k, v in c.items()}) for c in counters]

        # Normalize return value
        if as_matrix:
            return freqs_to_matrix(freqs, self.num_alleles)
        else:
            return freqs

    def genetic_drift(self, n_generations, population_size=None,
                      sample_size=None, label=None):
        """
        Applies a simple genetic drift model for mutating the frequencies.
        Returns a new empty population initialized with the new frequencies.

        Args:
            n_generations:
                Number of generations to apply the drift model
            population_size:
                Effective size of population. Defaults to the current population
                size.
            sample_size:
                Number of individuals of the newly initialized population.
                Defaults to the current population size.
            label:
                Name of new population.

        Example:
            >>> pop = Population.make_random(20, num_loci=50)
            >>> pop2 = pop.genetic_drift(20, population_size=10)

            After many generations, genetic drift in a small population tends
            to fixate alleles. The above case has a chance greater than 50%
            of fixating each allele.

            >>> freqs_binomial = pop2.freqs_vector
            >>> fixed = freqs_binomial * (1 - freqs_binomial) == 0
            >>> fixed.any() # it is random, but we can be pretty sure that at
            ...             # least 1 allele has been fixed!
            True

        """

        population_size = population_size or self.size
        eff_size = self.ploidy * population_size
        label = label or _sub_population_label(self.label, n_generations)
        if sample_size is None:
            sample_size = self.size

        if not self.is_biallelic:
            raise NotImplementedError
        else:
            freqs = self.freqs_vector
            for _ in range(n_generations):
                for j, f in enumerate(freqs):
                    freqs[j] = np.random.binomial(eff_size, f) / eff_size
            freqs = fill_frequencies(freqs)

        # Create population
        from kpop import Population

        pop = Population(freqs=freqs, label=label)
        pop.num_loci = self.num_loci
        pop.num_alleles = self.num_alleles
        pop.ploidy = self.ploidy
        pop.freqs_matrix = freqs

        # Fill with individuals
        if sample_size:
            pop.fill(sample_size)
        return pop

    def new(self, label=None, **kwargs):
        """
        Return a new random individual respecting the population allele
        frequencies.
        """

        kwargs['label'] = label or self._next_label()
        ind = Individual.from_freqs(self._freqs_fastest,
                                    population=self, **kwargs)
        ind.num_alleles = self.num_alleles
        return ind

    def new_offspring(self, i=None, j=None, **kwargs):
        """
        Return a new offspring created by breeding two random elements in the
        population.
        """

        size = len(self)
        if i is None:
            i = randrange(size)
        if j is None:
            j = randrange(size)
        return self[i].breed(self[j], population=self, **kwargs)

    def random(self):
        """
        Return a random individual from population.
        """

        i = randrange(len(self))
        return self[i]

    def breed(self, population, size, label='breed', parent=None, **kwargs):
        """
        Breed individuals from both populations.

        Args:
            population:
                Population to choose individuals from.
            size:
                Number of offspring.
            label:
                Label for the resulting population.
            **kwargs:

        Returns:
            A new population.
        """

        children = []
        for i in range(size):
            father, mother = self.random(), population.random()
            child = father.breed(mother, label='%s%s' % (label, i), **kwargs)
            children.append(child)
        return kpop.Population(children, parent=parent, label=label)

    def structure(self):
        """
        Not implemented
        """

    def structure_kmeans(self, k):
        """
        Not implemented
        """

    def structure_admixture(self, k, pop_labels=None, parental_labels=None):
        """
        Runs the ADMIXTURE program to detect structure in the population.

        Args:
            k: number of parental populations.

        Returns:
            A new Population object with all individuals classified with their
            respective admixture coefficients.
        """

        if not self.is_biallelic:
            raise ValueError('ADMIXTURE only supports biallelic populations')

        from kpop.external.admixture2 import run_admixture

        kwargs = {}
        if pop_labels:
            kwargs['supervised'] = pop_labels
        result = run_admixture(self, k, disp=0, **kwargs)
        parental = result.make_parental(labels=parental_labels)
        individuals = result.make_admixture_labels(self)
        out = self.transformed_copy(individuals, parent=parental)
        out.admixture_result = result
        return out

    def pca(self, k=2, *, method='flatten', norm=True):
        """
        Principal component analysis.

        PCA helps classifying individuals by reducing the problem's
        dimensionality from the number of loci to only the first K of principal
        axis. Usually this simple re-parametrization is able to detect population
        structure between data, i.e., the individuals of each sub-population are
        very well separated in the "principal axes" coordinates.

        The principal axes are chosen from a singular value decomposition (SVD)
        of the population matrix. This can be done in 2 different ways:

        1) The 'count' method simply counts the number of #1 alleles in each
           locus and creates a matrix C(n,j) of the number of #1 alleles for
           individual n in location j. This only works with biallelic data.
        2) The 'flatten' method shuffle the alleles at each loci, flatten, and
           creates a matrix C(n, j) with (N x 2J) elements. The rest of the
           algorithm proceeds identically.

        The matrix C is subjected to a SVD, and each individual is classified in
        terms of their coordinates in the K most significant principal axis.

        Args:
            k (int):
                Number of principal components to consider
            method : 'count', 'flatten'
                Genetic discrimination method. See description above.
            norm : bool
                If True (default), normalizes each j component of the C(n, j)
                matrix according to a measure of genetic drift. This procedure
                is described at Patterson et. al., "Population Structure and
                Eigenanalysis" and is recommended for SNPs subject to genetic
                drift. It is not recommended to use this normalization to
                micro-satellite data.

        Returns:
            A (size x k) matrix with the components of each individual in the
            principal axis.

        Examples:
            Consider a random synthetic population with two sub-populations with
            10 elements each. Each individual has 200 alleles.

            >>> popA = Population.make_random(10, 200, label='A')
            >>> popB = Population.make_random(10, 200, label='B')

            Usually the the principal axis alone will be enough to classify
            these individuals. Since the mean is zero, individuals of one
            population will have negative coordinates and the other population
            will have positive coordinates.

            >>> pop = popA + popB
            >>> coords = pop.pca()
            >>> sign_p1 = coords[10:] > 0
            >>> sign_p2 = coords[:10] > 0
            >>> (sign_p1 != sign_p2).all()
            True
        """

        # Compute covariance matrix for each feature
        pop = np.array(self)
        if method == 'count':
            cov_matrix = (pop == 1).sum(axis=2)
        elif method == 'flatten':
            pop = np.copy(pop)
            for ind in pop:
                for locus in ind:
                    np.random.shuffle(locus)
            cov_matrix = np.array([ind.flatten() for ind in pop],
                                  dtype=pop.dtype)
        else:
            raise ValueError(method)

        # Remove bias and re-normalize, if desired
        mu = cov_matrix.mean(axis=0)
        if norm:
            p = mu / 2
            norm = np.sqrt(p * (1 - p))
            # prevents problems with mean 0 or 1
            norm = np.where(norm, norm, 1)
            matrix = (cov_matrix - mu) / norm
        else:
            matrix = cov_matrix - mu

        # Compute the singular value decomposition of the rescaled matrix and
        # return the projections of each individual in this matrix
        _, _, pca_coords = np.linalg.svd(matrix, False)
        pca_coords_k = pca_coords[:k]
        return np.dot(cov_matrix, pca_coords_k.T)

    def non_biallelic(self, *, _data=None):
        """
        Return a list of all non biallelic loci.
        """

        data = np.array([ind.data for ind in self]) if _data is None else _data
        num_alleles = data.max(axis=(0, 2))
        return np.where(num_alleles > 2)[0]

    def drop_non_biallelic(self, *, _data=None):
        """
        Creates a new population that remove all non-biallelic loci.

        Returns:
            A (population, removed) tuple with the new population and a list of
            of all dropped locus indexes.
        """

        data = np.array([ind.data for ind in self]) if _data is None else _data
        bad_loci = self.non_biallelic(_data=data)
        return self.drop_markers(bad_loci, _data=data), bad_loci

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
        kwargs.setdefault('label', self.label)
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
                pop = kpop.Population(data, label=pop.label, parent=pop.parent)
                populations.append(pop)
            return self.multi_population_class(populations, **kwargs)

    def shuffle_loci(self):
        """
        Shuffle contents of each locus for each individual *inplace*.
        """

        for ind in self:
            ind.shuffle_loci()

    def _next_label(self):
        self._last_label_index += 1
        return '%s%s' % (self.label or 'ind', self._last_label_index)


def freqs_to_matrix(freqs, num_alleles=None):
    """
    Convert a list of Probs() into a frequency matrix.
    """

    num_loci = len(freqs)
    num_alleles = num_alleles or max(max(freq) for freq in freqs)
    data = np.zeros((num_loci, num_alleles), dtype=float)
    for i, distr in enumerate(freqs):
        for k, p in distr.items():
            data[i, k - 1] = p
    return data


def _sub_population_label(label, n_gen):
    """Produces a label for a sub-population."""

    if label is None:
        return None
    else:
        return '%s-gen%s' % (label, n_gen)

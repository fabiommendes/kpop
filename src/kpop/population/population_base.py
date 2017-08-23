import collections
import copy
import pickle
import operator as op
from functools import reduce
from random import randrange, shuffle

import numpy as np
from lazyutils import lazy
from sidekick import _, F

import kpop
from .io import RenderableMixin
from .utils import normalize_freqs_arg
from ..evolution import frequency_drift
from ..individual import Individual
from ..prob import Prob
from ..statistics import biallelic_pairwise_fst
from ..utils.frequencies import fill_freqs_vector, freqs_to_matrix
from ..utils.functional import fn_lazy, fn_property


class PopulationBase(RenderableMixin, collections.Sequence):
    """
    Base class for Population and MultiPopulation.

    Attrs:
        freqs:
            ??
        freqs_matrix:
            ??
        freqs_vector:
            Frequencies for allele 1. This is more useful for biallelic data, 
            since the frequency of the second allele is simply the complement.
        hfreqs_vector:
            Vector of frequencies of heterozygotes.

        
    """

    freqs_matrix = lazy(lambda _: freqs_to_matrix(_.freqs))
    freqs_vector = lazy(lambda _: np.ascontiguousarray(_.freqs_matrix[:, 0])) 

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
            vars_ = self.__dict__
            data = vars_.get('freqs_matrix') or vars_.get('freqs_vector')

            if data is not None:
                if data.ndim == 1:
                    data = freqs_to_matrix(data)
                self._freqs = [Prob(enumerate(locus)) for locus in data]

            else:
                self._freqs = self.empirical_freqs()

        return self._freqs

    @freqs.setter
    def freqs(self, value):
        self._freqs = normalize_freqs_arg(value)
        del self.freqs_vector
        del self.freqs_matrix

    @lazy
    def hfreqs_vector(self):
        if self.ploidy == 2:
            data = np.array(self)
            heterozygote = data[:, :, 0] != data[:, :, 1]
            return heterozygote.mean(axis=0)
        else:
            raise NotImplementedError('require diploid populations')

    # General shape
    size = property(len)
    num_loci = lazy(lambda _: _[0].num_loci)
    ploidy = lazy(lambda _: _[0].ploidy)
    shape = property(lambda _: (_.size, _.num_loci, _.ploidy))
    data_size = fn_property(_.size * _.num_loci * _.ploidy)
    dtype = np.dtype('uint8')

    # Allele statistics
    is_biallelic = fn_lazy(_.num_alleles == 2)
    num_alleles = lazy(lambda _: max(max(freq) for freq in _.freqs))

    # Multi population
    is_multi_population = False
    num_populations = fn_property(lambda _: len(_.populations))

    # Missing data
    has_missing = property(lambda _: any(ind.has_missing for ind in _))
    missing_total = property(lambda _: sum(ind.missing_total for ind in _))
    missing_ratio = fn_property(_.missing_total / _.data_size)

    #
    # Special attributes. These will be inserted later via monkey patching
    #
    populations = ()
    plot = object()
    clustering = object()
    classification = object()
    admixture = object()
    projection = object()
    stats = object()
    io = object()

    #
    # Population constructors
    #
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
                self.freqs_matrix = fill_freqs_vector(self.freqs_vector)
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

    def __eq__(self, other):
        # Treat other as a sequence and compare each item.
        try:
            N = len(other)
        except TypeError:
            return NotImplemented

        if len(self) != N:
            return False
        return all(x == y for (x, y) in zip(self, other))

    def as_array(self, copy=False):
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

    #
    # Statistics
    #
    def allele_count(self, allele=1):
        """
        Return an array of (size, num_loci) with the counts for the number of
        times the given allele appears in each individual at each locus.
        """

        data = self.as_array()
        return (data == allele).sum(axis=2)

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

    def non_biallelic(self, *, _data=None):
        """
        Return a list of all non biallelic loci.
        """

        data = np.array([ind.data for ind in self]) if _data is None else _data
        num_alleles = data.max(axis=(0, 2))
        return np.where(num_alleles > 2)[0]

    def fst_matrix(self, sizes=None):
        """
        Return the Fst divergence matrix for all sub-populations.

        See Also:
            :func:`kpop.statistics.biallelic_pairwise_fst`
        """

        pops = self.populations
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
        names = [pop.label for pop in self.populations]
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

        loci_names = ['L%s' % (i + 1) for i in range(self.num_loci)]
        columns = [pop.label for pop in self.populations]
        columns.insert(0, 'locus')
        freqs = [pop.freqs_vector for pop in self.populations]
        data = list(zip(loci_names, *freqs))

        # Call csv_lines with the prepared data
        kwargs['align'] = align
        kwargs['sep'] = sep
        kwargs['decimal_places'] = decimal_places
        lines = csv_lines(data, columns=columns, **kwargs)
        return ''.join(lines)

    #
    # Simulation
    #
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
            >>> pop = Population.random(20, num_loci=50)
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
        sample_size = self.size if sample_size is None else sample_size

        # Compute drift
        freqs = frequency_drift(self.freqs_matrix, n_generations, eff_size)

        # Create population
        from kpop import Population

        pop = Population(freqs=freqs, label=label)
        pop.num_loci = self.num_loci
        pop.num_alleles = self.num_alleles
        pop.ploidy = self.ploidy
        pop.freqs_matrix = freqs
        pop.freqs_vector = freqs[:, 0]

        # Fill with individuals
        if sample_size:
            pop.fill(sample_size)
        return pop

    #
    # New individuals
    #
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
            father, mother = self.random_individual(), population.random_individual()
            child = father.breed(mother, label='%s%s' % (label, i), **kwargs)
            children.append(child)
        return kpop.Population(children, parent=parent, label=label)

    def new_individual(self, label=None, **kwargs):
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

        This individual is not added to the population. If you want that,
        just do:

        >>> pop.add(pop.new_offspring())                        # doctest: +SKIP
        """

        size = len(self)
        if i is None:
            i = randrange(size)
        if j is None:
            j = randrange(size)
        return self[i].breed(self[j], population=self, **kwargs)

    def random_individual(self):
        """
        Return a random individual from population.
        """

        i = randrange(len(self))
        return self[i]

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

        from kpop.external.admixture.base import run_admixture

        kwargs = {}
        if pop_labels:
            kwargs['supervised'] = pop_labels
        result = run_admixture(self, k, disp=0, **kwargs)
        parental = result.make_parental(labels=parental_labels)
        individuals = result.make_admixture_labels(self)
        out = self.transformed_copy(individuals, parent=parental)
        out.admixture_result = result
        return out


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

    def shuffled_loci(self):
        """
        Return a copy with shuffled contents of each locus.
        """
        
        pop = self.copy()
        for ind in pop:
            for loci in ind.data:
                np.random.shuffle(loci)
        return pop

    def copy(self):
        """
        Return a copy of population.
        """
        return copy.deepcopy(self)

    def _next_label(self):
        self._last_label_index += 1
        return '%s%s' % (self.label or 'ind', self._last_label_index)


def _sub_population_label(label, n_gen):
    """Produces a label for a sub-population."""

    if label is None:
        return None
    else:
        return '%s-gen%s' % (label, n_gen)

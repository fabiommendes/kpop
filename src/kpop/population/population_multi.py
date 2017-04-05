import numpy as np
from lazyutils import lazy

import kpop
from kpop import admixture

from kpop.admixture import likelihood
from kpop.population.population_base import PopulationBase
from kpop.population.population_single import Population
from kpop.population.attr_populations import PopulationsList
from kpop.prob import Prob


class MultiPopulation(PopulationBase):
    """
    A population formed by several sub-populations.
    """

    is_multi_population = True

    @lazy
    def freqs_vector_pop(self):
        return np.array([pop.freqs_vector for pop in self.populations])

    @lazy
    def freqs_matrix_pop(self):
        return np.array([pop.freqs_matrix for pop in self.populations])

    @lazy
    def freqs_pop(self):
        return [pop.freqs for pop in self.populations]

    def prior(self):
        """
        Prior probability for each sub-population.
        """

        size = len(self.populations)
        if size == 0:
            return np.array([], dtype=float)
        return np.zeros(size, dtype=float) + 1 / size

    def __init__(self, populations=None, freqs=None, **kwargs):
        if freqs is not None:
            raise ValueError('cannot specify frequencies on MultiPopulation '
                             'objects')
        self.populations = PopulationsList()
        if populations:
            for population in populations:
                self.add_population(population)
        super().__init__(**kwargs)

    def __len__(self):
        return sum(len(x) for x in self.populations)

    def __getitem__(self, idx):
        if isinstance(idx, int):
            i = idx
            for pop in self.populations:
                size = len(pop)
                if i >= size:
                    i -= size
                else:
                    return pop[i]
            raise IndexError(idx)

    def __iter__(self):
        for pop in self.populations:
            yield from pop

    def __add__(self, other):
        self_populations = self.populations
        if isinstance(other, Population):
            populations = list(self_populations)
            if other not in populations:
                populations.append(other)
            return MultiPopulation(populations)
        elif isinstance(other, MultiPopulation):
            populations = list(self_populations)
            populations.extend([x for x in other.populations
                                if x not in populations])
            return MultiPopulation(populations)
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Population):
            populations = [other]
            populations.extend([x for x in self.populations if x is not other])
            return MultiPopulation(populations)
        return NotImplemented

    def add_population(self, population):
        """
        Adds a new sub-population.

        Args:
            population: A :class:`Population` instance.
        """

        if population in self.populations:
            raise ValueError('sub-population already present!')
        self.populations.append(population)

    def classify(self, ind, prior=None):
        """
        Classify individual in one of the parent populations.

        Args:
            ind: a :class:`Individual` instance.

        Returns:
            The label or index (if no label is defined) of the selected
            population.
        """

        return self.prob_classify(ind, prior=prior).mode()

    def prob_classify(self, ind, prior=None):
        """
        Classify individual in one of parental populations and return a
        probability distribution over population labels.

        Args:
            ind: a :class:`Individual` instance.

        Returns:
            A :class:`Prob` mapping from population labels to probabilities.
        """

        if prior is None:
            prior = self.prior()

        if self.is_biallelic:
            probs = likelihood.bayes(ind.data, self.freqs_vector_pop,
                                     prior=prior)
            data = {}
            for i, pop in enumerate(self.populations):
                data[pop.label or i] = probs[i]
            return Prob(data)
        else:
            raise NotImplementedError

    def admixture(self, ind, method='maxlike'):
        """
        Compute the admixture coefficients for individual

        Args:
            ind: a :class:`Individual` instance.

        Returns:
            A :class:`Prob` mapping from population labels to admixture
            coefficients.
        """

        if self.is_biallelic:
            f_pops = self.freqs_matrix_pop
            result = admixture.admixture(ind.data, f_pops, method=method)
            data = {}
            for i, pop in enumerate(self.populations):
                data[pop.label or i] = result[i]
            return Prob(data)
        else:
            raise NotImplementedError

    def admixture_classify(self, pop, classifier='admixture',
                           keep_parental=False, **kwargs):
        """
        Compute admixture coefficients for individuals in the given population
        by considering that the current MultiPopulation object represents the
        known parental populations.

        Args:
            pop:
                A population or multi-population that should be classified.
            classifier:
                The method used to classify individuals:
                    1) 'admixture' -- uses the ADMIXTURE program.
                    2) 'structure' -- uses the STRUCTURE program.
                    3) 'maxlike' -- maximum likelihood scheme
            keep_parental:
                If True, keep the parental individuals in the resulting
                population.

        Returns:
            A new classified population with admixture coefficients.
        """

        if classifier == 'admixture':
            result = self._mix_classifier__admixture(pop, **kwargs)
        elif classifier == 'admixture':
            result = self._mix_classifier__structure(pop, **kwargs)
        elif classifier == 'freqs':
            result = self._mix_classifier__freqs(pop, classifier, **kwargs)
        elif classifier == 'maxlike':
            from scipy.optimize import minimize
            from kpop.admixture.likelihood import loglike_mixture

            def admixture(x, f_pops):
                g = (x == 1).sum(axis=-1)
                gbar = 2 - g
                q0 = np.zeros(len(f_pops), dtype=float)
                q0 += 1 / len(f_pops)

                def f(q):
                    return -loglike_mixture(x, q, f_pops)

                def fprime(q):
                    f = (q[:, None] * f_pops).sum(0)
                    aux = (g / f - gbar / (1 - f))[None, :] * f_pops
                    return -aux.sum(axis=-1)

                def cons(q):
                    return q.sum() - 1

                bounds = [(0, 1) for _ in q0]
                cons = {'type': 'eq', 'fun': cons}
                res = minimize(f, q0, bounds=bounds, constraints=[cons])

                if not res.success:
                    raise ValueError(res.message)
                return res.x

            kwargs['function'] = admixture
            result = self._mix_classifier__freqs(pop, **kwargs)
        elif classifier in ('maxlike', 'bayes'):
            result = self._mix_classifier__kpop(pop, classifier, **kwargs)
        elif classifier in ('sharp'):
            result = self._mix_classifier__posterior(pop, sharp=True, **kwargs)
        elif classifier in ('posterior'):
            result = self._mix_classifier__posterior(pop, **kwargs)
        else:
            raise ValueError('invalid classification method: %r' % classifier)

        result.parent = self
        if keep_parental:
            for i, pop in enumerate(self.populations):
                result.populations.insert(0, pop.copy())
        return result

    def _mix_classifier__freqs(self, pop, function, alpha=0.0):
        f_pops = np.array(self.freqs_vector_pop)
        f_pops += alpha
        f_pops /= 1 + 2 * alpha

        inds = [ind.copy() for ind in pop]
        for ind in inds:
            distrib = function(ind.data, f_pops)
            ind.admixture_q = Prob(enumerate(distrib))
        return pop.transformed_copy(inds, parent=self)

    def _mix_classifier__kpop(self, pop, method):
        from kpop.admixture import admixture

        f_pops = np.array(self.freqs_matrix_pop)
        inds = [ind.copy() for ind in pop]
        for ind in inds:
            distrib = admixture(ind.data, f_pops, method=method)
            ind.admixture_q = Prob(enumerate(distrib))
        return pop.transformed_copy(inds, parent=self)

    def _mix_classifier__posterior(self, pop, alpha=0.5, sharp=False):
        from kpop.admixture.likelihood import bayes

        f_pops = np.array(self.freqs_vector_pop)
        f_pops += alpha
        f_pops /= 1 + 2 * alpha
        inds = [ind.copy() for ind in pop]
        for ind in inds:
            values = bayes(ind.data, f_pops)
            ind.admixture_q = Prob(enumerate(values))
            if sharp:
                ind.admixture_q = ind.admixture_q.sharp()
        return pop.transformed_copy(inds, parent=self)

    def _mix_classifier__structure(self, pop, keep_parental):
        raise NotImplementedError

    def _mix_classifier__admixture(self, pop):
        # Create supervised pop_labels
        pop_labels = []
        for i, label in enumerate(self.populations.labels()):
            pop_labels.extend([label] * self.populations[i].size)
        pop_labels.extend([None] * pop.size)

        # Create resulting population
        pop_total = self + pop
        classified_pop = pop_total.structure_admixture(
            k=len(self.populations),
            parental_labels=self.populations.labels(),
            pop_labels=pop_labels,
        )
        del classified_pop.populations[0: self.num_populations]
        return classified_pop

    def fill_missing(self):
        for population in self.populations:
            population.fill_missing()

    def new_admixed_population(self, coeffs, size=0, label=None, **kwargs):
        """
        Creates a new empty admixed population with the given mixture
        coefficients. If size is given, fill it with the given number of
        individuals.

        Args:
            coeffs:
                Admixture coefficients. Can be a list of values or a Prob()
                object mapping sub-population labels or indexes to admixture
                coefficients.
            size:
                Optional size of the new population.
            label:
                Label for the new population.
            **kwargs:
                Additional coefficients passed to the Population constructor.

        Returns:
            A new Population.
        """

        # Convert Probs to a flat list of probabilities
        if isinstance(coeffs, Prob):
            coeffs, probs = [], dict(coeffs)
            for i, pop in enumerate(self.populations):
                if pop.label in coeffs:
                    coeffs.append(probs.pop(pop.label))
                elif i in probs:
                    coeffs.append(probs.pop(i))
                else:
                    coeffs.append(0.0)
            if probs:
                key = probs.popitem()[0]
                raise ValueError('%r is not a sub-population' % key)

        # Normalize
        coeffs = np.array(coeffs)
        coeffs /= coeffs.sum()

        # Initialize freqs using freqs_matrix data for a low number of alleles
        if self.num_alleles is not None and self.num_alleles <= 5:
            freqs = np.zeros((self.num_loci, self.num_alleles), dtype=float)
            for q, pop in zip(coeffs, self.populations):
                freqs += q * pop.freqs_matrix

        # Full freqs initialization
        else:
            freqs = []
            for i in range(self.num_loci):
                probs = [pop.freqs[i] for pop in self.populations]
                prob = Prob.mixture(coeffs, probs)
                freqs.append(prob)

        # Create population object
        kwargs['label'] = label
        kwargs.setdefault('parent', self)
        kwargs.setdefault('ploidy', self.ploidy)
        kwargs.setdefault('num_alleles', self.num_alleles)
        pop = kpop.Population(freqs=freqs, **kwargs)
        pop.fill(size)
        return pop

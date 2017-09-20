import numpy as np

from .attr import Attr
from ..admixture import admixture
from ..prob import Prob


class Admixture(Attr):
    def structure_admixture(self, k, pop_ids=None, parental_ids=None):
        """
        Runs the ADMIXTURE program to detect structure in the population.

        Args:
            k: number of parental populations.

        Returns:
            A new Population object with all individuals classified with their
            respective admixture coefficients.
        """

        from kpop.external.admixture import run_admixture
        pop = self._population

        if not pop.is_biallelic:
            raise ValueError('ADMIXTURE only supports biallelic populations')

        kwargs = {}
        if pop_ids:
            kwargs['supervised'] = pop_ids
        result = run_admixture(pop, k, disp=0, **kwargs)
        parental = result.make_parental(ids=parental_ids)
        individuals = result.make_admixture_ids(pop)
        out = pop.transformed_copy(individuals, parent=parental)
        out.admixture_result = result
        return out

    def admixture(self, ind, method='maxlike'):
        """
        Compute the admixture coefficients for individual

        Args:
            ind: a :class:`Individual` instance.

        Returns:
            A :class:`Prob` mapping from population ids to admixture
            coefficients.
        """

        if self.is_biallelic:
            f_pops = self.freqs_matrix_pop
            result = admixture.admixture(ind.data, f_pops, method=method)
            data = {}
            for i, pop in enumerate(self.populations):
                data[pop.id or i] = result[i]
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
        # Create supervised pop_ids
        pop_ids = []
        for i, id in enumerate(self.populations.population_ids()):
            pop_ids.extend([id] * self.populations[i].size)
        pop_ids.extend([None] * pop.size)

        # Create resulting population
        pop_total = self + pop
        classified_pop = pop_total.structure_admixture(
            k=len(self.populations),
            parental_ids=self.populations.population_ids(),
            pop_ids=pop_ids,
        )
        del classified_pop.populations[0: self.num_populations]
        return classified_pop

    def new_admixed_population(self, coeffs, size=0, id=None, **kwargs):
        """
        Creates a new empty admixed population with the given mixture
        coefficients. If size is given, fill it with the given number of
        individuals.

        Args:
            coeffs:
                Admixture coefficients. Can be a list of values or a Prob()
                object mapping sub-population ids or indexes to admixture
                coefficients.
            size:
                Optional size of the new population.
            id:
                Label for the new population.
            **kwargs:
                Additional coefficients passed to the Population constructor.

        Returns:
            A new Population.
        """
        import kpop

        # Convert Probs to a flat list of probabilities
        if isinstance(coeffs, Prob):
            coeffs, probs = [], dict(coeffs)
            for i, pop in enumerate(self.populations):
                if pop.id in coeffs:
                    coeffs.append(probs.pop(pop.id))
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
        kwargs['id'] = id
        kwargs.setdefault('parent', self)
        kwargs.setdefault('ploidy', self.ploidy)
        kwargs.setdefault('num_alleles', self.num_alleles)
        pop = kpop.Population(freqs=freqs, **kwargs)
        pop.fill(size)
        return pop

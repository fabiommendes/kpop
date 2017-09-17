from random import randrange

from .attr import Attr
from .individual import Individual
from ..simulation import frequency_drift


class Simulation(Attr):
    def genetic_drift(self, n_generations, population_size=None,
                      sample_size=None, id=None):
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
            id:
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

        pop = self._population
        population_size = population_size or pop.size
        eff_size = pop.ploidy * population_size
        id = id or _sub_population_id_labels(pop.id, n_generations)
        sample_size = pop.size if sample_size is None else sample_size

        # Compute drift
        freqs = frequency_drift(pop.freqs_matrix, n_generations, eff_size)

        # Create population
        new = pop._population(freqs=freqs, id=id)
        new.num_loci = pop.num_loci
        new.num_alleles = pop.num_alleles
        new.ploidy = pop.ploidy
        new.freqs_matrix = freqs
        new.freqs_vector = freqs[:, 0]

        # Fill with individuals
        if sample_size:
            new.fill(sample_size)
        return new

    def breed(self, population, size, id='breed', **kwargs):
        """
        Breed individuals from both populations.

        Args:
            population:
                Population to choose individuals from.
            size:
                Number of offspring.
            id:
                Label for the resulting population.
            **kwargs:

        Returns:
            A new population.
        """

        children = []
        for i in range(size):
            father, mother = self.random_individual(), self.random_individual()
            child = father.breed(mother, id='%s%s' % (id, i), **kwargs)
            children.append(child)
        return self._new_population(children, id=id)

    def new_individual(self, id=None, **kwargs):
        """
        Return a new random individual respecting the population allele
        frequencies.
        """

        pop = self._population
        kwargs['id'] = id or pop._next_id()
        ind = Individual.from_freqs(pop.freqs_matrix, population=pop, **kwargs)
        ind.num_alleles = pop.num_alleles
        return ind

    def new_offspring(self, i=None, j=None, **kwargs):
        """
        Return a new offspring created by breeding two random elements in the
        population.

        This individual is not added to the population. If you want that,
        just do:

        >>> pop.add(pop.new_offspring())                        # doctest: +SKIP
        """

        pop = self._population
        size = len(pop)
        if i is None:
            i = randrange(size)
        if j is None:
            j = randrange(size)
        return pop[i].breed(pop[j], population=pop, **kwargs)

    def random_individual(self):
        """
        Return a random individual from population.
        """

        i = randrange(len(self._population))
        return self._population[i]


def _sub_population_id_labels(label, n_gen):
    "Produces an id label for a sub-population."

    if label is None:
        return None
    else:
        return '%s-gen%s' % (label, n_gen)

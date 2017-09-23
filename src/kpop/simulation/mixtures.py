from ..libs import np

from ..population.utils import random_individual_data
from kpop.libs import lazy_module

kpop = lazy_module('kpop')


def linear_mixture(pop1, pop2, size, progression=lambda x: (x, 1 - x), id=None):
    """
    Creates a new population of given size in which each individual has a
    different pattern of mixture between the given reference populations.

    Args:
        pop1, pop2:
            Reference populations.
        size (int):
            Size of the new admixed population.
        progression:
            A function of a parameter between 0 and 1 that return a pair with
            admixture proportions between the two populations.
        id:
            Id for the resulting population.

    Returns:
        A new population object.
    """

    alphas = np.linspace(0, 1, size)
    proportions = [progression(alpha) for alpha in alphas]
    freq1 = getattr(pop1, 'freqs_matrix', pop1)
    freq2 = getattr(pop2, 'freqs_matrix', pop2)
    ploidy = pop1.ploidy
    data = []

    for p, q in proportions:
        z = p + q
        p, q = p / z, q / z
        freqs = p * freq1 + q * freq2
        individual_data = random_individual_data(freqs, ploidy=ploidy)
        data.append(individual_data)

    return kpop.Population(data, id=id)

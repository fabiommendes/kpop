import numpy as np


def fill_frequencies(freqs):
    """
    Return a Nx2 array with N tuples of (f, 1-f) for each element f of
    freqs_binomial.
    """

    if len(np.shape(freqs)) == 2:
        return np.asarray(freqs)

    out = np.ones((len(freqs), 2), dtype=float)
    out[:, 0] = freqs
    out[:, 1] -= freqs
    return out


def flatten_frequencies(freqs):
    """
    Flatten an Nx2 array of (f, 1 - f) frequency pairs to a Nx1 array of f's.
    """

    if len(np.shape(freqs)) == 2:
        return np.asarray(freqs)[:, 0]
    else:
        return np.asarray(freqs)


def freqs_binomial(population, alpha=0.0, flat=False):
    """
    Return the inferred frequencies for population with biallelic data.

    Each frequency is computed as (n + alpha)/(N + 2*alpha) in which n is the
    number of occurrences of the "1" allele, and N is the total number of
    occurrences. The parameter "alpha" appears in the prior Dirichlet
    distribution for the frequencies parameters and must be greater than 0.0.
    A good choice is the geometric prior which corresponds to alpha=0.5.

    Examples:
        Given some simple population, compute its allele frequencies

        >>> pop = [[[1, 2], [1, 1], [2, 2]],
        ...        [[2, 1], [1, 1], [1, 2]],
        ...        [[1, 2], [1, 1], [2, 2]],
        ...        [[2, 1], [1, 1], [1, 2]]]
        >>> freqs_binomial(pop, filled=False)
        array([ 0.5 ,  0.  ,  0.75])
    """

    pop = np.asarray(population, dtype=np.int8)
    size = pop.shape[0] * pop.shape[2]
    freqs = np.asarray((pop == 1).sum(axis=(0, 2)), dtype=float)
    if alpha:
        freqs += alpha
    freqs /= (size + 2 * alpha)
    return freqs if flat else fill_frequencies(freqs)

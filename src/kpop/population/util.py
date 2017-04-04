import collections

import numpy as np

from kpop.freqmatrix import fill_frequencies
from kpop.prob import Prob


def normalize_freqs_arg(freqs):
    """
    Normalize frequencies to be a sequence of probability distributions.
    """

    if freqs is None:
        return None

    if len(freqs) == 0:
        raise ValueError('cannot initialize from empty frequencies')

    if isinstance(freqs[0], collections.Mapping):
        result = []
        for prob in freqs:
            if not isinstance(prob, collections.Mapping):
                prob = enumerate(prob)
            prob = Prob(prob)
            result.append(prob)
        return result
    else:
        freqs = np.asarray(freqs)
        result = []
        if freqs.ndim == 1:
            for prob in freqs:
                result.append(Prob({1: prob, 2: 1 - prob}))
        else:
            for prob in freqs:
                data = enumerate(prob, start=1)
                result.append(Prob(data))
        return result


def random_frequencies(num_loci, alleles=2, clip=0.0, seed=None):
    """
    Return random frequency distributions over loci.
    """

    if alleles <= 1:
        raise ValueError('needs at least 2 different alleles')
    uniform = np.random.uniform
    if seed:
        np.random.seed(seed)
    if alleles == 2:
        return fill_frequencies(uniform(clip, 1 - clip, size=num_loci))
    else:
        data = uniform(0, 1, size=num_loci * alleles)
        data = data.reshape((num_loci, alleles))
        data /= data.sum(axis=1)[:, None]
        return data


def random_pop_data(size, freqs, ploidy=2, dtype=np.int8):
    """
    Create a random population data of bi-allele individuals from the given
    frequency.

    freqs_binomial can be a list of frequencies f[j] or a list of (f[j], 1-f[j]) tuples.
    In both cases, the frequency f[j] represents the probability that the
    j-th feature has a value of 1.
    """

    num_loci = len(freqs)
    data = np.random.uniform(size=num_loci * ploidy * int(size))
    data = data.reshape((size, num_loci, ploidy))
    alleles = freqs.shape[1]
    if alleles == 2:
        return np.asarray(data < freqs[:, 0, None], dtype=dtype) + 1
    else:
        data = np.zeros((size, num_loci, ploidy), dtype=dtype)
        mask = data >= freqs[:, 0, None]
        for i in range(1, alleles):
            data[mask] = i
            mask &= data >= freqs[:, i, None]
        return data + 1

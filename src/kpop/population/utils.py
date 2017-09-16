import collections

import numpy as np

from kpop.prob import Prob
from kpop.utils import freqs_to_matrix


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


def freqs_fastest(pop):
    try:
        return pop.freqs_vector
    except ValueError:
        pass
    if pop.num_alleles <= 5:
        return pop.freqs_matrix
    return pop.freqs


def get_freqs(pop):
    if pop._freqs is None:
        vars_ = pop.__dict__
        data = vars_.get('freqs_matrix') or vars_.get('freqs_vector')

        if data is not None:
            if data.ndim == 1:
                data = freqs_to_matrix(data)
            pop._freqs = [Prob(enumerate(locus)) for locus in data]

        else:
            pop._freqs = pop.stats.empirical_freqs()

    return pop._freqs


def set_freqs(pop, value):
    pop._freqs = normalize_freqs_arg(value)
    del pop.freqs_vector
    del pop.freqs_matrix


def hfreqs_vector(pop):
    if pop.ploidy == 2:
        data = np.array(pop)
        heterozygote = data[:, :, 0] != data[:, :, 1]
        return heterozygote.mean(axis=0)
    else:
        raise NotImplementedError('require diploid populations')

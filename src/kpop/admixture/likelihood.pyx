cimport numpy as np
cimport cython
import numpy
import numpy as np

from kpop.freqmatrix import flatten_frequencies, fill_frequencies


# Functions from C math lib
cdef extern from "math.h":
    double log(double)
    double sqrt(double)


def loglike(x, freqs):
    """
    Compute the log-like log(P(x|freqs)) for a given individual x and
    the frequencies freqs. The frequencies can be given both in the
    form of a list of f_j's or a list of pairs (f_j, 1 - f_j).

    The like model is

        Pr(x|freqs_binomial) = Prod_j fj^(uj + vj) * (1 - fj)^(uj' + vj')

    in which we defined xj = (uj, vj) and the primes as being the boolean
    negation of the variable (e.g., if u0 = 1, then u0' = 0). This function
    returns the logarithm of this model.

    Examples:
        Let us compute some values manually

        >>> x = [[1, 1], [2, 2], [1, 2]]
        >>> f1 = [0.5, 0.5, 0.5]
        >>> f2 = [0.25, 0.75, 0.5]

        The log like for the first population is just ``6 ln(0.5) =
        -4.1488...``, since it corresponds to 6 contributions with probability
        = 0.5.

        >>> loglike(x, f1)                                  # doctest: +ELLIPSIS
        -4.1588...

        The first two loci in the second population contributes with 4 terms of
        ln(0.75) and the last locus contributes with ln(0.5) twice. The gives an
        overall value of -2.5370...

        >>> loglike(x, f2)                                  # doctest: +ELLIPSIS
        -2.5370...

        The loglike_logf() function can speed up calculations a little bit if
        the logarithms of f's are precomputed. It accepts an argument in the
        form of a list of (log(f), log(1-f)) terms. The _results should be the
        same.

        >>> loglike_logf(x, np.log(fill_frequencies(f2)))   # doctest: +ELLIPSIS
        -2.5370...

        Another convenience is loglike_mix(), which automatically mix the two
        populations with some arbitrary proportions.

        >>> fmix = [0.375, 0.625, 0.5]  # 0.5 * f1 + 0.5 * f2
        >>> loglike(x, fmix) == loglike_mix(x, [0.5, 0.5], [f1, f2])
        True

    See also:
        loglike_logf() --- useful when the log's of freqs_binomial are precomputed
        loglike_mix()  --- like of a mixture between 2 or more populations
        like()         --- like function, without logarithms
    """

    x = np.asarray(x, dtype=np.uint8)
    freqs = np.ascontiguousarray(flatten_frequencies(freqs))
    return _loglike(x, freqs)


def loglike_logf(x, logfreqs):
    """
    Compute the log-like log(P(x|freqs)) for a given individual x and
    the logarithm of the frequencies logfreqs.

    Args:
        x:
            Sequence of loci
        logfreqs:
            Sequence of [(log(f0), log(1-f0)), ..., (log(fJ), log(1 - fJ))]

    See also:
        loglike() --- Similar, but expect non-logarithm frequencies. The
            docstring contains examples and more information on the like
            model.
    """

    x = np.asarray(x, dtype=np.uint8)
    freqs = np.ascontiguousarray(logfreqs)
    return _loglike_logf(x, freqs)


def like(x, freqs):
    """
    Compute the like P(x|freqs), instead of its logarithm.

    Using log-likelihoods is more stable. This method exists since direct
    likelihood values are easier to interpret.

    See also:
        loglike() --- Examples and more information on the like model.
    """

    x = np.asarray(x, dtype=np.uint8)
    freqs = np.ascontiguousarray(flatten_frequencies(freqs))
    return _like(x, freqs)


def loglike_mix(x, admix_q, f_pops):
    """
    Compute the log-like log(P(x|freqs_binomial)) for an individual x and
    freqs_binomial given by a particular mix of populations in f_pops with
    proportions admix_q.

    See also:
        loglike() --- Examples and more information on the like model.
    """

    x = np.asarray(x, dtype=np.uint8)
    admix_q = np.asarray(admix_q)
    f_pops = np.asarray(f_pops)
    return _loglike_mix(x, admix_q, f_pops)


def loglike_mixture(x, mixture_q, f_pops):
    """
    Compute the log-like log(P(x|freqs_binomial)) for an individual x and
    freqs_binomial given by a particular mix of populations in f_pops with
    proportions admix_q.

    See also:
        loglike() --- Examples and more information on the like model.
    """

    # freqs = (mixture_q[:, None] * f_pops).sum(axis=0)
    # return loglike(x, freqs)

    x = np.asarray(x, dtype=np.uint8)
    mixture_q = np.asarray(mixture_q)
    f_pops = np.asarray(f_pops)
    return _loglike_mix(x, mixture_q, f_pops)


def logbayes(x, f_pops, *, prior=None):
    """
    Return a list of the logarithm of Bayes factors for individual x belonging
    to each population defined by the frequencies in f_pops.

    An optional prior probability can be assigned to each population.


    Examples:
        Consider this individual and some simple ancestral populations

        >>> x = np.array([[2, 2], [2, 1]], dtype=np.uint8)
        >>> f1 = [0.75, 0.75]
        >>> f2 = [0.25, 0.25]

        We call logbayes() with frequencies for each population

        >>> B1, B2 = logbayes(x, [f1, f2])

        Let us compare with manually computed values

        >>> B1 == 3 * np.log(0.75) + np.log(0.25)
        True
        >>> B2 == 3 * np.log(0.25) + np.log(0.75)
        True

        It is more likely that x was drawn from f1 that from f2. This can be
        seen by bayes_classify(),

        >>> bayes_classify(x, [f1, f2])
        0

        Sometimes is more instructive to see the normalized Bayes factors in
        order to have a better idea of how much one option is favored over the
        other

        >>> bayes(x, [f1, f2])
        array([ 0.9,  0.1])

    See also:
        bayes_classify() --- select the index of the most probable population
            for the given individual.
        bayes()  --- return a list of normalized Bayes factors, rather
            than logarithms.
        ilogbayes_logf() --- the log-frequencies version of this function.
    """

    num_pops = len(f_pops)
    if prior is None:
        log_probs = np.zeros(num_pops, dtype=float)
        log_probs -= np.log(num_pops)
    else:
        log_probs = np.log(prior)

    x = np.asarray(x, dtype=np.uint8)
    log_probs += [loglike(x, f) for f in f_pops]
    return log_probs


def bayes(ind, f_pops, *, prior=None):
    """
    Compute the normalized Bayes factors from a list of population frequencies
    and optionally prior probabilities.

    Bayes factors are calculated from

        P(freq_i|ind) = P(freq_i) P(ind|freq_i) / P(ind)

    The term P(ind) is just a normalization constant whereas P(freq_i) is the
    prior probability defined by "prior".
    """

    log_probs = logbayes(ind, f_pops, prior=prior)
    #max_log_prob = log_probs.max()
    #log_probs -= max_log_prob
    probs = np.exp(log_probs)
    probs /= probs.sum()
    return probs


def bayes_classify(ind, f_pops, *, prior=None):
    """
    Return the index of the population which has the greatest Bayes factor for
    the individual "ind".
    """

    return logbayes(ind, f_pops, prior=prior).argmax()


#
# C implementations
#
cdef double _loglike_logf(unsigned char[:, :] x, double[:, :] logfreqs):
    """
    Fast implementation of core.loglike_logf().
    """

    cdef int j, a
    cdef double loglike = 0.0
    cdef int num_loci = x.shape[0]
    cdef int ploidy = x.shape[1]

    for j in range(num_loci):
        for a in range(ploidy):
            if x[j, a] == 1:
                loglike += logfreqs[j, 0]
            elif x[j, a] == 2:
                loglike += logfreqs[j, 1]
            else:
                raise ValueError('invalid allele value: %s' % x[j, a])
    return loglike


cdef double _loglike(unsigned char[:, :] x, double[:] freqs):
    """
    Fast implementation of core.loglike().
    """

    cdef int j, a
    cdef double loglike = 0.0
    cdef int num_loci = x.shape[0]
    cdef int ploidy = x.shape[1]

    for j in range(num_loci):
        for a in range(ploidy):
            if x[j, a] == 1:
                loglike += log(freqs[j])
            elif x[j, a] == 2:
                loglike += log(1 - freqs[j])
            else:
                raise ValueError('invalid allele value: %s' % x[j, a])
    return loglike


cdef double _like(unsigned char[:, :] x, double[:] freqs):
    """
    Fast implementation of core.like()
    """

    cdef int j, a
    cdef double like = 1.0
    cdef int num_loci = x.shape[0]
    cdef int ploidy = x.shape[1]

    for j in range(num_loci):
        for a in range(ploidy):
            if x[j, a] == 1:
                like *= freqs[j]
            elif x[j, a] == 2:
                like *= 1 - freqs[j]
            else:
                raise ValueError('invalid allele value: %s' % x[j, a])
    return like


cdef double _loglike_mix(unsigned char[:, :] x, double[:] mixture, double[:, :] f_pops):
    """
    Fast implementation of core.loglike_mix().
    """

    cdef int j, k, a
    cdef double fj, loglike = 0.0
    cdef int num_loci = x.shape[0]
    cdef int k_max = mixture.shape[0]
    cdef int ploidy = x.shape[1]

    for j in range(num_loci):
        fj = 0
        for k in range(k_max):
            fj += mixture[k] * f_pops[k, j]

        for a in range(ploidy):
            if x[j, a] == 1:
                loglike += log(fj)
            elif x[j, a] == 2:
                loglike += log(1 - fj)
            else:
                raise ValueError('invalid allele value: %s' % x[j, a])
    return loglike

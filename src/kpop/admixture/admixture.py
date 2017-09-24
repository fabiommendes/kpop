from kpop.admixture.admixture_maxlike import admixture_maxlike
from .likelihood import loglike, logbayes
from ..libs import np, sp_integrate


def admixture_em_single(x, f_pops):
    f_pops = np.asarray(f_pops)
    k = len(f_pops)
    g = (np.asarray(x) == 1).sum(axis=1)
    q = np.ones(k, dtype=float) / k
    q_ = q[:, None]
    a = q_ * f_pops / (q_ * f_pops).sum(axis=0)
    b = q_ * (1 - f_pops) / (q_ * (1 - f_pops)).sum(axis=0)
    num_loci = len(x)
    return (g * a + (2 - g) * b).sum(axis=1) / (2 * num_loci)


def admixture(x, f_pops, method='bayes'):
    """
    Return the reference admixture proportions for the given individual.

    These are the allowed methods:
        'maxlike':
            Return the proportion that maximizes the individual likelihood.
            (This is the fastest method and usually yields similar values to
            Bayesian methods)
        'bayes':
            Compute the proportion as the mean value over the posterior
        'bayesfull':
            Similar to 'bayes', but return a list of pairs
            (proportion[i], delta[i]) with delta computed as the standard
            deviation of the i-th proportion.

    Returns:
        An array with the admixture coefficients.
    """
    if method == 'maxlike':
        return admixture_maxlike(x, f_pops)

    elif method in ('bayes', 'fullbayes'):
        full = method == 'fullbayes'
        alpha = 0.5
        f1, f2 = f_pops.T
        samples = np.array([loglike(x, p * f1 + (1 - p) * f2) for p in
                            np.linspace(0, 1, 5)])
        log_max = samples.max()

        def pdf(p):
            ln_like = loglike(x, p * f1 + (1 - p) * f2)
            return np.exp(ln_like - log_max) * (p * (1.0 - p)) ** (alpha - 1.0)

        def p_func(p):
            return p * pdf(p)

        def p2_func(p):
            return (p - pmean) * (p - pmean) * pdf(p)

        # Compute the normalization factor
        normalization, _ = sp_integrate.quad(pdf, 0, 1)
        Ip, _ = sp_integrate.quad(p_func, 0, 1)
        p = pmean = Ip / normalization

        if full:
            Ip2, _ = sp_integrate.quad(p2_func, 0, 1)
            dp = np.sqrt(Ip2 / normalization)
            return np.array([[p, dp], [1 - p, dp]])
        else:
            return np.array([p, 1 - p])
    else:
        raise ValueError('invalid method, %s' % method)


def admixture_logbayes(x, f_pops, normalize=False):
    """
    Return a dictionary with the log Bayes factors for each possible combination
    of elements in f_pops.

    ::

        D[i]    --> log(P("belongs to population i"| x))
        D[i, j] --> log(P("admixed between i and j"| x))


    See also:
        admixture_bayes: returns the normalized Bayes factors.
    """

    map = {}
    for i, logBF in enumerate(logbayes(x, f_pops)):
        map[i] = logBF

    for i, fi in enumerate(f_pops):
        for j, fj in enumerate(f_pops, i + 1):
            samples = np.array([loglike(x, p * fi + (1 - p) * fj) for p in
                                np.arange(0, 1.1, 0.25)])
            log_min = samples.min()
            if np.isinf(log_min):
                log_min = min(x for x in samples if not np.isinf(x))

            def func(p):
                return np.exp(
                    loglike(x, p * fi + (1 - p) * fj) - log_min) / np.sqrt(
                    p * (1 - p))

            Z, _ = sp_integrate.quad(func, 0, 1)
            map[i, j] = np.log(Z) + log_min

    if normalize:
        min_v = min(x for x in map.values() if not np.isinf(x))
        map = {k: v - min_v for (k, v) in map.items()}
    return map


def admixture_bayes(x, f_pops):
    """
    Similar to admixture_logbayes, but return a dictionary with the normalized
    Bayes factors for all combinations of mixtures. Hence::

        D[i]    --> P("belongs to population i"| x)
        D[i, j] --> P("admixed between i and j"| x)

    """
    map = admixture_logbayes(x, f_pops, normalize=True)
    map = {k: np.exp(v) for (k, v) in map.items()}
    norm = sum(map.values())
    map = {k: v / norm for (k, v) in map.items()}
    return map

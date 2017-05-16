import math

import numpy as np
import scipy.integrate
import scipy.optimize

from kpop.admixture.likelihood import loglike, loglike_mix, loglike_logf, \
    logbayes


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

    elif method == 'bayes' or method == 'fullbayes':
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
        normalization, _ = scipy.integrate.quad(pdf, 0, 1)
        Ip, _ = scipy.integrate.quad(p_func, 0, 1)
        p = pmean = Ip / normalization

        if full:
            Ip2, _ = scipy.integrate.quad(p2_func, 0, 1)
            dp = np.sqrt(Ip2 / normalization)
            return np.array([[p, dp], [1 - p, dp]])
        else:
            return np.array([p, 1 - p])
    else:
        raise ValueError('invalid method, %s' % method)


def admixture_maxlike(x, f_pops, optimizer='slsqrp', q0=None):
    """
    Uses maximum likelihood in order to estimate the admixture parameters

    Usage:
        Let us create some random populations and a mixed individual

        >>> from kpop import Population, Individual
        >>> p1 = Population.make_random(10, 1000)
        >>> p2 = Population.make_random(10, 1000)
        >>> x = Individual.from_freqs(0.5 * p1.freqs + 0.5 * p2.freqs)

        We correctly infer that the admixture coefficients are close to the correct
        value of 0.5

        >>> p, q = admixture_maxlike(x, [f1, f2])
        >>> assert abs(p - 0.5) < 0.1

        This also work with more than 2 ancestral populations

        >>> p3 = Population.make_random(10, 1000)
        >>> q1, q2, q3 = admixture_maxlike(x, [f1, f2, f3])
    """

    if len(f_pops) == 2 and optimizer in ['bound', 'brent', 'cg', 'ncg',
                                          'bfgs']:
        return admixture_maxlike_1d(x, f_pops, method=optimizer)

    # Cache some values
    f_pops = np.asarray(f_pops)
    f_pops_single = np.ascontiguousarray(f_pops[:, :, 0])
    xbar = np.asarray(x == 1, dtype=np.int8)
    g = xbar.sum(axis=1)
    gbar = x.shape[-1] - g
    k = len(f_pops)

    # Define the objective function, its gradient and hessian
    # TODO: probably better to implement in Cython!
    def func(q):
        return -loglike_mix(x, q, f_pops_single)

    def fprime(q):
        f = (q[:, None] * f_pops_single).sum(0)
        aux = (g / f - gbar / (1 - f))[None, :] * f_pops_single
        return -aux.sum(axis=-1)

    def fhess(q):
        hess = np.empty((k, k), dtype=float)
        fp = f_pops_single
        f = (q[:, None] * fp).sum(0)
        aux = (g / f ** 2 - gbar / (1 - f) ** 2)
        for i in range(k):
            for j in range(i, k):
                ff = fp[i, :] * fp[j, :]
                hess[i, j] = (aux * ff).sum()
                hess[j, i] = hess[i, j]
        return hess

    # Define constraints, gradient and hessian
    neg_ones_vector = -np.ones(k, dtype=float)
    eqcons_result = neg_ones_vector[None, :].T
    null_elements_matrix = np.zeros((k, k), dtype=float)
    identity = np.zeros((k, k), dtype=float)
    for i in range(k):
        identity[i, i] = 1

    def eqcons(q):
        return 1 - q.sum()

    def eqcons_fprime(q):
        return eqcons_result

    def eqcons_hess(q):
        return null_elements_matrix

    # Define initial condition and bounds
    q0 = np.ones(k) / k if q0 is None else np.array(q0)
    bounds = np.array([(0, 1) for _ in range(k)])

    # Minimize using sequential least squares programming
    fmin = scipy.optimize.fmin_slsqp
    if optimizer == 'slsqrp':
        return fmin(
            # Function and gradient definition
            func, q0,
            fprime=fprime,

            # Equality constraints
            eqcons=[eqcons],
            fprime_eqcons=eqcons_fprime,

            # Equality constraints
            bounds=bounds,  # Bounds on 0<=qi<=1
            disp=1,
        )

    # Like before, but does not use gradients
    elif optimizer == 'slsqrp_f':
        return fmin(func, q0, eqcons=[eqcons], bounds=bounds, disp=1)
    else:
        raise ValueError('unknown method: %s' % optimizer)


def admixture_maxlike_1d(x, f_pops, method='bound'):
    """
    Find the maximum likelihood using one of many possible methods: 'bound',
    'brent', 'cg', 'ncg', 'bfgs'. Usually 'bound' is the fastest method.
    """

    # Define populations
    f_pops = np.asarray(f_pops)
    f1, f2 = f_pops
    like = loglike_logf
    arctan = math.atan
    pi = math.pi

    # Prepare data
    xbar = np.asarray(x == 0, dtype=np.int8)
    g = x.sum(1)
    gbar = xbar.sum(1)

    # Define minimization function using q as argument
    def func_q(q):
        return -like(x, np.log(q * f1 + (1 - q) * f2))

    # Define minimization function and derivative using unconstrained x=1/(1-q)
    def func_u(u):
        q = 0.5 + arctan(u) / pi
        return -like(x, np.log(q * f1 + (1 - q) * f2))

    def fprime_u(u):
        q = 0.5 + arctan(u) / pi
        f = q * f1[:, 0] + (1 - q) * f2[:, 0]
        aux = (g / f - gbar / (1 - f)) * (f1[:, 0] - f2[:, 0])
        fprime = -aux.sum() / (1 + u ** 2) / pi
        return fprime

    # Optimize function using any of the desired methods
    if method == 'bound':
        q = scipy.optimize.fminbound(func_q, 0, 1)

    elif method == 'brent':
        u = scipy.optimize.brent(func_u)
        q = 0.5 + arctan(u) / pi

    elif method == 'bfgs':
        u = scipy.optimize.fmin_bfgs(func_u, np.array(0), fprime_u, disp=0)
        q = 0.5 + arctan(u) / pi

    elif method == 'cg':
        u = scipy.optimize.fmin_cg(func_u, np.array(0), fprime_u, disp=0)
        q = 0.5 + arctan(u) / pi

    elif method == 'ncg':
        u = scipy.optimize.fmin_ncg(func_u, np.array(0), fprime_u, disp=0)
        q = 0.5 + arctan(u) / pi

    else:
        raise ValueError('invalid method: %s' % method)

    return np.array([q, 1 - q])


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

    for i in range(len(f_pops)):
        for j in range(i + 1, len(f_pops)):
            f1, f2 = f_pops[i], f_pops[j]
            samples = np.array([loglike(x, p * f1 + (1 - p) * f2) for p in
                                np.arange(0, 1.1, 0.25)])
            log_min = samples.min()
            if np.isinf(log_min):
                log_min = min(x for x in samples if not np.isinf(x))

            def func(p):
                return np.exp(
                    loglike(x, p * f1 + (1 - p) * f2) - log_min) / np.sqrt(
                    p * (1 - p))

            Z, _ = scipy.integrate.quad(func, 0, 1)
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

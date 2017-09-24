import math

from .likelihood import loglike_mix, loglike_logf
from ..libs import np, sp_optimize


class AdmixtureMaxlike:
    """
    A namespace that implements the admixture_maxlike function.
    """

    _1d_methods = ['bound', 'brent', 'cg', 'ncg', 'bfgs']

    def __init__(self, x, f_pops, optimizer, q0):
        self.x = x
        self.f_pops = np.asarray(f_pops)
        self.f_pops_single = np.ascontiguousarray(self.f_pops[:, :, 0])
        self.optimizer = optimizer
        self.q0 = q0

        # Cache some values
        self.xbar = np.asarray(x == 1, dtype=np.int8)
        self.g = self.xbar.sum(axis=1)
        self.gbar = x.shape[-1] - self.g
        self.k = k = len(f_pops)
        self.eqcons_result = (-np.ones(k, dtype=float))[None, :].T
        self.null_elements_matrix = np.zeros((k, k), dtype=float)

    def run(self):
        x, f_pops, optimizer, q0 = self.x, self.f_pops, self.optimizer, self.q0

        if len(f_pops) == 2 and optimizer in self._1d_methods:
            return admixture_maxlike_1d(x, f_pops, method=optimizer)

        # Define constraints, gradient and hessian
        k = self.k
        identity = np.zeros((k, k), dtype=float)
        for i in range(k):
            identity[i, i] = 1

        # Define initial condition and bounds
        q0 = np.ones(k) / k if q0 is None else np.array(q0)
        bounds = np.array([(0, 1) for _ in range(k)])

        # Minimize using sequential least squares programming
        fmin = sp_optimize.fmin_slsqp
        if optimizer == 'slsqrp':
            return fmin(
                # Function and gradient definition
                self.func, self.q0,
                fprime=self.fprime,

                # Equality constraints
                eqcons=[self.eqcons],
                fprime_eqcons=self.eqcons_fprime,

                # Equality constraints
                bounds=bounds,  # Bounds on 0<=qi<=1
                disp=1,
            )

        # Like before, but does not use gradients
        elif optimizer == 'slsqrp_f':
            return fmin(self.func, q0, eqcons=[self.eqcons], bounds=bounds,
                        disp=1)
        else:
            raise ValueError('unknown method: %s' % optimizer)

    #
    # Objective function
    #
    def func(self, q):
        return -loglike_mix(self.x, q, self.f_pops_single)

    def fprime(self, q):
        f = (q[:, None] * self.f_pops_single).sum(0)
        aux = (self.g / f - self.gbar / (1 - f))[None, :] * self.f_pops_single
        return -np.sum(aux, axis=-1)

    def fhess(self, q):
        hess = np.empty((self.k, self.k), dtype=float)
        fp = self.f_pops_single
        f = (q[:, None] * fp).sum(0)
        aux = (self.g / f ** 2 - self.gbar / (1 - f) ** 2)
        for i in range(self.k):
            for j in range(i, self.k):
                ff = fp[i, :] * fp[j, :]
                hess[i, j] = np.sum(aux * ff)
                hess[j, i] = hess[i, j]
        return hess

    #
    # Constraints
    #
    def eqcons(self, q):
        return 1 - q.sum()

    def eqcons_fprime(self, q):
        return self.eqcons_result

    def eqcons_hess(self, q):
        return self.null_elements_matrix


def admixture_maxlike(x, f_pops, optimizer='slsqrp', q0=None):
    """
    Uses maximum likelihood in order to estimate the admixture parameters

    Usage:
        Let us create some random populations and a mixed individual

        >>> from kpop import Population, Individual
        >>> p1 = Population.random(10, 1000)
        >>> p2 = Population.random(10, 1000)
        >>> x = Individual.from_freqs(0.5 * p1.freqs + 0.5 * p2.freqs)

        We correctly infer that the admixture coefficients are close to the correct
        value of 0.5

        >>> p, q = admixture_maxlike(x, [f1, f2])
        >>> assert abs(p - 0.5) < 0.1

        This also work with more than 2 ancestral populations

        >>> p3 = Population.random(10, 1000)
        >>> q1, q2, q3 = admixture_maxlike(x, [f1, f2, f3])
    """


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
        q = sp_optimize.fminbound(func_q, 0, 1)

    elif method == 'brent':
        u = sp_optimize.brent(func_u)
        q = 0.5 + arctan(u) / pi

    elif method == 'bfgs':
        u = sp_optimize.fmin_bfgs(func_u, np.array(0), fprime_u, disp=0)
        q = 0.5 + arctan(u) / pi

    elif method == 'cg':
        u = sp_optimize.fmin_cg(func_u, np.array(0), fprime_u, disp=0)
        q = 0.5 + arctan(u) / pi

    elif method == 'ncg':
        u = sp_optimize.fmin_ncg(func_u, np.array(0), fprime_u, disp=0)
        q = 0.5 + arctan(u) / pi

    else:
        raise ValueError('invalid method: %s' % method)

    return np.array([q, 1 - q])

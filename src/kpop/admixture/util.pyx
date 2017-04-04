cimport cython
import itertools
import numpy as np


cpdef force_bounds(double[:, :] F, double x0=0, double x1=1):
    """
    Take a matrix F and transform all values to be inside the given range.

    Changes are made *inplace*.

    Args:
        F (n x m ndarray):
            A frequency matrix
        x0, x1:
            Interval of valid values.
    """

    cdef double f
    cdef int i

    for i in range(F.shape[0]):
        for j in range(F.shape[1]):
            f = F[i, j]
            if f > x1:
                F[i, j] = x1
            elif f < x0:
                F[i, j] = x0


cpdef force_probs(double[:, :] Q, double x0=0, double x1=1):
    """
    Take a matrix of probability distributions and for all values to be bounded
    in the given range and force normalization to 1.

    Changes are made *inplace*.

    Args:
        Q (n x m ndarray):
            A frequency matrix. All rows should sum 1.
        x0, x1:
            Interval of valid values.
    """

    cdef double q, norm
    cdef int i, k

    for i in range(Q.shape[0]):
        norm = 0.0

        for k in range(Q.shape[1]):
            q  = Q[i, k]
            if q > x1:
                Q[i, k] = x1
                norm += x1
            elif q < x0:
                Q[i, k] = x0
                norm += x0
            else:
                norm += q

        if norm == 0.0:
            for k in range(Q.shape[1]):
                Q[i, k] = 1.0 / k
        elif norm != 1.0:
            for k in range(Q.shape[1]):
                Q[i, k] /= norm


def k_reordering(self, *args):
    """
    Reorder populations by the new given indexes.
    """

    if sorted(args) != list(range(self.k)) or len(set(args)) != self.k:
        raise ValueError('incomplete/invalid permutation: ' + str(args))

    q = np.array([self.q[:, k] for k in args]).T
    self.q = np.ascontiguousarray(q)

    f = np.array([self.freqs[:, k] for k in args]).T
    self.freqs = np.ascontiguousarray(f)


def reorder_k_to_best_q_fit(self, obj_q, disp=0):
    """
    Reorder populations such as the best fit with the objective admixture
    proportions is found.
    """

    best = None
    best_dist = float('inf')
    q = self.q
    freqs = self.freqs
    for perm in itertools.permutations(range(self.k), self.k):
        self.q = q
        self.freqs = freqs
        self.k_reordering(*perm)
        distance = abs(obj_q - self.q).sum()
        if distance < best_dist:
            best = self.q, self.freqs
            best_dist = distance
    self.q, self.freqs = best
    if disp:
        print('best q-fit distance:', best_dist)


def plot(self, show=True, **kwargs):
    """
    Plot all found admixture coefficients.
    """

    from kpop.plots import admixture_bars
    from matplotlib import pyplot as plt

    admixture_bars(self.q, **kwargs)
    if show:
        plt.show()

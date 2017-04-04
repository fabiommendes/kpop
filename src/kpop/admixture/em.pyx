#cython: boundscheck=False, cdivision=True, wraparound=False, initializedcheck=False, overflowcheck=False

cimport cython
import numpy as np
from .util cimport force_bounds, force_probs


cdef extern from "alloca.h":
    void* alloca(size_t) nogil

cdef extern from "math.h":
    double sqrt(double) nogil


#
# EM ALGORITHM
#
# These functions implements the EM iterations for solution of the
# max-likelihood problem. EM is a slow algorithm and probably should be used
# with some kind of acceleration.
#
# Speed comparison
#
# The EM step is about 2-3 times *faster* then a single objective function
# evaluation. This is mainly due to the log() terms in phi vs a loop in K in EM.
# For a large K, em(Q, F, G) might be slower than phi(Q, F, G).
#
# In either case, we don't want to evaluate phi() at every step of an EM update
# since that incurs in a noticeable overhead. For this reason, the accelerated
# algorithm accepts the new solution unconditionally, even if it would increase
# the value of the objective function.
#
def em(Q, F, G, accelerate=True, inplace=False, epsilon=1e-50):
    """
    Implements a single non-accelerated EM step and updates the Q and F
    matrices.

    Args:
        Q (array):
            (num_ind, k) array of admixture fractions
        F (array):
            (num_loci, k) array of parental frequencies
        G (array):
            (num_ind, num_loci) array of first allele counts
        accelerate (bool):
            if True, uses Cython implementation. There is no reason to disable
            acceleration unless you are testing kpop.
        inplace (bool):
            if True, reuses the Q, F matrices in the output
        epsilon (float):
            stabilization parameter added to logarithm arguments (should be
            close to zero)

    Returns:
        The updated Q and F matrices
    """
    if accelerate:
        if not inplace:
            return _em(Q, F, G, epsilon)
        else:
            Q[:], F[:] = _em(Q, F, G, epsilon)
            return Q, F

    num_loci = F.shape[0]

    # Update A and B matrices
    A = Q[:, None, :] * F[None, ...]
    A /= A.sum(2)[..., None]

    B = Q[:, None, :] * (1 - F)[None, ...]
    B /= B.sum(2)[..., None]

    # Update frequencies and q's
    GA = A * G[:, :, None]
    GB = B * (2 - G)[:, :, None]
    sum_ga_i = GA.sum(0)
    sum_gb_i = GB.sum(0)
    Fout = sum_ga_i / (sum_ga_i + sum_gb_i)
    Qout = (GA.sum(1) + GB.sum(1)) / (2 * num_loci)
    if inplace:
        Q[:], F[:] = Qout, Fout
    return Qout, Fout


def em_squarem(Q, F, G, inplace=False, stepper=em, args=(), kwargs=None, max_step=5):
    """
    SQUAREM acceleration for the EM iteration. This performs two EM evaluations
    per step, but usually converges much faster than a double evaluation.

    Args:
        Q, F, G, inplace:
            same as in :func:`em`
        stepper:
            function that computes a single iteration (defaults to :func:`em`)
        args, kwargs:
            arguments to the stepper function
        max_step:
            maximum size of an acceleration step (None for unbounded). A step
            size of 1 means a non-accelerated algorithm.
    """
    Q0, F0 = Q, F
    kwargs = kwargs or {}
    Q1, F1 = stepper(Q0, F0, G, *args, **kwargs)
    Q2, F2 = stepper(Q1, F1, G, *args, **kwargs)

    # r and v parameters
    rQ = Q1 - Q0
    rF = F1 - F0
    vQ = (Q2 - Q1) - rQ
    vF = (F2 - F1) - rF

    # Step size
    vv = np.sum(vQ * vQ) + np.sum(vF * vF)
    rr = np.sum(rQ * rQ) + np.sum(rF * rF)
    step_sqs3 = sqrt(rr / vv)
    step = max(1.0, step_sqs3)

    # Alternative step sizes (not implemented)
    # rv = np.sum(rQ * vQ) + np.sum(rF * vF)
    # step_sqs1 = -rv / vv
    # step_sqs2 = -rr / rv
    # steps = [step_sqs1, step_sqs2, step_sqs3]
    # step = max(min(s for s in steps if s > 0), 1.0)

    if np.isfinite(step):
        # If step == 1, there is no acceleration. We do not allow step
        # size to be very large since this will drive solution away
        # from the normalization and positivity constraints
        if max_step is None:
            max_step = float('inf')
        step = min(step, max_step)

        # Acceleration step
        Qacc = Q0 + 2 * rQ * step + vQ * step**2
        Facc = F0 + 2 * rF * step + vF * step**2

        # Project back into allowed manifold
        force_probs(Qacc, 0, 1)
        force_bounds(F, 0, 1)
    else:
        Qacc, Facc = Q2, F2

    if inplace:
        Q[:], F[:] = Qacc, Facc

    return Qacc, Facc


cdef _em(double[:, :] Q, double[:, :] F, int[:, :] G, double e=1e-50):
    # About 7-8x faster than the pure-Python version

    # Shapes
    cdef int ii = Q.shape[0], jj = F.shape[0], kk = Q.shape[1]

    # Matrices
    Q_out = np.zeros((ii, kk), dtype=float)
    F_out = np.zeros((jj, kk), dtype=float)
    cdef double[:, :] _Q = Q_out, _F = F_out

    # Stack allocate temporary vectors for freqs. This should be safe since
    # k is small.
    cdef double* f_j = <double*> alloca(sizeof(double) * kk)
    cdef double* F_j = <double*> alloca(sizeof(double) * kk)

    # Variables
    cdef int i, j, k
    cdef int n, n_comp
    cdef double g_ij, mean_f_ij, mean_F_ij, f_jk, q_ik, ga_ijk, gb_ijk
    cdef double Q_norm = 1.0 / (2.0 * jj)

    with nogil:
        # We iterate in j before i since usually jj > ii and this
        # gives a more efficient access to the _F[j, k] matrices
        for j in range(jj):
            for k in range(kk):
                f_j[k] = F_j[k] = 0.0

            for i in range(ii):
                mean_f_ij = mean_F_ij = e
                n = G[i, j]
                n_comp = 2 - n

                for k in range(kk):
                    f_jk = F[j, k]
                    q_ik = Q[i, k]
                    mean_f_ij += q_ik * f_jk
                    mean_F_ij += q_ik * (1.0 - f_jk)

                for k in range(kk):
                    q_ik = Q[i, k]
                    f_jk = F[j, k]
                    ga_ijk = n * q_ik * f_jk / mean_f_ij
                    gb_ijk = n_comp * q_ik * (1.0 - f_jk) / mean_F_ij
                    _Q[i, k] += ga_ijk + gb_ijk
                    f_j[k] += ga_ijk
                    F_j[k] += gb_ijk

            # Save results in F
            for k in range(kk):
                _F[j, k] = f_j[k] / (f_j[k] + F_j[k] + e)

        # Normalize Q
        for i in range(ii):
            for k in range(kk):
                _Q[i, k] *= Q_norm

    return Q_out, F_out

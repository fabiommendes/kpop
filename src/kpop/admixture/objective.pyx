# cython: boundscheck=False, cdivision=True, wraparound=False, initializedcheck=False, overflowcheck=False, language_level=3

"""
Definitions for the objective function phi(Q, F; G) and its decompositions and
their respective gradients, hessians, etc.
"""

cimport cython
cimport numpy as np
import numpy as np


cdef extern from "math.h" nogil:
    double log(double)
    double sqrt(double)


cdef extern from "alloca.h" nogil:
    void* alloca(size_t)


cdef class State:
    cdef int ii, jj, kk
    cdef double[:, :] Q
    cdef double[:, :] F
    cdef int[:, :] G

    def __cinit__(self, double[:, :] Q, double[:, :] F, int[:, :] G):
        """
        Return state from input data.
        """

        self.Q = Q
        self.F = F
        self.G = G

        if G is None:
            raise ValueError('G matrix is required')
        self.ii = G.shape[0]
        self.jj = G.shape[1]

        if Q is None and F is None:
            raise ValueError('Q or F matrix is required')
        if (Q is not None and Q.shape[0] != self.ii or
                        F is not None and F.shape[0] != self.jj or
                (Q is not None and F is not None and Q.shape[1] != F.shape[1])):
            raise ValueError('Invalid shapes: Q[%s, %s], F[%s, %s], G[%s, %s]' % (
                Q.shape[0], Q.shape[1], F.shape[0], F.shape[1], self.ii, self.jj
            ))
        if Q is not None:
            self.kk = Q.shape[1]
        else:
            self.kk = F.shape[1]


def validate_input(Q, F, G):
    """
    Validate input for Q, F and G matrices.
    """

    State(Q, F, G)


#
# Objective function
#
def phi(Q, F, G, accelerate=True, epsilon=1e-50):
    """
    Compute the full objective minus log-likelihood for the model.

    Args:
        Q (array):
            (num_ind, k) array of admixture fractions
        F (array):
            (num_loci, k) array of parental frequencies
        G (array):
            (num_ind, num_loci) array of first allele counts
        accelerate (bool):
            if True, uses Cython implementation
        epsilon (float):
            stabilization parameter added to logarithm arguments (should be
            close to zero)

    Returns:
        The scalar value of the minus log-likelihood.
    """

    if accelerate:
        return _phi(Q, F, G, epsilon)
    fbar = np.dot(Q,  F.T)
    fbar_comp = np.dot(Q, 1 - F.T)
    return -np.sum(G * np.log(fbar + epsilon) +
                   (2 - G) * np.log(fbar_comp + epsilon))


cdef _phi(double[:, :] Q, double[:, :] F, int[:, :] G, double e=1e-50):
    cdef double phi, fbar, fbar_
    cdef double q_ik, f_jk, n_ij
    cdef int i, j, k

    # Check input shapes
    cdef int ii=Q.shape[0], jj=F.shape[0], kk=Q.shape[1]
    if kk != F.shape[1] or ii != G.shape[0] or jj != G.shape[1]:
        raise ValueError('invalid input shapes')

    phi = 0
    for i in range(ii):
        for j in range(jj):
            n_ij = G[i, j]

            if n_ij == 0:
                fbar = e
                for k in range(kk):
                    fbar += Q[i, k] * (1 - F[j, k])
                phi -= 2 * log(fbar)
            elif n_ij == 2:
                fbar = e
                for k in range(kk):
                    fbar += Q[i, k] * F[j, k]
                phi -= 2 * log(fbar)
            else:
                fbar = fbar_ = 0
                for k in range(kk):
                    q_ik = Q[i, k]
                    f_jk = F[j, k]
                    fbar += q_ik * f_jk
                    fbar_ += q_ik * (1 - f_jk)
                phi -= log(fbar * fbar_ + e)

    return phi


#
# Block reductions of objective function
#
def phi_qi(i, Qi, F, G, epsilon=1e-50, accelerate=True):
    """
    Contribution from a single individual i to the objective function.

    The objective function can be factorized as a
    ``sum(func_qi[j] for i in range(num_ind))``.

    Args:
        i: individual index.
        Qi: array of ancestral frequencies for individual.
    """

    if accelerate:
       return _phi_qi(i, Qi, F, G, epsilon)

    mean_f = np.dot(F, Qi)
    mean_F = np.dot(1 - F, Qi)
    gi = G[i]
    return -np.dot(gi, np.log(mean_f)) \
           - np.dot(2 - gi, np.log(mean_F))


cdef double _phi_qi(int i, double[:] Qi,
                    double[:, :] F, int[:, :] G, double e=1e-50) nogil:
    # This implementation is accelerated with 16 partial multiplications before
    # calling log(x). Log is an expensive function and this approach cuts the
    # number of log(x) calls by a factor of 16.

    cdef int j, k, n
    cdef double mean_f_ij, mean_F_ij, f_jk, q_ik, logf, result
    cdef int ii = G.shape[0], jj = G.shape[1], kk = Qi.shape[0]

    result = 0.0
    partmul = 1.0

    for j in range(jj):
        mean_f_ij = mean_F_ij = 0.0

        for k in range(kk):
            f_jk = F[j, k]
            q_ik = Qi[k]
            mean_f_ij += q_ik * f_jk
            mean_F_ij += q_ik * (1 - f_jk)

        n = G[i, j]
        if n == 2:
            result = mean_f_ij * mean_f_ij
        elif n == 0:
            result = mean_F_ij * mean_F_ij
        else:
            result = mean_f_ij * mean_F_ij

        if i % 16 == 0:
            result -= log(partmul + e)
            partmul = 1.0

    if partmul != 1.0:
        result -= log(partmul + e)

    return result


def phi_fj(int j, np.ndarray[double, ndim=1] Fj,
           np.ndarray[double, ndim=2] Q, np.ndarray[int, ndim=2] G,
           epsilon=1e-100, accelerate=True):
    """
    Contribution from a single locus point j to the objective function.

    The objective function can be factorized as a
    ``sum(func_fj[j] for i in range(num_loci))``.

    Args:
        j: locus index.
        fj: array of ancestral frequencies at locus j.
    """

    cdef double result

    if accelerate:
        return _phi_fj(j, &Fj[0], State(Q, None, G), epsilon)

    mean_f = np.dot(Q, Fj)
    mean_F = np.dot(Q, 1 - Fj)
    gj = G[:, j]
    return -np.dot(gj, np.log(mean_f)) \
           - np.dot(2 - gj, np.log(mean_F))


cdef double _phi_fj(int j, double* Fj, State st, double e=1e-100):
    # This implementation is accelerated with 16 partial multiplications before
    # calling log(x). Log is an expensive function and this approach cuts the
    # number of log(x) calls by a factor of 16.

    cdef int i, k, n
    cdef int ii = st.ii, jj = st.jj, kk = st.kk
    cdef double mean_f_ij, mean_F_ij, f_jk, q_ik, result, partmul
    cdef double[:, :] Q = st.Q
    cdef int[:, :] G = st.G

    result = 0.0
    partmul = 1.0

    for i in range(ii):
        mean_f_ij = mean_F_ij = 0
        for k in range(kk):
            f_jk = Fj[k]
            q_ik = Q[i, k]
            mean_f_ij += q_ik * f_jk
            mean_F_ij += q_ik * (1 - f_jk)

        n = G[i, j]
        if n == 2:
            partmul *= mean_f_ij * mean_f_ij
        elif n == 0:
            partmul *= mean_F_ij * mean_F_ij
        else:
            partmul *= mean_f_ij * mean_F_ij

        if i % 16 == 0:
            result -= log(partmul + e)
            partmul = 1.0

    if partmul != 1.0:
        result -= log(partmul + e)

    return result


#
# Block reductions of gradient
#
def grad_qi(int i, np.ndarray[double, ndim=1] Qi, F, G, epsilon=1e-50, accelerate=True):
    """
    Gradient in the f variables for the i-th individual.
    """

    cdef np.ndarray[double, ndim=1] grad

    if accelerate:
        grad = np.empty(Qi.shape[0], dtype=float)
        _grad_qi(i, &Qi[0], State(None, F, G), &grad[0], epsilon)
        return grad

    F_comp = 1 - F
    G_comp = 2 - G

    aux1 = G[..., None] * F
    aux1 /= np.dot(F, Qi)[..., None]
    aux1 = aux1.sum(0)

    aux2 = G_comp[..., None] * F_comp
    aux2 /= np.dot(F_comp, Qi)[..., None]
    aux2 = aux2.sum(0)

    return -(aux1 + aux2)


cdef void _grad_qi(int i, double* Qi, State st, double* grad, double e=1e-50) nogil:
    cdef int j, k, n
    cdef int ii = st.ii, jj = st.jj, kk = st.kk
    cdef double q_ik, f_jk, mean_f_ij, mean_F_ij, result, aux
    cdef double[:, :] F = st.F
    cdef int[:, :] G = st.G

    for j in range(jj):
        mean_f_ij = mean_F_ij = e
        for k_ in range(kk):
            f_jk = F[j, k_]
            q_ik = Qi[k_]
            mean_f_ij += q_ik * f_jk
            mean_F_ij += q_ik * (1 - f_jk)

        n = G[i, j]
        for k in range(kk):
            f_jk = F[j, k]
            grad[k] -= (n * f_jk / mean_f_ij + (2 - n) * (1 - f_jk) / mean_F_ij)


def grad_fj(int j, np.ndarray[double, ndim=1] Fj, Q, G, epsilon=1e-50, accelerate=True):
    """
    Gradient in the f variables for the j-th locus.
    """

    cdef np.ndarray[double, ndim=1] grad

    if accelerate:
        grad = np.empty(Fj.shape[0], dtype=float)
        _grad_fj(j, &Fj[0], State(Q, None, G), &grad[0], epsilon)
        return grad

    gj = G[:, j]
    mean_f = np.dot(Q, Fj) + epsilon
    mean_F = np.dot(Q, 1 - Fj) + epsilon
    aux1 = gj[..., None] * Q / mean_f[..., None]
    aux2 = (2 - gj)[..., None] * Q / mean_F[..., None]
    return -(aux1.sum(0) - aux2.sum(0))


cdef void _grad_fj(int j, double* Fj, State st, double* grad, double e=1e-50) nogil:
    cdef int i, k, n
    cdef int ii = st.ii, jj = st.jj, kk = st.kk
    cdef double q_ik, f_jk, aux, mean_F_ij, mean_f_ij
    cdef double[:, :] Q = st.Q
    cdef int[:, :] G = st.G

    for k in range(kk):
        grad[k] = 0.0

    for i in range(ii):
        mean_f_ij = mean_F_ij = e
        for k in range(kk):
            f_jk = Fj[k]
            q_ik = Q[i, k]
            mean_f_ij += q_ik * f_jk
            mean_F_ij += q_ik * (1 - f_jk)

        n = G[i, j]
        if n == 0:
            aux = -2.0 / mean_F_ij
        elif n == 2:
            aux = 2.0 / mean_f_ij
        elif n == 1:
            aux = 1.0 / mean_f_ij - 1.0 / mean_F_ij
        else:
            aux = 0.0

        for k in range(kk):
            grad[k] -= Q[i, k] * aux


#
# Block reductions of hessian
#
def hess_qi(i, Qi, F, G, epsilon=1e-50, accelerate=True):
    """
    Hessian matrix for the i-th individual at given point.
    """

    if accelerate:
        return _hess_qi(i, Qi, F, G, epsilon)

    F_comp = 1 - F
    gi = G[i]
    gi_comp = 2 - gi
    mean_f = np.dot(F, Qi)
    mean_F = np.dot(F_comp, Qi)

    aux1 = gi[..., None, None] * F[:, :, None] * F[:, None, :]
    aux1 /= (mean_f * mean_f)[..., None, None]

    aux2 = gi_comp[..., None, None] * F_comp[:, :, None] * F_comp[:, None, :]
    aux2 /= (mean_F * mean_F)[..., None, None]

    return aux1.sum(0) + aux2.sum(0)


def hess_fj(int j, np.ndarray[double, ndim=1] Fj, Q, G, epsilon=1e-50, accelerate=True):
    """
    Hessian matrix for the j-th locus at given point.
    """

    cdef np.ndarray[double, ndim=2] hess
    cdef int k

    if accelerate:
        k = Fj.shape[0]
        hess = np.empty((k, k), dtype=float)
        _hess_fj(j, &Fj[0], State(Q, None, G), &hess[0, 0], epsilon)
        return hess

    gj = G[:, j]
    gj_comp = 2 - gj
    Fj_comp = 1 - Fj

    aux1 = gj[..., None, None] * Q[:, :, None] * Q[:, None, :]
    aux1 /= (np.dot(Q, Fj) ** 2)[..., None, None]
    aux1 = aux1.sum(0)

    aux2 = gj_comp[..., None, None] * Q[:, :, None] * Q[:, None, :]
    aux2 /= (np.dot(Q, Fj_comp) ** 2)[..., None, None]
    aux2 = aux2.sum(0)

    return aux1 + aux2


cpdef _hess_qi(int i, double[:] Qi,
              double[:, :] F, int[:, :] G, double e=1e-50):
    cdef int j, k, n
    cdef int ii = G.shape[0], jj = G.shape[1], kk = Qi.shape[0]
    cdef double r, s, q_ik, f_jk, f_jl, mean_f_ij, mean_F_ij, result
    hess = np.zeros((kk, kk), dtype=float)
    cdef double[:, :] hess_view = hess

    with nogil:
        for k in range(kk):
            for l in range(kk):
                result = 0.0

                for j in range(jj):
                    mean_f_ij = mean_F_ij = 0.0
                    for k_ in range(kk):
                        f_jk = F[j, k_]
                        q_ik = Qi[k_]
                        mean_f_ij += q_ik * f_jk
                        mean_F_ij += q_ik * (1 - f_jk)

                    n = G[i, j]
                    f_jk = F[j, k]
                    f_jl = F[j, l]
                    result += n * f_jk * f_jl / (mean_f_ij * mean_f_ij + e)
                    result += (2 - n) * (1 - f_jk) * (1 - f_jl) / (mean_F_ij * mean_F_ij + e)
                hess_view[k, l] = result
    return hess


cdef void _hess_fj(int j, double* Fj, State st, double* hess, double e=1e-50) nogil:
    cdef int i, k, l, n
    cdef int ii = st.ii, jj = st.jj, kk = st.kk
    cdef double q_ik, f_jk, mean_f_ij, mean_F_ij, aux
    cdef double[:, :] Q = st.Q
    cdef int[:, :] G = st.G

    for k in range(kk * kk):
        hess[k] = 0.0

    for i in range(ii):
        mean_f_ij = mean_F_ij = 0
        for k in range(kk):
            f_jk = Fj[k]
            q_ik = Q[i, k]
            mean_f_ij += q_ik * f_jk
            mean_F_ij += q_ik * (1 - f_jk)

        n = G[i, j]
        if n == 0:
            aux = 2.0 / (mean_F_ij * mean_F_ij + e)
        elif n == 2:
            aux = 2.0 / (mean_f_ij * mean_f_ij + e)
        elif n == 1:
            aux = 1.0 / (mean_F_ij * mean_F_ij + e) + 1.0 / (mean_f_ij * mean_f_ij + e)
        else:
            aux = 0.0

        for k in range(kk):
            q_ik = Q[i, k]
            for l in range(k, kk):
                hess[k * kk + l] += q_ik * Q[i, l] * aux

    # Compute symmetric terms in Hessian
    for k in range(kk):
        for l in range(k + 1, kk):
            hess[l * kk + k] += hess[k * kk + l]


#
# Batch evaluators
#
def phi_F(double[:, :] Q, double[:, :] F, int[:, :] G, epislon=1e-50):
    """
    Evaluate the F components of the objective function in batch.

    Args:
        Q, F, G, episilon:
            same meaning as in :func:`phi`.

    Returns:
        A (num_loci) array with each phi_j value.
    """

    cdef State st = State(Q, F, G)
    cdef np.ndarray[double, ndim=1] out = np.empty(st.jj, dtype=float)
    cdef int j, k = st.kk

    for j in range(st.jj):
        out[j] = _phi_fj(j, &F[j, 0], st, epislon)

    return out


def grad_F(double[:, :] Q, double[:, :] F, int[:, :] G, epislon=1e-50):
    """
    Evaluate F gradients in batch.

    Args:
        Q, F, G, episilon:
            same meaning as in :func:`phi`.

    Returns:
        A (num_loci, k) array with each F[j] gradient.
    """

    cdef State st = State(Q, F, G)
    cdef np.ndarray[double, ndim=2] out = np.empty((st.jj, st.kk), dtype=float)
    cdef int j, k = st.kk

    for j in range(st.jj):
        _grad_fj(j, &F[j, 0], st, &out[j, 0], epislon)

    return out


def hess_F(double[:, :] Q, double[:, :] F, int[:, :] G, epislon=1e-50):
    """
    Evaluate F hessians in batch.

    Args:
        Q, F, G, episilon:
            same meaning as in :func:`phi`.

    Returns:
        A (num_loci, k, k) array with the j-th Heassian at each out[j].
    """

    cdef State st = State(Q, F, G)
    cdef int j = st.jj, k = st.kk
    cdef np.ndarray[double, ndim=3] out = np.empty((j, k, k), dtype=float)

    for j in range(st.jj):
        _hess_fj(j, &F[j, 0], st, &out[j, 0, 0], epislon)

    return out


def full_objective_F(double[:, :] Q, double[:, :] F, int[:, :] G, epislon=1e-50):
    """
    Evaluate F func/gradient/hessians in batch.

    Args:
        Q, F, G, epsilon:
            same meaning as in :func:`phi`.

    Returns:
        func, grad and hessian matrices.
    """

    cdef State st = State(Q, F, G)
    cdef int j = st.jj, k = st.kk

    cdef np.ndarray[double, ndim=1] phi = np.empty(j, dtype=float)
    cdef np.ndarray[double, ndim=2] grad = np.empty((j, k), dtype=float)
    cdef np.ndarray[double, ndim=3] hess = np.empty((j, k, k), dtype=float)

    for j in range(st.jj):
        _full_objective_fj(j, &F[j, 0], st, &phi[j], &grad[j, 0], &hess[j, 0, 0], epislon)
    return phi, grad, hess


def derivatives_F(double[:, :] Q, double[:, :] F, int[:, :] G, epislon=1e-50):
    """
    Evaluate F gradient/hessians in batch.

    Args:
        Q, F, G, epsilon:
            same meaning as in :func:`phi`.

    Returns:
        grad and hessian matrices.
    """

    cdef State st = State(Q, F, G)
    cdef int j = st.jj, k = st.kk

    cdef np.ndarray[double, ndim=2] grad = np.empty((j, k), dtype=float)
    cdef np.ndarray[double, ndim=3] hess = np.empty((j, k, k), dtype=float)

    for j in range(st.jj):
        _derivatives_fj(j, &F[j, 0], st, &grad[j, 0], &hess[j, 0, 0], epislon)
    return grad, hess


#
# Shared evaluators
#
def full_objective_fj(int j, double[:] Fj, double[:, :] Q, int[:, :] G, epsilon=1e-50):
    """
    Compute the function value, gradient and hessian at single step.

    This function is more efficient than computing each quantity separately.

    Returns:
        A tuple of (function value, gradient, hessian)
    """

    cdef State st = State(Q, None, G)
    cdef int k = st.kk

    cdef double phi
    cdef np.ndarray[double, ndim=1] grad = np.empty(k, dtype=float)
    cdef np.ndarray[double, ndim=2] hess = np.empty((k, k), dtype=float)

    _full_objective_fj(j, &Fj[0], st, &phi, &grad[0], &hess[0, 0], epsilon)
    return phi, grad, hess


cdef void _full_objective_fj(int j, double* Fj, State st, double* phi, double* grad, double* hess, double e=1e-50) nogil:
    cdef int i, k, l, n
    cdef int ii = st.ii, jj = st.jj, kk = st.kk
    cdef double q_ik, f_jk, mean_f_ij, mean_F_ij, aux_hess, aux_grad, result
    cdef double partmul, aux_phi
    cdef double[:, :] Q = st.Q
    cdef int[:, :] G = st.G
    result = 0.0
    partmul = 1.0

    # Fill hessian and gradient with null components
    for k in range(kk * kk):
        hess[k] = 0.0
    for k in range(kk):
        grad[k] = 0.0

    for i in range(ii):
        mean_f_ij = mean_F_ij = e
        for k in range(kk):
            f_jk = Fj[k]
            q_ik = Q[i, k]
            mean_f_ij += q_ik * f_jk
            mean_F_ij += q_ik * (1 - f_jk)

        n = G[i, j]
        if n == 0:
            aux_phi = mean_F_ij * mean_F_ij
            partmul *= aux_phi
            aux_grad = -2.0 / mean_F_ij
            aux_hess = 2.0 / (aux_phi + e)
        elif n == 2:
            aux_phi = mean_f_ij * mean_f_ij
            partmul *= aux_phi
            aux_grad = 2.0 / mean_f_ij
            aux_hess = 2.0 / (aux_phi + e)
        elif n == 1:
            partmul *= mean_f_ij * mean_F_ij
            aux_grad = 1.0 / mean_f_ij - 1.0 / mean_F_ij
            aux_hess = 1.0 / (mean_F_ij * mean_F_ij + e) + 1.0 / (mean_f_ij * mean_f_ij + e)
        else:
            aux_grad = aux_hess = 0.0

        # Update gradient and hessian
        for k in range(kk):
            q_ik = Q[i, k]
            grad[k] -= q_ik * aux_grad
            for l in range(k, kk):
                hess[k * kk + l] += q_ik * Q[i, l] * aux_hess

        # Feed partial log arguments to function every 16 frames
        if i % 16 == 0:
            result -= log(partmul + e)
            partmul = 1.0

    # Save resulting function value
    if partmul != 1.0:
        result -= log(partmul + e)
    phi[0] = result

    # Compute symmetric terms in Hessian
    for k in range(kk):
        for l in range(k + 1, kk):
            hess[l * kk + k] += hess[k * kk + l]


def derivatives_fj(int j, double[:] Fj, double[:, :] Q, int[:, :] G, epsilon=1e-50):
    """
    Compute the gradient and hessian at single step.

    This function is more efficient than computing each quantity separately.

    Returns:
        A tuple of (function value, gradient, hessian)
    """

    cdef State st = State(Q, None, G)
    cdef int k = st.kk

    cdef np.ndarray[double, ndim=1] grad = np.empty(k, dtype=float)
    cdef np.ndarray[double, ndim=2] hess = np.empty((k, k), dtype=float)

    _derivatives_fj(j, &Fj[0], st, &grad[0], &hess[0, 0], epsilon)
    return grad, hess


cdef void _derivatives_fj(int j, double* Fj, State st, double* grad, double* hess, double e=1e-50) nogil:
    cdef int i, k, l, n
    cdef int ii = st.ii, jj = st.jj, kk = st.kk
    cdef double q_ik, f_jk, mean_f_ij, mean_F_ij, aux_hess, aux_grad
    cdef double[:, :] Q = st.Q
    cdef int[:, :] G = st.G

    # Fill hessian and gradient with null components
    for k in range(kk * kk):
        hess[k] = 0.0
    for k in range(kk):
        grad[k] = 0.0

    for i in range(ii):
        mean_f_ij = mean_F_ij = e
        for k in range(kk):
            f_jk = Fj[k]
            q_ik = Q[i, k]
            mean_f_ij += q_ik * f_jk
            mean_F_ij += q_ik * (1 - f_jk)

        n = G[i, j]
        if n == 0:
            aux_grad = -2.0 / mean_F_ij
            aux_hess = 2.0 / (mean_F_ij * mean_F_ij + e)
        elif n == 2:
            aux_grad = 2.0 / mean_f_ij
            aux_hess = 2.0 / (mean_f_ij * mean_f_ij + e)
        elif n == 1:
            aux_grad = 1.0 / mean_f_ij - 1.0 / mean_F_ij
            aux_hess = 1.0 / (mean_F_ij * mean_F_ij + e) + 1.0 / (mean_f_ij * mean_f_ij + e)
        else:
            aux_grad = aux_hess = 0.0

        # Update gradient and hessian
        for k in range(kk):
            q_ik = Q[i, k]
            grad[k] -= q_ik * aux_grad
            for l in range(k, kk):
                hess[k * kk + l] += q_ik * Q[i, l] * aux_hess

    # Compute symmetric terms in Hessian
    for k in range(kk):
        for l in range(k + 1, kk):
            hess[l * kk + k] += hess[k * kk + l]
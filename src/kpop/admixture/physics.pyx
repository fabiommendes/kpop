cimport cython
import numpy as np
cimport numpy as np
from cython.view cimport array as cvarray


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.overflowcheck(False)
def _force_q(double[:, :] Q, double[:, :] F, int[:, :] G, double e=1e-50):
    cdef double phi = 0.0, fbar_ij, fbar_ij_
    cdef double q_ik, f_jk, n_ij, force_qi
    cdef int i, j, k, k_
    cdef int ii=Q.shape[0], jj=F.shape[0], kk=Q.shape[1]
    result = np.empty((ii, kk), dtype=float)
    cdef double[:, :] force = result
    cdef double[:, :] fbar = cvarray(shape=(ii, jj), itemsize=sizeof(double), format='d')
    cdef double[:, :] fbar_ = cvarray(shape=(ii, jj), itemsize=sizeof(double), format='d')

    # Cache fbar and fbar_
    for i in range(ii):
        for j in range(jj):
            fbar_ij = fbar_ij_ = e
            for k in range(kk):
                q_ik = Q[i, k]
                f_jk = F[j, k]
                fbar_ij += q_ik * f_jk
                fbar_ij_ += q_ik * (1 - f_jk)
            fbar[i, j] = fbar_ij
            fbar_[i, j] = fbar_ij_

    for k in range(kk):
        for i in range(ii):
            force_qi = 0
            for j in range(jj):
                n_ij = G[i, j]
                f_jk = F[j, k]
                force_qi += n_ij * f_jk / fbar[i, j] + (2 - n_ij) * (1 - f_jk) / fbar_[i, j]
            force[i, k] = force_qi
    return result


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.overflowcheck(False)
def _force_f(double[:, :] Q, double[:, :] F, int[:, :] G, double e=1e-50):
    cdef double phi = 0.0, fbar_ij, fbar_ij_
    cdef double q_ik, f_jk, n_ij, force_fj
    cdef int i, j, k
    cdef int ii=Q.shape[0], jj=F.shape[0], kk=Q.shape[1]
    result = np.empty((jj, kk), dtype=float)
    cdef double[:, :] force = result
    cdef double[:, :] fbar = cvarray(shape=(ii, jj), itemsize=sizeof(double), format='d')
    cdef double[:, :] fbar_ = cvarray(shape=(ii, jj), itemsize=sizeof(double), format='d')

    # Cache fbar and fbar_
    for i in range(ii):
        for j in range(jj):
            fbar_ij = fbar_ij_ = e
            for k in range(kk):
                q_ik = Q[i, k]
                f_jk = F[j, k]
                fbar_ij += q_ik * f_jk
                fbar_ij_ += q_ik * (1 - f_jk)
            fbar[i, j] = fbar_ij
            fbar_[i, j] = fbar_ij_

    for j in range(jj):
        force_fj = 0
        for i in range(ii):
            n_ij = G[i, j]
            q_ik = Q[i, k]
            force_fj += q_ik * (n_ij / fbar[i, j] - (2 - n_ij) / fbar_[i, j])
        force[j] = force_fj

    return result

@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.overflowcheck(False)
def _forces(double[:, :] Q, double[:, :] F, int[:, :] G, double e=1e-50):
    cdef double phi = 0.0, fbar_ij, fbar_ij_
    cdef double q_ik, f_jk, n_ij, force_qi
    cdef int i, j, k, k_
    cdef int ii=Q.shape[0], jj=F.shape[0], kk=Q.shape[1]
    result = np.empty((ii, kk), dtype=float), np.empty((jj, kk), dtype=float)
    cdef double[:, :] force_q = result[0]
    cdef double[:, :] force_f = result[1]
    cdef double[:, :] chi = cvarray(shape=(ii, jj), itemsize=sizeof(double), format='d')
    cdef double[:, :] chi_ = cvarray(shape=(ii, jj), itemsize=sizeof(double), format='d')

    # Cache chi and chi_
    for i in range(ii):
        for j in range(jj):
            fbar_ij = fbar_ij_ = e
            for k in range(kk):
                q_ik = Q[i, k]
                f_jk = F[j, k]
                fbar_ij += q_ik * f_jk
                fbar_ij_ += q_ik * (1 - f_jk)
            chi[i, j] = G[i, j] / fbar_ij
            chi_[i, j] = (2 - G[i, j]) / fbar_ij_

    # Compute force_q
    for k in range(kk):
        for i in range(ii):
            force_qi = 0
            for j in range(jj):
                f_jk = F[j, k]
                force_qi += chi[i, j] * f_jk + chi_[i, j] * (1 - f_jk)
            force_q[i, k] = force_qi

    # Compute force_f
    for k in range(kk):
        for j in range(jj):
            force_fj = 0
            for i in range(ii):
                force_fj += Q[i, k] * (chi[i, j] - chi_[i, j])
            force_f[j, k] = force_fj

    return result

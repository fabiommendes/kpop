#cython: cdivision=True, wraparound=False, initializedcheck=False, overflowcheck=False

cimport cython
cimport numpy as np
import numpy as np
from numpy.linalg import inv, solve


cdef void _sinv2(double* M, double* out) nogil:
    """
    Inverts a symmetric 2x2 matrix using Crammer's method and save it in out.

    Since matrix is symmetric, it can be in either FORTRAN or C order.
    """

    cdef double a = M[0], b = M[1], c = M[3]
    cdef double Z = 1.0 / (a * c - b * b)

    out[0] = c * Z
    out[1] = out[2] = -b * Z
    out[3] = a * Z


cdef void _sinv3(double* M, double* out) nogil:
    """
    Inverts a symmetric 3x3 matrix using Crammer's method and save it in out.

    Since matrix is symmetric, it can be in either FORTRAN or C order.
    """

    cdef double a = M[0], b = M[1], c = M[2], d = M[4], e = M[5], f = M[8]
    cdef double A = d*f - e*e, B = e*c - b*f, C = b*e - c*d
    cdef double Z = 1.0 / (a * A +  b * B + c * C)

    out[0] = A * Z
    out[1] = out[3] = B * Z
    out[2] = out[6] = C * Z
    out[4] = (a*f - c*c) * Z
    out[5] = out[7] = -(a*e - b*c) * Z
    out[8] = (a*d - b*b) * Z


cpdef sinv2(double[:, :] M):
    """
    Return the inverse of a 2x2 symmetric matrix.
    """

    if M.shape[0] != 2 or M.shape[1] != 2:
        raise ValueError('input must be a 2x2 matrix!')
        
    out = np.empty((2, 2), dtype=float)
    cdef double[:, :] out_view = out
    _sinv2(&M[0, 0], &out_view[0, 0])
    return out


cpdef sinv3(double[:, :] M):
    """
    Return the inverse of a 3x3 symmetric matrix.
    """

    if M.shape[0] != 3 or M.shape[1] != 3:
        raise ValueError('input must be a 3x3 matrix!')
    
    out = np.empty((3, 3), dtype=float)
    cdef double[:, :] out_view = out
    _sinv3(&M[0, 0], &out_view[0, 0])
    return out


def sinv(double[:, :] M):
    """
    Invert symmetric matrix using a fast method.
    """
    
    cdef int m = M.shape[0], n = M.shape[1]
    if m != n:
        raise ValueError('matrix must be square!')
    if m == 2:
        return sinv2(M)
    elif m == 3:
        return sinv3(M)
    else:
        return inv(M)
    
    
def sinv_batch(data):
    """
    Invert all symmetric matrices in the data array.
    
    Args:
        data (size, n, n):
            A array with all matrices. Each (n, n) matrix must be symmetric. 
    """
    
    cdef double [:, :, :] in_ = data
    cdef int size = in_.shape[0], m = in_.shape[1], n = in_.shape[2]
    if n != m:
        raise ValueError('matrix must be square!')
    
    result = np.empty((size, n, n), dtype=float)
    cdef double[:, :, :] out = result
    
    cdef int k
    if n == 2:
        for k in range(size):
            _sinv2(&in_[k, 0, 0], &out[k, 0, 0])
    elif n == 3:
        for k in range(size):
            _sinv3(&in_[k, 0, 0], &out[k, 0, 0])
    else:
        for k in range(size):
            result[k] = inv(data[k])
    return result

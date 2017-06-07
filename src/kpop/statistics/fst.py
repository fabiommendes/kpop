import numpy as np


def biallelic_pairwise_fst(freqs: np.ndarray,
                           h_freqs: np.ndarray,
                           sizes: np.ndarray) -> np.ndarray:
    """
    Return a matrix of pairwise Fst components from frequency data. The Fst
    coefficients are computed using the method introduced by B. S. Weir and C.
    Clark Cockerham in *F-Statistics for the Analysis of Population
    Structure*, Evolution, Vol. 38, No. 6 (Nov., 1984), pp. 1358-1370".

    Args:
        freqs:
            A 2D array of frequencies. freqs[k, j] corresponds to the
            frequency of the first allele in the j-th locus for the k-th
            population.
        h_freqs:
            A 2D array of heterozygote frequencies. Represents the probability
            that the j-th locus in the k-th population is an heterozygote.
        sizes:
            List of sizes for each population.

    Returns:
        A squared matrix of fst coefficients.
    """
    freqs = np.asarray(freqs)
    h_freqs = np.asarray(h_freqs)
    k, j = freqs.shape
    sizes = np.asarray(sizes)
    fst_matrix = np.zeros((k, k), dtype=float)

    if freqs.shape != h_freqs.shape:
        shapes = freqs.shape , h_freqs.shape
        raise ValueError('inconsistent frequency shapes: %s' % (shapes,))
    if freqs.shape[0] != len(sizes):
        shapes = freqs.shape, sizes.shape
        raise ValueError('inconsistent shapes: %s and %s' % (shapes))

    for k1 in range(k):
        for k2 in range(k1 + 1, k):
            f1, f2 = freqs[k1], freqs[k2],
            h1, h2 = h_freqs[k1], h_freqs[k2],
            n1, n2 = sizes[k1], sizes[k2]
            fst = biallelic_fst(f1, f2, h1, h2, [n1, n2])
            fst_matrix[k1, k2] = fst_matrix[k2, k1] = fst

    return fst_matrix


def biallelic_fst(freqs1, freqs2, h_freqs1, h_freqs2, sizes):
    """
    Compute the Fst statistics for a pair of populations.

    Args:
        freqs1, freqs2:
            Arrays of frequencies for the first allele in each population.
        h_freqs1, h_freqs2:
            Arrays of heterozygote frequencies in each population.
        sizes:
            Size of each sub-population.
    Returns:

    """

    f1, f2 = map(np.asarray, (freqs1, freqs2))
    h1, h2 = map(np.asarray, (h_freqs1, h_freqs2))
    n1, n2 = sizes

    # Compute auxiliary variables
    n = (n1 + n2) / 2
    nc = (2 * n - (n1**2 + n2**2) / (2 * n))
    p = (n1 * f1 + n2 * f2) / (2 * n)
    s2 = (n1 * (f1 - p)**2 + n2 * (f2 - p)**2) / n
    h = (n1 * h1 + n2 * h2) / (2 * n)

    # Compute a, b, c
    a = n / nc * (
        s2 - 1 / (n- 1) * (p * (1 - p) - s2 / 2- h / 4)
    )
    b = n/ (n- 1) * (
        p * (1 - p) - s2 / 2 - (2 * n - 1) / (4 * n) * h
    )
    c = h / 2

    # Compute the Fst
    e = 1e-50
    return a.sum() / ((a + b + c).sum() + e)


def repr_pairwise_matrix(matrix):
    """
    Return a string representation of a square matrix that represents pairwise
    information as a triangular matrix displaying each pair.

    This function ignores the lower triangular part of the matrix. Usually the
    input matrix should be symmetric, but some applications may simply choose to
    omit the redundant lower triangular part.

    Examples:
        >>> M = [[0, 2, 3],
        ...      [2, 0, 4],
        ...      [3, 4, 0]]
        >>> print(repr_pairwise_matrix(M))
            0  1
        1 - 2,
        2 - 3, 4,

    Args:
        matrix: n x n matrix data

    Returns:
        A string representing the input matrix.
    """

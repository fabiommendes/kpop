import numpy as np
from lazyutils import lazy
from sklearn import manifold

from .population_base import PopulationBase
NOT_GIVEN = object()


class Projection:
    """
    Implements all projection methods (such as PCA) applied to Kpop populations.
    """

    _data = property(lambda _: np.array(_.population))

    def __init__(self, population):
        self.population = population

    def as_data(self, method='flatten', norm=NOT_GIVEN):
        """
        Convert genotype data into a numpy array. This is a basic pre-processing
        step in many dimensionality reduction algorithms.

        Genotypes are categorical data and usually it doesn't make sense to
        treat the integer encoding used in kpop as ordinal data (there is
        no ordering implied when treating say, allele 1 vs allele 2 vs allele
        3).

        There are a few basic strategies to convert categorical data in a 
        genotype to something that can be handled by machine-learning 
        algorithms.

        1. The 'count' method simply counts the number of #1 alleles in each
           locus and creates a matrix C(n,j) of the number of #1 alleles for
           individual n in location j. This only makes sense when working with
           biallelic data.
        2. The 'flatten' method shuffles the alleles at each loci, flatten
           it, and creates a matrix C(n, j) with (N x 2J) elements. The rest
           of the algorithm proceeds identically. This pre-processing must be
           used with care for non-biallelic data.
        3. We still don't implemented methods that can be

        Args:
            method : 'count', 'flatten'
                Genetic discrimination method. See description above.
            norm : dfs

        Returns:
            An ndarray with transformed data.
        """

        data = as_data(self.population, method)
        if norm is NOT_GIVEN:
            norm = {'count': 'snp', 'flatten': 'std'}.get(method)
            data = whiten(data, norm)
        elif norm is not None:
            data = whiten(data, norm)
        return data
    
    def pca(self, k=2, *, method='count', norm=True):
        """
        Principal component analysis.

        PCA helps classifying individuals by reducing the problem's
        dimensionality from the number of loci to only the first K of principal
        axis. Usually this simple re-parametrization is able to detect population
        structure between data, i.e., the individuals of each sub-population are
        very well separated in the "principal axes" coordinates.

        The principal axes are chosen from a singular value decomposition (SVD)
        of the population matrix. This can be done in 2 different ways:


        The matrix C is subjected to a SVD, and each individual is classified in
        terms of their coordinates in the K most significant principal axis.

        Args:
            k (int):
                Number of principal components to consider
            method:
                See the :method:`as_data` method.
            norm (bool):
                If True (default), normalizes each j component of the C(n, j)
                matrix according to a measure of genetic drift. This procedure
                is described at Patterson et. al., "Population Structure and
                Eigenanalysis" and is recommended for SNPs subject to genetic
                drift. It is not recommended to use this normalization to
                micro-satellite data.

        Returns:
            A (size x k) matrix with the components of each individual in the
            principal axis.

        Example:

            Consider a random synthetic population with two sub-populations with
            10 elements each. Each individual has 200 alleles.

            >>> popA = Population.make_random(10, 200, label='A')
            >>> popB = Population.make_random(10, 200, label='B')

            Usually the the principal axis alone will be enough to classify
            these individuals. Since the mean is zero, individuals of one
            population will have negative coordinates and the other population
            will have positive coordinates.

            >>> pop = popA + popB
            >>> coords = pop.projection.pca()
            >>> sign_p1 = coords[10:] > 0
            >>> sign_p2 = coords[:10] > 0
            >>> (sign_p1 != sign_p2).all()
            True
        """
        norm = {True: 'snp', False: 'mean'}.get(norm, norm)
        data = self.as_data(method, norm=norm)

        # Compute the singular value decomposition of the rescaled matrix and
        # return the projections of each individual in this matrix
        _, _, pca_coords = np.linalg.svd(data, False)
        pca_coords_k = pca_coords[:k]
        return np.dot(data, pca_coords_k.T)

    def tsne(self, k=2, *, method='count', norm=True, **kwargs):
        """
        t-distributed stochastic neighbor embedding.

        It is a dimensionality reduction method that tries to preserve 
        similarities between different features. Differently from PCA, which is
        a simple rotation, t-SNE tends to preserve similarities and lend itself
        for better clustering performance in the reduced coordinates.

        Args:
            Same as in :method:`pca`.

        Notes:
            See: http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html#sklearn.manifold.TSNE
        """

        norm = {True: 'snp', False: 'mean'}.get(norm, norm)
        data = self.as_data(method, norm=norm)
        return sklearn_result(manifold.TSNE(k, **kwargs), data, **kwargs)

    def mds(self, k=2, *, method='count', norm=False, **kwargs):
        """
        Multidimensional scaling.

        This is a low dimensional representation of data that tries to preserve
        the same distance matrix as in the original high dimensionality space.
        """

        norm = {True: 'snp', False: 'mean'}.get(norm, norm)
        data = self.as_data(method, norm=norm)
        return sklearn_result(manifold.MDS(k, **kwargs), data, **kwargs)


# Patch Population class
PopulationBase.projection = lazy(lambda _: Projection(_))




#
# Utility functions
#
def as_data(population, method):
    """
    Convert raw genotype data to a format usable in dimensionality reduction
    algorithms.
    """

    if method == 'count':
        data = (np.array(population) == 1).sum(axis=2)
    elif method == 'flatten':
        pop = population.shuffled_loci()
        data = [ind.flatten() for ind in pop]
        data = np.array(data, dtype=pop.dtype)
    else:
        raise ValueError('invalid method: %r' % method)
    
    return data


def whiten(data, method):
    "Remove bias and re-normalize, if desired"

    mu = data.mean(axis=0)
    if method == 'snp':
        p = mu / 2
        norm = np.sqrt(p * (1 - p))
    elif method == 'std':
        norm = data.std(axis=0)
    elif method == 'mean':
        norm = 1.0
    else:
        raise ValueError('invalid normalization method: %r' % method)

    # Prevents problems with mean 0 or 1.
    norm = np.where(norm, norm, 1)
    return (data - mu) / norm

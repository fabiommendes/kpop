import numpy as np
from sklearn import manifold, decomposition


from ..result import result, transform_result

NOT_GIVEN = object()


class projection:
    """
    Implements all projection methods (such as PCA) applied to Kpop populations.
    """

    _data = property(lambda _: np.array(_.population))
    _methods = {'pca', 'tsne', 'lle', 'mds', 'isomap', 'spectral'}

    def __init__(self, population):
        self._population = population

    def as_data(self, method='count'):
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

           The default strategy is to normalize each component of the C(n, j)
           matrix according to a measure of genetic drift. This procedure is
           described at Patterson et. al., "Population Structure and
           Eigenanalysis" and is recommended for SNPs subject to genetic
           drift. It is not recommended to use this normalization for
           micro-satellite data.
        2. The 'flatten' method shuffles the alleles at each loci, flatten
           it, and creates a matrix C(n, j) with (N x 2J) elements. The rest
           of the algorithm proceeds identically. This pre-processing must be
           used with care for non-biallelic data.

        Args:
            method : 'count', 'flatten'
                Conversion method. See description above.

        Returns:
            An ndarray with transformed data.
        """

        data = as_raw_data(self._population, method)
        if method == 'count':
            data = whiten(data, 'snp')
        elif method == 'flatten':
            data = whiten(data, 'std')
        return data

    def project(self, which, k=2, **kwargs):
        """
        Alias to population.project(which, k, ...)
        """

        return self(which, k, **kwargs)

    def __call__(self, which='pca', k=2, **kwargs):
        """
        Embed population data into a low dimensionality sub-space using one of
        the flowing methods:

            * PCA: a classic method. Performs a simple rotation that
            removes correlations from data and keeps the k coordinates with the
            highest variance.
            * t-SNE: ...

        For more details, see the corresponding method name in the class
        documentation.
        """
        if callable(which):
            norm = kwargs.pop('norm', True)
            method = kwargs.pop('method', 'count')
            norm = {True: 'snp', False: 'mean'}.get(norm, norm)
            data = self.as_data(method)
            return transform_result(which, data, k, **kwargs)

        elif isinstance(which, str):
            which_ = which.lower().replace('-', '')
            method = getattr(self, which_)
            return method(k, **kwargs)

        raise ValueError('invalid method: %r' % which)

    def pca(self, k=2, *, method='count', svd=False, **kwargs):
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
            method ('count' or 'flatten'):
                See the :method:`as_data` method.
            svd (bool):
                If True, performs the PCA using an explicit single value
                decomposition. If SVD is true and k=None, then it performs the
                full decomposition of data. It is generally slower than the
                standard PCA for low k, but can be faster for the full
                decomposition.

            Extra keyword arguments are passed to the
            :cls:`sklearn.decomposition.PCA` constructor before creating the
            transformation object.


        Returns:
            A (size x k) matrix with the components of each individual in the
            principal axis.

        See Also:
            http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html

        Example:

            Consider a random synthetic population with two sub-populations with
            10 elements each. Each individual has 200 alleles.

            >>> popA = Population.random(10, 200, label='A')
            >>> popB = Population.random(10, 200, label='B')

            Usually the the principal axis alone will be enough to classify
            these individuals. Since the mean is zero, individuals of one
            population will have negative coordinates and the other population
            will have positive coordinates.

            >>> pop = popA + popB
            >>> coords = pop.scatter_coords.pca()
            >>> sign_p1 = coords[10:] > 0
            >>> sign_p2 = coords[:10] > 0
            >>> (sign_p1 != sign_p2).all()
            True
        """
        # Compute the singular value decomposition of the rescaled matrix and
        # return the projections of each individual in this matrix
        if svd:
            data = self.as_data(method)

            _, _, pca_coords = np.linalg.svd(data, False)
            if k is None or k == 0:
                pca_coords_k = pca_coords[:]
            else:
                pca_coords_k = pca_coords[:k]
            return result(np.dot(data, pca_coords_k.T))
        else:
            pca = decomposition.PCA
            return sklearn_manifold(pca, self, k, method, **kwargs)

    def tsne(self, k=2, *, method='count', max_dimensions=15,
             perplexity=30, **kwargs):
        """
        t-distributed stochastic neighbor embedding.

        It is a dimensionality reduction method that tries to preserve 
        similarities between different features. Differently from PCA, which is
        a simple rotation, t-SNE tends to preserve similarities and lend itself
        for better clustering performance in the reduced coordinates.

        Args:
            k (int):
                Number of principal components to consider
            method ('count' or 'flatten'):
                See the :method:`as_data` method.


            Extra keyword arguments are passed to the
            :cls:`sklearn.manifold.TSNE` constructor before creating the
            transformation object.


        Returns:
            A (size x k) matrix with the components of each individual in the
            principal axis.

        See Also:
            http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html#sklearn.manifold.TSNE
        """
        data = self.as_data(method)
        if max_dimensions is not None and self._population.size > max_dimensions:
            data = transform_result(decomposition.PCA, data, max_dimensions)

        # Extra keyword arguments
        kwargs['perplexity'] = perplexity

        return transform_result(manifold.TSNE, data, k, **kwargs)

    def mds(self, k=2, *, method='count', **kwargs):
        """
        Multidimensional scaling.

        This is a low dimensional representation of data that tries to preserve
        the same distance matrix as in the original high dimensionality space.
        """
        return sklearn_manifold(manifold.MDS, self, k, method, **kwargs)

    def isomap(self, k=2, *, method='count', n_neighbors=None, **kwargs):
        """
        Isomap
        """
        if n_neighbors is None:
            n_neighbors = min(int(self._population.size / 3), 5)
        kwargs['n_neighbors'] = n_neighbors

        def isomap(k, **kwargs):
            kwargs['n_components'] = k
            return manifold.Isomap(**kwargs)

        return sklearn_manifold(isomap, self, k, method, **kwargs)

    def lle(self, k=2, *, method='count', n_neighbors=10, **kwargs):
        # if n_neighbors is None:
        #     n_neighbors = min(int(self._population.size / 3), 5)
        kwargs['n_neighbors'] = n_neighbors
        kwargs.setdefault('reg', 1.0)

        def lle(k, **kwargs):
            kwargs['n_components'] = k
            return manifold.LocallyLinearEmbedding(**kwargs)

        return sklearn_manifold(lle, self, k, method, **kwargs)

    # def ltsa(self, k=2, *, method='count', **kwargs):
    #     return sklearn_manifold(manifold.LTSA, self, k, method, **kwargs)
    #
    # def hessian_lle(self, k=2, *, method='count', **kwargs):
    #     return sklearn_manifold(manifold.HLLE, self, k, method, **kwargs)

    def spectral(self, k=2, *, method='count', **kwargs):
        kwargs.setdefault('affinity', 'rbf')
        return sklearn_manifold(manifold.SpectralEmbedding, self, k, method, **kwargs)


#
# Utility functions
#
def sklearn_manifold(_sk_method, proj: projection, k: int, _method:
str,
                     **kwargs):
    """
    A projection method from a scikit method for manifold learning.
    """
    data = proj.as_data(_method)
    return transform_result(_sk_method, data, k, **kwargs)


def as_raw_data(population, method):
    """
    Convert raw genotype data to a format usable in a dimensionality reduction
    algorithms. It does not normalize the results.
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

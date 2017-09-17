import numpy as np
from sklearn import manifold, decomposition

from .attr import Attr
from ..result import transform_result
from ..utils import is_transformer

NOT_GIVEN = object()


class Projection(Attr):
    """
    Implements the population.projection attribute with all projection
    methods (such as PCA) applied to Kpop populations.
    """

    _methods = {'pca', 'tsne', 'lle', 'mds', 'isomap', 'spectral',
                'kernel_pca', 'nmf', 'ica'}

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
        if is_transformer(which):
            return self.sklearn(which, k=k, **kwargs)
        elif callable(which):
            raise NotImplementedError('do not accept function transformers')
        elif isinstance(which, str):
            which_ = which.lower().replace('-', '_')
            if which_ in self._methods:
                method = getattr(self, which_)
                return method(k, **kwargs)

        raise ValueError('invalid method: %r' % which)

    def _as_array(self, data='count'):
        return self._population.as_array(data)

    def project(self, which, k=2, **kwargs):
        """
        Alias to population.project(which, k, ...)
        """

        return self(which, k, **kwargs)

    def pca(self, k=2, data='count-unity', **kwargs):
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
            data:
                Data conversion method. See the
                :method:`kpop.Population.as_data` method.

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

            >>> popA = Population.random(10, 200, id='A')
            >>> popB = Population.random(10, 200, id='B')

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
        if 'pca' in kwargs:
            raise TypeError('invalid argument: pca')

        pca = decomposition.PCA
        data = self._as_array(data)
        return transform_result(pca, data, n_components=k, **kwargs)

    def mds(self, k=2, *, data='count-unity', **kwargs):
        """
        Multidimensional scaling.

        This is a low dimensional representation of data that tries to preserve
        the same distance matrix as in the original high dimensionality space.
        """
        return self.sklearn(manifold.MDS, k, data, **kwargs)

    def kernel_pca(self, k=2, *, data='count-unity', **kwargs):
        """
        Kernel PCA transformation of data.

        Args:
            k:
                Number of dimensions of the reduced space.
            data:
                Data conversion method. See the
                :method:`kpop.Population.as_data` method.

        """
        kernel_pca = decomposition.KernelPCA
        kwargs = with_defaults(
            kwargs,
            kernel='sigmoid',
        )
        return self.sklearn(kernel_pca, k, data, **kwargs)

    def tsne(self, k=2, *, data='count-unity', **kwargs):
        """
        t-distributed stochastic neighbor embedding.

        It is a dimensionality reduction method that tries to preserve
        similarities between different features. Differently from PCA, which is
        a simple rotation, t-SNE tends to preserve similarities and lend itself
        for better clustering performance in the reduced coordinates.

        Args:
            k (int):
                Number of principal components to consider
            data:
                Data conversion method. See the
                :method:`kpop.Population.as_data` method.
            pca (int):
                Apply PCA targeting the given number of dimensions before
                applying TSNE. This accelerates computation, but also can lead
                to better results since t-SNE tends to behave somewhat badly
                in very high dimension spaces.

        Notes:
            Extra keyword arguments are passed to the
            :cls:`sklearn.manifold.TSNE` constructor before creating the
            transformation object.

        Returns:
            A (size x k) matrix with the components of each individual in the
            principal axis.

        See Also:
            http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html#sklearn.manifold.TSNE
        """
        pca = manifold.TSNE
        pca_dims = None if self._population.num_loci < 15 else 15
        kwargs = with_defaults(kwargs, perplexity=30, pca=pca_dims)
        return self.sklearn(pca, k, data, **kwargs)

    def ica(self, k=2, *, data='count-unity', **kwargs):
        ica = decomposition.FastICA
        return self.sklearn(ica, k, data, **kwargs)

    def factor(self, k=2, *, data='count-unity', **kwargs):
        factor = decomposition.FactorAnalysis
        return self.sklearn(factor, k, data, **kwargs)

    def lda_projection(self, k=2, n_populations=5, **kwargs):
        """
        Compute the admixture distribution assuming n_populations parental
        populations and project the resulting Q-matrix to k dimensions using
        PCA. This method can be used to visualize data in a triangle (for
        n_populations=3) or in a flattened simplex.
        """
        if k >= n_populations:
            raise ValueError('k must be greater than n_populations')

        lda = decomposition.LatentDirichletAllocation
        data = self._as_array('count')

        if is_sklearn_version_gt(19):
            kwargs.update(n_components=n_populations)
        else:
            kwargs.update(n_topics=n_populations)

        q_matrix = transform_result(lda, data, learning_method='batch', **kwargs)
        return decomposition.PCA(k).fit_transform(q_matrix)

    def nmf(self, k=2, *, data='count', **kwargs):
        if 'pca' in kwargs:
            raise TypeError('nmf does not support PCA pre-processing')
        nmf = decomposition.NMF
        return self.sklearn(nmf, k, data, **kwargs)

    def isomap(self, k=2, *, data='count-unity', n_neighbors=None, **kwargs):
        """
        Isomap Embedding.

        Non-linear dimensionality reduction through Isometric Mapping.

        Args:
            n_neighbors (int):
                Number of neighbors to consider for each point.

        Returns:

        """
        isomap = manifold.Isomap
        kwargs = with_defaults(
            kwargs,
            n_neighbors=n_neighbors or
                        min(n_neighbors or 25, self._population.size - 1),
        )
        return self.sklearn(isomap, k, data, **kwargs)

    def lle(self, k=2, *, data='count', n_neighbors=None, reg=1.0, **kwargs):
        """
        Locally linear embedding.

        This method is implemented here for sake of completeness, but it
        seems to perform very poorly in genotype data.

        See Also:
            :cls:`sklearn.manifold.LocallyLinearEmbedding`
        """
        lle = manifold.LocallyLinearEmbedding
        kwargs_ = with_defaults(
            kwargs,
            n_neighbors=n_neighbors or
                        min(n_neighbors or 25, self._population.size - 1),
        )
        return self.sklearn(lle, k, data, **kwargs_)

    def spectral(self, k=2, *, data='count-unity', **kwargs):
        spectral = manifold.SpectralEmbedding
        kwargs = with_defaults(kwargs, affinity='rbf')
        return self.sklearn(spectral, k, data, **kwargs)

    def sklearn(self, transformer, k=2, data='count-unity', pca=None,
                **kwargs):
        """
        Apply a scikit-learn dimensionality reduction method to data. Users
        must pass the transformer class associated with the method (e.g.
        sklearn.manifold.TSNE).

        Args:
            transformer (callable):
                Any callable that return a transformer object.
            k (int):
                Number of dimensions of target space.
            data ('count-?' or 'flat-?'):
                See the :method:`as_array` method. Defaults to 'count-unity'.
            pca (int or None):
                If given, preprocess data with PCA and reduce it to the given
                number of dimensions. Many algoritms perform better (both in
                terms of speed and quality) if data is reduced to a medium
                dimensionality space before treatment.
            """
        if pca is None:
            arr_data = self._as_array(data)
        else:
            arr_data = self.pca(pca, data)
        return transform_result(transformer, arr_data, n_components=k, **kwargs)


#
# Utility functions
#
def with_defaults(dic, **kwargs):
    """
    Return a copy of dictionary where the positional keyword arguments are
    inserted as default values.
    """
    dic = dict(dic)
    for k, v in kwargs.items():
        dic.setdefault(k, v)
    return dic


def is_sklearn_version_gt(version=19):
    """
    Test Scikit Learn minor version.
    """
    import sklearn
    return int(sklearn.__version__.split('.')[1]) >= version

import numpy as np

from .attr import Attr
from ..plot import admixture_scatter, admixture_bars
from ..plot.utils import group_individuals, _colors
from ..utils import lazy_module

plt = lazy_module('matplotlib.pyplot')


class Plot(Attr):
    """
    Implements the Population.plots attribute.
    """

    _populations = property(lambda self: self._population.populations)

    def _pop_sizes(self):
        return [len(pop) for pop in self._populations]

    def _pop_labels(self, population=None):
        populations = (population or self._population).populations
        labels = []
        for i, pop in enumerate(populations):
            label = pop.id if pop.id else 'pop-%s' % i
            labels.append(label)
        return labels

    def freqs(self, merge=False, title='Alelle frequencies',
              colors=None, ylabel=None, scatter=False,
              show=True, sorted=False):
        """
        Plot allele frequencies.

        Args:
            title:
                Optional plot title.
            colors:
                A colormap or a list of colors for each allele.
            ylabel:
                Label for the y-axis.
            scatter:
                If True, produces a scatter plot of frequencies. This only
                makes sense if the number of alleles is greater than 2.
            merge (bool):
                If True, treats a MultiPopulation as a single population.
                Otherwise separate data for each sub-population in the graph.
            show (bool):
                If True (default), immediately displays graph on screen.
            sorted (bool):
                If True, reorganize allele frequencies so that the most likely
                frequency is always associated with the first allele. This might
                help visualization of biallelic data.
        """
        if merge:
            populations = [self._population]
        else:
            populations = self._populations

        data = np.concatenate([pop.freqs_matrix for pop in populations], 0)
        pop_sizes = [pop.num_loci for pop in populations]
        pop_labels = [pop.id or i for i, pop in enumerate(populations)]
        alleles = data.shape[1]
        labels = ['allele-%s' % (i + 1) for i in range(alleles)]
        if sorted:
            data *= -1
            data.sort(axis=-1)
            data *= -1

        if scatter:
            if alleles <= 2:
                raise ValueError('scatter plots requires at least 3 alleles.')
            ax = admixture_scatter(
                data,
                title=title,
                colors=colors,
                parental_labels=labels,
                pop_sizes=pop_sizes,
                pop_labels=pop_labels,
            )
        else:
            ax = admixture_bars(
                data,
                title=title,
                colors=colors,
                ylabel=ylabel or 'frequencies',
                parental_labels=labels,
                pop_sizes=pop_sizes,
                pop_labels=pop_labels,
            )
        if show:
            plt.show()
        return ax

    #
    # Scatter plots as projections
    #
    def scatter(self, which='pca', projection_kwargs=None, **kwargs):
        """
        A scatter plot of genetic data embeded into a 2D manifold.

        Args:
            which (str):
                The projection method used to reduce dimensionality. Can be
                anyone of 'pca', 'tsne', 'mds', or 'isomap'.
            data:
                Method used to constructinitial raw data for visualization. This
                is the same argument as in :cls:`kpop.Population.as_array`
            title:
                Optional plot title.
            colors:
                A colormap or a list of colors for each allele.
            merge (bool):
                If True, treats a MultiPopulation as a single population.
                Otherwise separate data for each sub-population in the graph.
        """

        scatter = {
            'merge', 'colors', 'title', 'legend', 'axes', 'alpha', 'fuzzy',
            'labels',
        }
        coords_kwargs = {k: v for k, v in kwargs.items() if k not in scatter}
        scatter_kwargs = {k: v for k, v in kwargs.items() if k in scatter}
        if projection_kwargs:
            coords_kwargs.update(projection_kwargs)

        kwargs.pop('self', None)
        pop = self._population
        coords = pop.projection(2, which, **coords_kwargs)
        return self.scatter_coords(coords, **scatter_kwargs)

    def scatter_coords(self, coords, merge=False, colors=None, title=None,
                       legend=True, axes=None, alpha=1.0, fuzzy=False,
                       labels=None):
        """
        Plot a 2D projection of population data from an array of 2D coordinates.
        This method is used internally by all scatter plot methods. It can also
        be useful to test projection methods not supported by kpop.

        Args:
            coords (2d matrix):
                The dataset in the reduced set of coordinates.
            title:
                Optional plot title.
            colors:
                A colormap or a list of colors for each allele.
            merge (bool):
                If True, treats a MultiPopulation as a single population.
                Otherwise separate data for each sub-population in the graph.
            show (bool):
                If True (default), immediately displays graph on screen.

        Returns:
             A matplotlib axes.
        """

        if fuzzy:
            coords = fuzzyfy_coords(coords, fuzzy)

        if merge:
            coords_list = [coords]
            pop_ids = [self._population.id or None]
        elif labels is not None:
            pop_ids = sorted(set(labels))
            coords_list = [np.array(
                [pt for id_, pt in zip(labels, coords) if id_ == id]
            ) for id in pop_ids]
        else:
            pop_sizes = [len(pop) for pop in self._populations]
            coords_list = group_individuals(coords, pop_sizes)
            pop_ids = [pop.id or i for i, pop in
                       enumerate(self._populations)]

        colors = _colors(colors, len(coords_list))

        # Plot scattered elements
        ax = axes or plt.axes()
        for i, coords in enumerate(coords_list):
            X, Y = np.asarray(coords).T
            label = pop_ids[i]
            ax.plot(X, Y, 'o', color=colors[i], label=label, alpha=alpha)

        # Additional plot elements
        if title:
            plt.title(title, axes=ax)
        if legend:
            ax.legend()

        return ax

    def pca(self, title='Principal component analysis', **kwargs):
        """
        A 2D principal component analysis plot for population.

        Args:
            data:
                Parameters of PCA. See :meth:`PopulationBase.pca`.
            title:
                Optional plot title.
            colors:
                A colormap or a list of colors for each allele.
            merge (bool):
                If True, treats a MultiPopulation as a single population.
                Otherwise separate data for each sub-population in the graph.
            show (bool):
                If True, immediately displays graph on screen.

        Returns:
             A matplotlib axes.

        See Also:
            :method:`kpop.population.projection.pca`
        """

        return self.scatter('pca', title=title, **kwargs)

    def kernel_pca(self, title='Kernel PCA', **kwargs):
        """
        Kernel PCA.

        See Also:
            :method:`kpop.population.projection.projection.kernel_pca`
        """

        return self.scatter('kernel_pca', title=title, **kwargs)

    def nmf(self, title='Non-negative matrix factorization', regularization=0,
            **kwargs):
        """
        Non-negative matrix factorization.

        See Also:
            :method:`kpop.population.projection.projection.nmf`
        """

        return self.scatter(
            'nmf',
            title=title,
            projection_kwargs={'alpha': regularization},
            **kwargs
        )

    def ica(self, title='Independent Components Analysis', **kwargs):
        """
        ICA

        See Also:
            :method:`kpop.population.projection.projection.nmf`
        """

        return self.scatter('ica', title=title, **kwargs)

    def tsne(self, title='t-SNE', **kwargs):
        """
        t-Distributed Stochastic Neighbor Embedding. This is a widely used
        method to project high dimensional data for visualization.

        See Also:
            :method:`kpop.population.projection.projection.tsne`
        """

        return self.scatter('tsne', title=title, **kwargs)

    def mds(self, title='Multidimensional scaling', **kwargs):
        """
        Multidimensional scaling.

        See Also:
            :method:`kpop.population.projection.projection.mds`
        """

        return self.scatter('mds', title=title, **kwargs)

    def isomap(self, title='Isomap', **kwargs):
        """
        Isomap

        See Also:
            :method:`kpop.population.projection.projection.isomap`
        """

        return self.scatter('isomap', title=title, **kwargs)

    def lle(self, title='Locally Linear Embedding', **kwargs):
        """
        Local Linear Embedding.

        See Also:
            :method:`kpop.population.projection.projection.lle`
        """

        return self.scatter('lle', title=title, **kwargs)

    def spectral(self, title='Spectral Embedding', **kwargs):
        """
        Spectral Embedding.

        See Also:
            :method:`kpop.population.projection.projection.spectral`
        """

        return self.scatter('spectral', title=title, **kwargs)

    def factor(self, title='Factor Analysis', **kwargs):
        return self.scatter('factor', title=title, **kwargs)

    def scatter_lda(self, n_populations=5,
                    title='Latent Dirichilet Allocation', **kwargs):
        """
        A scatter plot of LDA data. It computes the Q-matrix of admixture
        coefficients and project it in a plane using PCA. This is probably not
        as useful as an admixture plot but can give some insight on the
        structure of a population.

        Args:
            n_populations:
                Number of parental populations.
        """
        kwargs.setdefault('n_populations', n_populations)
        return self.scatter('lda', title=title, **kwargs)

    #
    # Admixture plots
    #
    def admixture(self, parental_labels=None, scatter=False, **kwargs):
        if self._population.parent is None:
            raise ValueError('population must define a parental population')

        data = np.array([ind.admixture_vector for ind in self._population])
        kwargs['pop_sizes'] = self._pop_sizes()
        kwargs['pop_labels'] = self._pop_labels()
        kwargs['parental_labels'] = self._pop_labels(self._population.parent)
        if scatter:
            return admixture_scatter(data, **kwargs)
        else:
            return admixture_bars(data, **kwargs)


def fuzzyfy_coords(coords, scale):
    """
    Fuzzyfy coords of list of points. Useful to display data points occluded
    in a scatter plot.
    """

    scale = float(scale)
    n = len(coords)
    xmin = coords[:, 0].min()
    xmax = coords[:, 0].max()
    ymin = coords[:, 1].min()
    ymax = coords[:, 1].max()
    dx = scale * (xmax - xmin) / 50
    dy = scale * (ymax - ymin) / 50
    random_dx = np.random.uniform(-dx, dx, size=n)
    random_dy = np.random.uniform(-dy, dy, size=n)
    coords = np.array(coords)
    coords[:, 0] += random_dx
    coords[:, 1] += random_dy
    return coords

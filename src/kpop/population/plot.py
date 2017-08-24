import numpy as np
from matplotlib import pyplot as plt

from kpop.plot import admixture_scatter, admixture_bars
from kpop.plot.utils import group_individuals, _colors


class PlotAttribute:
    """
    Implements the Population.plots attribute.
    """

    @property
    def _populations(self):
        return self._population.populations

    def __init__(self, population):
        self._population = population

    def _pop_sizes(self):
        return [len(pop) for pop in self._populations]

    def _pop_labels(self, population=None):
        populations = (population or self._population).populations
        labels = []
        for i, pop in enumerate(populations):
            label = pop.label if pop.label else 'pop-%s' % i
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
        pop_labels = [pop.label or i for i, pop in enumerate(populations)]
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
    def scatter(self, which='pca', **kwargs):
        """
        A scatter plot of genetic data embeded into a 2D manifold.

        Args:
            which (str):
                The projection method used to reduce dimensionality. Can be
                anyone of 'pca', 'tsne', 'mds', or 'isomap'.
            method, norm:
                Method used to represent data before the dimensionality
                reduction.
            title:
                Optional plot title.
            colors:
                A colormap or a list of colors for each allele.
            merge (bool):
                If True, treats a MultiPopulation as a single population.
                Otherwise separate data for each sub-population in the graph.
        """

        scatter = {'merge', 'colors', 'show', 'title', 'legend', 'axes'}
        coords_kwargs = {k: v for k, v in kwargs.items() if k not in scatter}
        scatter_kwargs = {k: v for k, v in kwargs.items() if k in scatter}

        kwargs.pop('self', None)
        pop = self._population
        coords = pop.projection(which, 2, **coords_kwargs)
        return self.scatter_coords(coords, **scatter_kwargs)

    def scatter_coords(self, coords, merge=False, colors=None, show=None,
                       title=None, legend=True, axes=None):
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

        if merge:
            coords_list = [coords]
            pop_labels = [self._population.label or None]
        else:
            pop_sizes = [len(pop) for pop in self._populations]
            coords_list = group_individuals(coords, pop_sizes)
            pop_labels = [pop.label or i for i, pop in
                          enumerate(self._populations)]

        colors = _colors(colors, len(coords_list))

        # Plot scattered elements
        ax = axes or plt.axes()
        for i, coords in enumerate(coords_list):
            X, Y = np.asarray(coords).T
            label = pop_labels[i]
            ax.plot(X, Y, 'o', color=colors[i], label=label)

        # Additional plot elements
        if title:
            plt.title(title, axes=ax)
        if legend:
            ax.legend()

        return ax

    def pca(self, method='count', merge=False, colors=None,
            show=None, title='Principal component analysis', legend=True,
            **kwargs):
        """
        A 2D principal component analysis plot for population.

        Args:
            method, norm:
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
            :method:`kpop.population.projection.projection.pca`
        """

        return self.scatter('pca', **fix_locals(locals()))

    def tsne(self, method='count', merge=False, colors=None,
            show=None, title='t-SNE', legend=True, **kwargs):
        """
        t-Distributed Stochastic Neighbor Embedding. This is a widely used
        method to project high dimensional data for visualization.

        See Also:
            :method:`kpop.population.projection.projection.tsne`
        """

        return self.scatter('tsne', **fix_locals(locals()))

    def mds(self, method='count', merge=False, colors=None,
            show=None, title='Multidimensional scaling', legend=True, **kwargs):
        """
        Multidimensional scaling.

        See Also:
            :method:`kpop.population.projection.projection.mds`
        """

        return self.scatter('mds', **fix_locals(locals()))

    def isomap(self, method='count', merge=False, colors=None,
            show=None, title='Isomap', legend=True, **kwargs):
        """
        Isomap

        See Also:
            :method:`kpop.population.projection.projection.isomap`
        """

        return self.scatter('isomap', **fix_locals(locals()))

    def lle(self, method='count', merge=False, colors=None,
            show=None, title='LLE', legend=True, **kwargs):
        """
        Local Linear Embedding.

        See Also:
            :method:`kpop.population.projection.projection.lle`
        """

        return self.scatter('lle', **fix_locals(locals()))

    def spectral(self, method='count', merge=False, colors=None,
            show=None, title='Spectral embedding', legend=True, **kwargs):
        """
        Spectral Embedding.

        See Also:
            :method:`kpop.population.projection.projection.spectral`
        """

        return self.scatter('spectral', **fix_locals(locals()))




    #
    # Admixture plots
    #
    def admixture(self, parental_labels=None,
                  scatter=False, show=True, **kwargs):
        if self._population.parent is None:
            raise ValueError('population must define a parental population')

        data = np.array([ind.admixture_vector for ind in self._population])
        kwargs['pop_sizes'] = self._pop_sizes()
        kwargs['pop_labels'] = self._pop_labels()
        kwargs['parental_labels'] = self._pop_labels(self._population.parent)
        if scatter:
            ax = admixture_scatter(data, **kwargs)
        else:
            ax = admixture_bars(data, **kwargs)
        if show:
            plt.show()
        return ax


def fix_locals(ns):
    "Pops the 'self' key from dictionary and return dictionary."

    ns.pop('self', None)
    ns.update(ns.pop('kwargs', {}))
    return ns
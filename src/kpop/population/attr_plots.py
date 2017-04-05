import numpy as np
from matplotlib import pyplot as plt

from kpop import plots as kpop_plots
from kpop.plots import admixture_scatter, admixture_bars, group_individuals


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

    def pca(self, method='flatten', norm=False, merge=False, colors=None,
            show=True, title='Principal component analysis', legend=True):
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
                If True (default), immediately displays graph on screen.

        Returns:
             A matplotlib axes.
        """

        pca_coords = self._population.pca(k=2, method=method, norm=norm)
        return self.projection(pca_coords, merge=merge, colors=colors,
                               show=show, title=title, legend=legend)

    def projection(self, coords, merge=False, colors=None, show=True,
                   title=None, legend=True):
        """
        Plot a 2D projection of population data.

        Args:
            coords (2d matrix):
                The dataset in the reduced set of coordinates.
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

        colors = kpop_plots._colors(colors, len(coords_list))

        # Plot scattered elements
        ax = plt.axes()
        for i, coords in enumerate(coords_list):
            X, Y = np.asarray(coords).T
            label = pop_labels[i]
            ax.plot(X, Y, 'o', color=colors[i], label=label)

        # Additional plot elements
        if title:
            plt.title(title, axes=ax)
        if legend:
            ax.legend()

        # Show and return element
        if show:
            plt.show()
        return ax

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

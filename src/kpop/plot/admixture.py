import numpy as np
from matplotlib import pyplot as plt

from .utils import normalize_coeffs, _colors, _pop_labels, _pop_sizes, \
    sort_coeffs, group_individuals, admixture_coords


def admixture_bars(coeffs, colors=None, legend=True,
                   ylabel='Ancestry proportions',
                   pop_sizes=None, pop_labels=None, parental_labels=None,
                   separaton_lines=True, sort=None,
                   title='Admixture coefficients', axes=None):
    """
    Makes an bar plot of admixture coefficients.

    It accepts all arguments as :func:`admixture_scatter`.

    Args:
        coeffs:
            An array with admixture coefficents. Each row is an individual and
            each column represents the admixture coefficients for parental
            populations. These coefficients do not need to be normalized.
        colors (str or Sequence):
            A colormap or a sequence of colors for each population.
        legend (bool):
            If True (default), display a legend with sub-populations.
        ylabel:
            Label of the y-axis.
        pop_sizes (list[int]):
            A list of population sizes. The sum of all pop_sizes must be equal
            to the number of individuals in the coeffs data.
        pop_labels (list[str]):
            A list of population labels. Must be of the same size as pop_sizes.
        parental_labels (list[str]):
            A list of labels for the parental populations. The bar plot shows
            the proportions for each parental population.
        separaton_lines (bool):
            If true, draws a thin line separating each population.
        sort (int):
            Corresponds to the index of a parental population that is used to
            sort individuals. If given, individuals within each population are
            sorted by ancestry coefficient for the given parental population.
        title (str):
            Plot's title.


    Returns:
        A matplotlib axes instance.
    """
    coeffs = normalize_coeffs(coeffs)
    num_individuals, num_parental = coeffs.shape
    colors = _colors(colors, num_parental)
    parental_labels = _pop_labels(parental_labels, num_parental, 'parental-%s')
    pop_sizes = _pop_sizes(pop_sizes, num_individuals)

    # Sort by coefficient
    if sort is True:
        sort = 0
    if sort not in (None, False):
        coeffs = sort_coeffs(coeffs, pop_sizes, sort)

    # Adjust pop sizes and names
    if pop_sizes is not None:
        pop_ticks = np.add.accumulate(pop_sizes) - np.asarray(pop_sizes) / 2
        pop_labels = _pop_labels(pop_labels, len(pop_sizes))
    else:
        pop_ticks = None

    # Create stacked bars
    ax = plt.axes() if axes is None else axes
    X = np.arange(num_individuals)
    bottom = np.zeros(num_individuals, dtype=float)
    for i, Y in enumerate(coeffs.T):
        color = colors[i]
        label = parental_labels[i]
        ax.bar(X, Y, bottom=bottom, color=color, width=1.0, linewidth=0,
               label=label)
        bottom += Y
    ax.axis([0, num_individuals, 0, 1])

    # Create vertical lines for separating populations
    if separaton_lines and pop_sizes is not None:
        for x in np.add.accumulate(pop_sizes)[:-1]:
            ax.plot([x, x], [0, 1], 'k--', lw=2)

    # Create additional graphic elements
    if legend:
        ax.legend()
    if ylabel:
        plt.ylabel(ylabel, axes=ax)
    if pop_ticks is not None:
        plt.xticks(pop_ticks, pop_labels, axes=ax)
    else:
        plt.xticks([], axes=ax)
    if title:
        plt.title(title, axes=ax)
    return ax


def admixture_scatter(coeffs, colors=None, legend=True,
                      pop_sizes=None, pop_labels=None, parental_labels=None,
                      title='Admixture coefficients'):
    """
    Makes an scatter plot of admixture coefficients using polygonal projections.
    For k=3, this is the traditional triangle plot and it completly encodes a
    probability distribution into a point inside the triangle.

    For more than 3 coordinates, the polygonal projection ceases to be bijective
    and each point in the graph is associated with many different probability
    distributions.

    Args:
        coeffs:
            An array with admixture coefficents. Each row is an individual and
            each column represents the admixture coefficients for parental
            populations. These coefficients do not need to be normalized.
        colors (str or Sequence):
            A colormap or a sequence of colors for each population.
        legend (bool):
            If True (default), admix legend.
        pop_sizes (list[int]):
            A list of population sizes. The sum of all pop_sizes must be equal
            to the number of individuals in the coeffs data.
        pop_labels (list[str]):
            A list of population labels. Must be of the same size as pop_sizes.
        parental_labels (list[str]):
            A list of labels for the parental populations.
        title (str):
            Plot's title.

    Examples:
        Consider we have 3 `Mocambo` individuals and 4 `Kalungas`. The
        admixture coefficients for 'african', 'european' and 'american' parental
        populations are:

        >>> mocambo = [
        ...     [0.8, 0.1, 0.1],
        ...     [0.7, 0.2, 0.1],
        ...     [0.8, 0.2, 0.0],
        ... ]
        >>> kalunga = [
        ...     [0.5, 0.2, 0.3],
        ...     [0.6, 0.2, 0.2],
        ...     [0.7, 0.1, 0.2],
        ...     [0.6, 0.1, 0.3],
        ... ]
        >>> coeffs = mocambo + kalunga
        >>> ax = admixture_scatter(
        ...     coeffs,
        ...     pop_sizes=[3, 4],
        ...     pop_labels=['mocambo', 'kalunga'],
        ...     parental_labels=['african', 'european', 'american'])
        >>> ax.show()

    Returns:
        A matplotlib axes instance.
    """
    data = normalize_coeffs(coeffs)
    num_individuals, num_parental = coeffs.shape
    pop_sizes = _pop_sizes(pop_sizes, num_individuals)
    populations = group_individuals(data, pop_sizes)
    num_pops = len(populations)
    parental_labels = _pop_labels(parental_labels, num_parental, 'parental-%s')
    colors = _colors(colors, num_pops)
    pop_labels = _pop_labels(pop_labels, num_pops)

    # Create vertices
    vertices = [[0, 1]]
    theta = 2 * np.pi / num_parental
    while len(vertices) < num_parental:
        x, y = vertices[-1]
        vertices.append([x * np.cos(theta) - y * np.sin(theta),
                         x * np.sin(theta) + y * np.cos(theta)])
    vertices = np.array(vertices)

    # Draw names of parental populations
    ax = plt.axes()
    ax.axis('scaled')
    if num_parental == 3:
        ax.axis([-1.1, 1.1, -0.7, 1.2])
        delta = (0, 0.1)
        top, left, right = vertices
        kwargs = {'verticalalignment': 'center',
                  'horizontalalignment': 'center'}
        ax.text(*(top + delta), s=parental_labels[0], **kwargs)
        ax.text(*(left - delta), s=parental_labels[1], **kwargs)
        ax.text(*(right - delta), s=parental_labels[2], **kwargs)
    else:
        ax.axis([-1.2, 1.2, -1.2, 1.2])
        for label, vertex in zip(parental_labels, vertices):
            x, y = vertex * 1.1
            ax.text(x, y, label,
                    verticalalignment='center',
                    horizontalalignment='center',
                    rotation=180 * np.arctan2(-x, y) / np.pi)

    # Create circuit path
    X, Y = np.array(list(vertices) + [vertices[0]]).T
    ax.plot(X, Y, '-', lw=2, color='k')

    # Create plots for each sub-population
    for i, pop in enumerate(populations):
        label = pop_labels[i]
        color = colors[i]
        X, Y = admixture_coords(pop, vertices).T
        ax.plot(X, Y, 'o', color=color, label=label)

    # Additional graphic elements
    if legend:
        ax.legend()
    if title:
        plt.title(title, axes=ax)
    plt.xticks([], axes=ax)
    plt.yticks([], axes=ax)
    return ax

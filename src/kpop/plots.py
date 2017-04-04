import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

__all__ = ['admixture_bars', 'admixture_scatter']


def unique_colors(n, colormap=cm.rainbow):
    """
    Iterator with n unique colors from colormap.
    """

    if isinstance(colormap, str):
        try:
            colormap = getattr(cm, colormap)
        except AttributeError:
            raise ValueError('invalid colormap: {0!r}'.format(colormap))
    size = colormap.N
    yield colormap(0)
    if n > 1:
        delta = size / (n - 1)
        x = 0
        for i in range(n - 1):
            x += delta
            yield colormap(int(x))


def normalize_coeffs(coeffs):
    """
    Normalize admixture coefficients.
    """

    coeffs = np.asarray(coeffs, dtype=float)
    norm = coeffs.sum(axis=1)
    return coeffs / norm[:, None]


def _colors(colors, n):
    if colors is None:
        colors = unique_colors(n)
    elif isinstance(colors, str):
        colors = unique_colors(n, colormap=colors)
    colors = list(colors)
    return colors


def _pop_sizes(pop_sizes, total):
    if pop_sizes is not None:
        pop_sizes = np.array(pop_sizes, dtype=int)
        if sum(pop_sizes) != total:
            fmt = total, sum(pop_sizes)
            raise ValueError('pop_sizes assumes {0!s} individuals, got {1!s}'.format(*fmt))
    return pop_sizes


def _pop_labels(pop_labels, n, fmt='pop-%s'):
    if pop_labels is None:
        pop_labels = [fmt % n for n in range(1, n + 1)]
    elif len(pop_labels) != n:
        raise ValueError('expect {0!s} populations, got {1!r}'.format(n, list(pop_labels)))
    return pop_labels


def group_individuals(data, chunk_sizes=None):
    """
    Group elements in sequence data in chunks of given sizes.

    Args:
        data: list of elements
        chunk_sizes: list of chunk sizes.
    """
    if chunk_sizes is None:
        chunk_sizes = [len(data)]

    if sum(chunk_sizes) != len(data):
        raise ValueError('expect {0!s} individuals, got {1!s}'.format(sum(chunk_sizes), len(data)))

    result = []
    data = iter(data)
    for size in chunk_sizes:
        chunk = [x for i, x in zip(range(size), data)]
        result.append(chunk)
    return result


def sort_coeffs(coeffs, pop_sizes, index):
    """
    Sort individuals in sub-populations by given ancestry index.
    """

    sub_pops = group_individuals(coeffs, pop_sizes)
    for group in sub_pops:
        group.sort(key=lambda x: x[index])
    result = np.concatenate(sub_pops)
    return result


def admixture_bars(coeffs, colors=None, legend=True,
                   ylabel='Ancestry proportions',
                   pop_sizes=None, pop_labels=None, parental_labels=None,
                   separaton_lines=True, sort=0,
                   title='Admixture coefficients', axes=None):
    """
    Makes an bar plot of admixture coefficients.

    It accepts all arguments as :func:`admixture_scatter`.

    Args:
        ylabel:
            Label of the y-axis.

    Returns:
        A matplotlib axes instance.
    """
    coeffs = normalize_coeffs(coeffs)
    num_individuals, num_parental = coeffs.shape
    colors = _colors(colors, num_parental)
    parental_labels = _pop_labels(parental_labels, num_parental, 'parental-%s')
    pop_sizes = _pop_sizes(pop_sizes, num_individuals)

    # Sort by coefficient
    if sort is not None or sort is not False:
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
    if separaton_lines:
        if pop_sizes is not None:
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


def _admix_coords(pop, vertices):
    result = []
    for proportions in pop:
        coords = np.array([0, 0], dtype=float)
        for p, v in zip(proportions, vertices):
            coords += p * v
        result.append(list(coords))
    return np.array(result)


def admixture_scatter(coeffs, colors=None, legend=True,
                      pop_sizes=None, pop_labels=None, parental_labels=None,
                      title='Admixture coefficients'):
    """
    Makes an scatter plot of admixture coefficients.

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
            A list of population sizes.
        pop_labels (list[str]):
            A list of population labels.
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
        X, Y = _admix_coords(pop, vertices).T
        ax.plot(X, Y, 'o', color=color, label=label)

    # Additional graphic elements
    if legend:
        ax.legend()
    if title:
        plt.title(title, axes=ax)
    plt.xticks([], axes=ax)
    plt.yticks([], axes=ax)
    return ax


if __name__ == '__main__':
    coeffs = np.array([
        [20, 2, 0],
        [2, 2, 3],
        [1, 0, 4],
        [5, 2, 1],
        [5, 3, 2],
        [5, 2, 2],
        [5, 2, 3],
    ])

    admixture_scatter(coeffs, pop_sizes=[3, 4], pop_labels=['foo', 'bar'],
                      parental_labels=['africa', 'europe', 'america'])
    plt.show()

    admixture_bars(coeffs, pop_sizes=[3, 4], pop_labels=['foo', 'bar'],
                   parental_labels=['africa', 'europe', 'america'])
    plt.show()

from ..libs import np

from kpop.libs import lazy_module

cm = lazy_module('matplotlib.cm')


def unique_colors(n, colormap=None):
    """
    Iterator with n unique colors from colormap.
    """
    colormap = colormap or cm.rainbow

    if isinstance(colormap, str):
        try:
            colormap = getattr(cm, colormap)
        except AttributeError:
            raise ValueError('invalid colormap: %r' % colormap)
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
            raise ValueError('pop_sizes assumes %s individuals, got %s' % fmt)
    return pop_sizes


def _pop_labels(pop_labels, n, fmt='pop-%s'):
    if pop_labels is None:
        pop_labels = [fmt % n for n in range(1, n + 1)]
    elif len(pop_labels) != n:
        raise ValueError('expect %s populations, got %r' %
                         (n, list(pop_labels)))
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
        raise ValueError('expect %s individuals, got %s' %
                         (sum(chunk_sizes), len(data)))

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


def admixture_coords(pop, vertices):
    result = []
    for proportions in pop:
        coords = np.array([0, 0], dtype=float)
        for p, v in zip(proportions, vertices):
            coords += p * v
        result.append(list(coords))
    return np.array(result)

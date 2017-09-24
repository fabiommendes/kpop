import functools

from .libs import np, sk_preprocessing


def normalized(func):
    """
    Normalize the results of func.
    """

    @functools.wraps(func)
    def decorated(self, *args, **kwargs):
        data = func(self, *args, **kwargs)
        return self._scale(data)

    if 'norm' not in func.__name__:
        name = func.__name__ + '_norm'
        decorated.__name__ = decorated.__qualname__ = name

    return decorated


class DataConverter:
    """
    Convert genetic data shaped as (size, num_loci, ploidy) into other useful
    formats.
    """

    def __init__(self, data, dtype=None):
        self.data = np.asarray(data, dtype=dtype or np.uint8)
        if self.data.ndim != 3:
            raise ValueError('expect data shaped as [size, num_loci, ploidy]')

    def __call__(self, method, **kwargs):
        def error():
            raise ValueError('invalid method: %r' % method)

        method_ = method.lower().replace('-', '_').replace(' ', '_')
        if method.startswith('_') or method == 'data':
            error()
        try:
            func = getattr(self, method_)
        except AttributeError:
            error()
        else:
            return func(**kwargs)

    def _scale(self, data):
        return sk_preprocessing.scale(data.astype(float))

    def randomized(self):
        data = self.data.copy()
        np.random.shuffle(data)
        return data

    def raw(self):
        return self.data

    def raw_norm(self):
        data = self.data
        data = data - data.mean(axis=0)
        std = data.std(axis=0)
        data /= np.where(std, std, 1)
        return data

    def flat(self):
        size, num_loci, ploidy = self.data.shape
        return self.data.reshape(size, num_loci * ploidy)

    def rflat(self):
        size, num_loci, ploidy = self.data.shape
        return self.randomized().reshape(size, num_loci * ploidy)

    def count(self):
        return (np.array(self.data) == 1).sum(axis=2)

    def count_snp(self):
        count = self.count()
        mu = count.mean(axis=0)
        p = mu / self.data.shape[2]
        norm = np.sqrt(p * (1 - p))
        norm = np.where(norm, norm, 1)
        return (count - mu) / norm

    def count_center(self):
        count = self.count()
        ploidy = self.data.shape[2]
        if ploidy % 2:
            return count - ploidy / 2
        else:
            return count - ploidy // 2

    # Normalized versions
    rflat_norm = normalized(rflat)
    count_norm = normalized(count)
    flat_norm = normalized(flat)

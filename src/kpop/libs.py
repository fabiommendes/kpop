import importlib


class LazyModule:
    """
    A lazy module object.
    """

    def __init__(self, name, lock=False):
        self.__path = name
        self.__mod = None
        self.__lock = lock

    def __load(self):
        if self.__lock:
            raise RuntimeError('import is locked')
        self.__mod = importlib.import_module(self.__path)

    def __getattr__(self, item):
        if self.__mod is None:
            self.__load()
        value = getattr(self.__mod, item)
        setattr(self, item, value)
        return value


def lazy_module(mod, lock=False):
    """
    Load a lazy module.

    Import is done only after first attribute access.
    """
    return LazyModule(mod, lock)


# We set this to True to detect any module that loads an external
# library upon import. Lock should obviously always be False on
# production.
_lock = False

#
# Lazy modules. This is used to speed up kpop initialization.
# Loding those massive pydata modules is done lazyly.
#
np = lazy_module('numpy', _lock)
pd = lazy_module('pandas', _lock)
plt = lazy_module('matplotlib.pyplot', _lock)

# Kpop modules
kpop = lazy_module('kpop', _lock)
kpop_result = lazy_module('kpop.result', _lock)
kpop_admixture = lazy_module('kpop.admixture', _lock)

# Scipy sub-modules
sp = lazy_module('scipy', _lock)
sp_integrate = lazy_module('scipy.integrate', _lock)
sp_optimize = lazy_module('scipy.optimize', _lock)

# Scikit learn sub-modules
sklearn = lazy_module('sklearn', _lock)
sk_base = lazy_module('sklearn.base', _lock)
sk_decomposition = lazy_module('sklearn.decomposition', _lock)
sk_manifold = lazy_module('sklearn.manifold', _lock)
sk_naive_bayes = lazy_module('sklearn.naive_bayes', _lock)
sk_preprocessing = lazy_module('sklearn.preprocessing', _lock)
sk_svm = lazy_module('sklearn.svm', _lock)

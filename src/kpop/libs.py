from .utils import lazy_module

np = lazy_module('numpy')
pd = lazy_module('pandas')
mlt = lazy_module('matplotlib.plt')

# Scikit learn sub-modules
sklearn = lazy_module('sklearn')
sk_naive_bayes = lazy_module('sklearn.naive_bayes')

kpop = lazy_module('kpop')

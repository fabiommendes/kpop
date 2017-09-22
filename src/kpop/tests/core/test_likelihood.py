from kpop.admixture.likelihood import *
from kpop.utils import fill_freqs_vector
from math import log


def test_loglike():
    x = [[1, 1], [1, 2], [2, 2]]
    freqs = [0.2, 0.3, 0.4]
    logfreqs = np.log(fill_freqs_vector(freqs))

    result = (
        2 * log(0.2) + 0 * log(0.8) +
        1 * log(0.3) + 1 * log(0.7) +
        0 * log(0.4) + 2 * log(0.6)
    )
    assert abs(result - loglike(x, freqs)) < 1e-6
    assert abs(result - loglike_logf(x, logfreqs)) < 1e-6
    assert abs(like(x, freqs) - np.exp(loglike(x, freqs))) < 1e-6

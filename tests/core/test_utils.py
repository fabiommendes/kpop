from kpop.libs import lazy_module


def test_lazy_module():
    math = lazy_module('math')
    assert math.sqrt(4) == 2

    plt = lazy_module('matplotlib.pyplot')
    plt.plot

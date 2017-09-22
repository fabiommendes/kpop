import importlib

from lazyutils import lazy
from sidekick import fn

fn_property = lambda x: property(fn(x)._)  # noqa: E731
fn_lazy = lambda x: lazy(fn(x)._)  # noqa: E731


class LazyModule:
    """
    A lazy module object.
    """

    def __init__(self, name):
        self.__path = name
        self.__mod = None

    def __load(self):
        self.__mod = importlib.import_module(self.__path)

    def __getattr__(self, item):
        if self.__mod is None:
            self.__load()
        value = getattr(self.__mod, item)
        setattr(self, item, value)
        return value


def lazy_module(mod):
    """
    Load a lazy module.

    Import is done only after first attribute access.
    """
    return LazyModule(mod)

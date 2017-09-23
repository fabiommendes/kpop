from lazyutils import lazy
from sidekick import fn

fn_property = lambda x: property(fn(x)._)  # noqa: E731
fn_lazy = lambda x: lazy(fn(x)._)  # noqa: E731

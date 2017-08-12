from sidekick import fn
from lazyutils import lazy

fn_property = lambda x: property(fn(x)._)
fn_lazy = lambda x: lazy(fn(x)._)
 

from lazyutils import lazy

from .population_base import PopulationBase as _base
from .population_single import Population
from .population_multi import MultiPopulation
from .plot import plot
from .projection import projection

# Attributes
_base.plot = lazy(lambda x: plot(x))
_base.projection = lazy(lambda _: projection(_))

# Additional monkey patching
Population._multi_population_class = MultiPopulation

del lazy
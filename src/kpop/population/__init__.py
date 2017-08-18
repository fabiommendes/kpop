from lazyutils import lazy

from .population_base import PopulationBase as _base
from .population_single import Population
from .population_multi import MultiPopulation
from .projectionattribute import ProjectionAttribute
from .plot import PlotAttribute as _plot
from .projectionattribute import ProjectionAttribute as _projection

# Attributes
_base.plot = lazy(lambda x: _plot(x))
_base.projection = lazy(lambda _: _projection(_))

# Additional mokey patching
Population._compose_class = MultiPopulation

del lazy
from lazyutils import lazy

from .individual import Individual
from .population_base import PopulationBase as _base
from .population_single import Population
from .population_multi import MultiPopulation

# Additional monkey patching
Population._multi_population_class = MultiPopulation

del lazy
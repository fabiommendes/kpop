from .population_single import Population
from .population_multi import MultiPopulation
from .projection import Projection
Population._compose_class = MultiPopulation

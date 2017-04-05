from .population_single import Population
from .population_multi import MultiPopulation

Population._compose_class = MultiPopulation

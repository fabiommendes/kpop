from .population import Population
from .multi_population import MultiPopulation

Population._compose_class = MultiPopulation

import numpy as np

from kpop import Population
from kpop.population.io import format_from_path


class TestRender:
    "Test functions that render a population to string"

    def test_render_population(self, popA):
        assert popA.io.render() == '''
A1: 11 22 12 12 12
A2: 11 22 11 22 21
A3: 11 22 11 22 21
A4: 11 22 11 21 12
A5: 11 22 21 21 21
A6: 11 22 21 22 21
A7: 11 22 11 12 12
A8: 11 22 21 22 12
    '''.strip()

    def test_load_file(self, path):
        file = path('popAB.pickle')
        pop = Population.io.load(file)
        assert len(pop) == 12

    def test_file_format_inference(self):
        for fmt in ['csv', 'pickle', 'ped']:
            assert format_from_path('foo.' + fmt) == fmt

    def test_render_hides_info_from_large_populations(self):
        pop = Population(np.ones([10, 2, 2]))
        out = pop.io.render(max_ind=5)
        assert '...' in out

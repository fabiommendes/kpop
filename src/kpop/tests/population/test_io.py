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

    def test_render_population_biallelic_freqs(self, popA, popB):
        freqs_render = (popA + popB).stats.render_biallelic_freqs()
        assert freqs_render == '''
locus         A         B
   L1  1.000000  0.000000
   L2  0.000000  0.250000
   L3  0.750000  0.500000
   L4  0.250000  0.375000
   L5  0.500000  0.375000
'''[1:]

        assert (popA + popB).stats.render_biallelic_freqs(sep=', ') == '''
locus,        A,        B
   L1, 1.000000, 0.000000
   L2, 0.000000, 0.250000
   L3, 0.750000, 0.500000
   L4, 0.250000, 0.375000
   L5, 0.500000, 0.375000
'''[1:]

    def _test_render_population_frequencies(self, popA, popB):
        assert (popA + popB).render_biallelic_frequencies(sep=', ') == '''

'''.strip()
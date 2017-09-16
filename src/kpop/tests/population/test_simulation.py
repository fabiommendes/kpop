from kpop.simulation.genetic_drift import frequency_drift_1D, frequency_drift


class TestFrequencyDrift:
    def test_fixed_alleles_stay_fixed(self):
        data = [0, 1, 0.5]
        new = frequency_drift_1D(data, 2, 50)
        assert new[0] == 0
        assert new[1] == 1
        assert new[2] != 0.5

    def test_small_fluctuations_with_large_populations(self):
        data = [0.5, 0.5]
        new = frequency_drift(data, 1, 500)
        assert (0.4 <= new).all() and (new <= 0.6).all()

    def test_accepts_2d_data(self):
        data = [[0.5, 0.5], [0.5, 0.5]]
        new = frequency_drift(data, 1, 500)
        assert (0.4 <= new).all() and (new <= 0.6).all()

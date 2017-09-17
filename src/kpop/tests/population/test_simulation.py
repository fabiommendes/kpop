from kpop.simulation.genetic_drift import frequency_drift_1D, frequency_drift


class TestFrequencyDrift:
    "Test genetic drift functions"

    def test_fixed_alleles_stay_fixed(self):
        data = [0, 1, 0.5]
        new = frequency_drift_1D(data, 2, 50)
        assert new[0] == 0
        assert new[1] == 1
        assert abs(new[2] - 0.5) < 1e-1

    def test_small_fluctuations_with_large_populations(self):
        data = [0.5, 0.5]
        new = frequency_drift(data, 1, 500)
        assert (0.4 <= new).all() and (new <= 0.6).all()

    def test_accepts_2d_data(self):
        data = [[0.5, 0.5], [0.5, 0.5]]
        new = frequency_drift(data, 1, 500)
        assert (0.4 <= new).all() and (new <= 0.6).all()


class TestEvolutionAndOffspring:
    def test_create_individual_labels(self, popA):
        assert popA[0].id == 'A1'

    def test_create_offspring(self, popA):
        ind3 = popA.simulation.random_individual()
        ind2 = popA.simulation.new_offspring(id='Achild')
        ind1 = popA.simulation.new_individual(id='Arandom')

        for ind in ind1, ind2, ind3:
            print(ind)
            assert ind.population is popA
            assert ind.id is not None
            assert ind[0, 0] == 1
            assert ind[1, 0] == 2

    def test_population_genetic_drift(self, popA):
        size = 1e6  # we need a very large size to avoid going back to the same
        # frequence by pure chance

        # Frequencies do not change with evolution of zero generations
        pop2 = popA.simulation.genetic_drift(0, sample_size=0,
                                             population_size=size)
        assert popA.freqs == pop2.freqs

        # Genetic drift
        pop2 = popA.simulation.genetic_drift(10, sample_size=0,
                                             population_size=size)
        assert popA.freqs != pop2.freqs

        # Fixed alleles do not change
        assert popA.freqs[0] == pop2.freqs[0]
        assert popA.freqs[1] == pop2.freqs[1]

        # Non-fixed alleles should always change
        for i in range(2, 5):
            assert popA.freqs[i] != pop2.freqs[i]

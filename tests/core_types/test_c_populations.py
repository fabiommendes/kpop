import pytest
from sidekick import record

from kpop.population.population_list import PopulationList, \
    ImmutablePopulationList


class TestPopulationsList:
    @pytest.fixture
    def p1(self):
        return record(id='a')

    @pytest.fixture
    def p2(self):
        return record(id='b')

    def test_population_list_data_access(self, p1, p2):
        plist = PopulationList([p1, p2])
        assert plist['a'] == p1
        assert plist['b'] == p2
        assert repr(plist) == "[record(id='a'), record(id='b')]"

    def test_population_list_indexes(self, p1, p2):
        plist = PopulationList([p1, p2])
        assert plist.population_ids() == ['a', 'b']

    def test_population_list_mutation(self, p1, p2):
        plist = PopulationList([p1, p2])
        plist['a'] = p2
        assert plist == [p2, p2]

    def test_population_list_bad_index(self, p1, p2):
        plist = PopulationList([p1, p2])
        with pytest.raises(KeyError):
            plist['c']

    def test_immutable_population_list(self):
        p1 = record(id='a')
        p2 = record(id='b')
        plist = ImmutablePopulationList([p1, p2])

        with pytest.raises(TypeError):
            plist[0] = p1

        with pytest.raises(TypeError):
            del plist[0]
        with pytest.raises(TypeError):
            plist.insert(0, p1)

from kpop import Individual


def test_can_create_individual_from_data():
    ind = Individual([[1, 2], [3, 4]], label='foo')
    assert ind.render() == 'foo: 12 34'


def test_can_create_individual_from_string():
    ind = Individual('123 456', label='foo')
    assert ind.num_loci == 2
    assert ind.ploidy == 3
    assert ind.render() == 'foo: 123 456'


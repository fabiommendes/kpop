import numpy as np
import pytest
from numpy.testing import assert_equal

from kpop.freqmatrix import fill_frequencies
from kpop.individual import Individual, label_from_parents


@pytest.fixture()
def ind():
    return Individual('12 11 22 21')


@pytest.fixture()
def ind_data():
    return [[1, 2], [1, 1], [2, 2], [2, 1]]


def test_can_create_individual_from_data():
    ind = Individual([[1, 2], [3, 4]], label='foo')
    assert ind.render() == 'foo: 12 34'


def test_can_create_individual_from_string():
    ind = Individual('123 456', label='foo')
    assert ind.num_loci == 2
    assert ind.ploidy == 3
    assert ind.render() == 'foo: 123 456'


def test_individual_properties(ind, ind_data):
    assert ind.num_loci == 4
    assert ind.ploidy == 2
    assert ind.num_alleles == 2
    assert ind.is_biallelic
    assert not ind.has_missing
    assert ind.missing_ratio == 0
    assert_equal(ind.data, ind_data)
    assert ind.dtype == np.uint8


def test_render_individual(ind):
    assert repr(ind) == "Individual('ind: 12 11 22 21')"
    assert str(ind) == 'ind: 12 11 22 21'
    assert ind.render(limit=2) == 'ind: 12 ... 21'


def test_individual_as_sequence(ind, ind_data):
    assert_equal(list(ind), ind_data)
    assert_equal(ind[0], [1, 2])
    assert len(ind) == ind.num_loci


def test_breed(ind):
    child = ind.breed(ind)
    assert_equal(child[1], [1, 1])
    assert_equal(child[2], [2, 2])


def test_breed_monoic():
    ind = Individual('1 1 2 2', label='foo')
    child = ind.breed(ind)
    assert_equal(child.data, ind.data)
    assert child.label == 'foo'


def test_copy_individual(ind):
    cp = ind.copy()
    assert_equal(cp.data, ind.data)
    assert cp.label == ind.label
    assert cp.allele_names == ind.allele_names


def test_shuffle_loci(ind):
    ind_orig = ind.copy()
    sum_values = ind.data.sum(axis=-1)
    for _ in range(10):
        ind.shuffle_loci()
        assert (ind.data.sum(axis=-1) == sum_values).all()
        if (ind.data != ind_orig.data).any():
            break
    else:
        raise ValueError('individual not changed')


def test_haplotypes(ind):
    h1, h2 = ind.haplotypes()
    assert_equal(h1, [1, 1, 2, 2])
    assert_equal(h2, [2, 1, 2, 1])


def test_label_from_parents():
    f = label_from_parents
    assert f(None, None) is None
    assert f('foo', None) == 'foo_'
    assert f(None, 'foo') == 'foo_'
    assert f('foo', 'foo') == 'foo'
    assert f('fooA', 'fooB') == 'foo-A,B'
    assert f('foo', 'fooB') == 'foo-_,B'
    assert f('fooA', 'foo') == 'foo-A,_'


def test_random_individual():
    f1 = [1.0, 0.0, 0.5]
    ind1 = Individual.from_freqs(f1, label='rand1')
    ind2 = Individual.from_freqs(fill_frequencies(f1), label='rand2')

    for ind in ind1, ind2:
        print(ind)
        assert ind.ploidy == 2
        assert ind.num_loci == 3
        assert ind[0, 0] == 1
        assert ind[1, 0] == 2


def test_make_individual_from_data():
    ind = Individual('foo: ab aa bb')
    assert ind.label == 'foo'
    assert ind[0, 1] == 2

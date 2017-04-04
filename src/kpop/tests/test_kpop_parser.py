from kpop.kpop_parser import *
import pytest


def test_invalid_kpop_data():
    with pytest.raises(ValueError):
        str_to_data('01 22 11')

    with pytest.raises(IndexError):
        str_to_data('11 22 33', [{}, {}])


def test_missing_data():
    data, names = str_to_data('-1 22 11')
    assert data[0, 0] == 0
    assert data[0, 1] == 1


def test_alpha_data():
    data, names = str_to_data('ab 12 dc')
    assert data[0, 0] == 1
    assert data[2, 0] == 2


def test_has_same_names():
    data, names = str_to_data('ab aa bb', {})
    assert data[2, 0] == 2

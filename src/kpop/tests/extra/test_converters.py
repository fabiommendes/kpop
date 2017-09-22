import io

import pytest

from kpop import Population
from kpop.io import export_arlequin, csv_lines
from kpop.tests import load_data


@pytest.fixture
def popA():
    return Population('11 22 11 22\n12 12 12 12', id='popA')


@pytest.fixture
def popB():
    return Population('22 11 22 11\n12 12 12 12', id='popB')


@pytest.fixture
def popC():
    return Population('-- 11 22 11\n12 -- 12 12', id='popB')


class TestArlequin:
    "Test creation of Arlequin files"

    @pytest.fixture
    def arldata_AB(self):
        return load_data('arldata_AB.arl').read()

    def test_export_population_to_arlequin(self, popA, popB, arldata_AB):
        F = io.StringIO()
        export_arlequin(popA + popB, F)
        data = F.getvalue()
        assert data == arldata_AB


class TestCsvWriter:
    "Test creation of CSV files from population objects"

    def test_simple_csv_file(self):
        lines = csv_lines([[1, 2], [3, 4]], columns=['A', 'B'])
        lines = list(lines)
        a, b, c = lines
        assert a == "A,B\n"
        assert b == "1,2\n"
        assert c == "3,4\n"

    def test_align_to_the_right(self):
        lines = csv_lines([[1, 2], [10, 20]], columns=['A', 'B'], align=True)
        lines = list(lines)
        print(lines)
        a, b, c = lines
        assert a == " A, B\n"
        assert b == " 1, 2\n"
        assert c == "10,20\n"

    def test_align_to_the_left(self):
        lines = csv_lines([[1, 2], [10, 20]], columns=['A', 'B'], align='left')
        lines = list(lines)
        a, b, c = lines
        assert a == "A ,B \n"
        assert b == "1 ,2 \n"
        assert c == "10,20\n"

    def test_large_data(self):
        data = [('L1', 1.0, 0.0), ('L2', 0.0, 0.25), ('L3', 0.75, 0.5),
                ('L4', 0.25, 0.375), ('L5', 0.5, 0.375)]

        csv = ''.join(csv_lines(data, sep=' '))
        assert csv == (
            'L1 1.0 0.0\n'
            'L2 0.0 0.25\n'
            'L3 0.75 0.5\n'
            'L4 0.25 0.375\n'
            'L5 0.5 0.375\n'
        )

    def test_large_data_with_titles(self):
        data = [('L1', 1.0, 0.0), ('L2', 0.0, 0.25), ('L3', 0.75, 0.5),
                ('L4', 0.25, 0.375), ('L5', 0.5, 0.375)]

        csv = ''.join(csv_lines(data, sep=' ', columns=['a', 'b', 'c']))
        assert csv == (
            'a b c\n'
            'L1 1.0 0.0\n'
            'L2 0.0 0.25\n'
            'L3 0.75 0.5\n'
            'L4 0.25 0.375\n'
            'L5 0.5 0.375\n'
        )

import io

import pytest

from kpop import Population
from kpop.io import export_arlequin, csv_lines


@pytest.fixture
def popA():
    return Population('11 22 11 22\n12 12 12 12', label='popA')


@pytest.fixture
def popB():
    return Population('22 11 22 11\n12 12 12 12', label='popB')


@pytest.fixture
def popC():
    return Population('-- 11 22 11\n12 -- 12 12', label='popB')


arldata_AB = """[Profile]
Title="Sample Arlequin data"
NbSamples=2
DataType=STANDARD
GenotypicData=1
LocusSeparator=WHITESPACE
GameticPhase=0
RecessiveData=0
MissingData='-'


[Data]
[[Samples]]
SampleName="popA"
SampleSize=2
SampleData={
    popA1 1 1 2 1 2
            1 2 1 2
    popA2 1 1 1 1 1
            1 1 1 1
}

SampleName="popB"
SampleSize=2
SampleData={
    popB1 1 2 1 2 1
            2 1 2 1
    popB2 1 1 1 1 1
            1 1 1 1
}

"""


class TestArlequin:
    def test_export_population_to_arlequin(self, popA, popB):
        F = io.StringIO()
        export_arlequin(popA + popB, F)
        data = F.getvalue()
        assert data == arldata_AB


class TestCsvWriter:
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
        assert csv == '''
L1 1.0 0.0
L2 0.0 0.25
L3 0.75 0.5
L4 0.25 0.375
L5 0.5 0.375
'''[1:]

    def test_large_data_with_titles(self):
        data = [('L1', 1.0, 0.0), ('L2', 0.0, 0.25), ('L3', 0.75, 0.5),
                ('L4', 0.25, 0.375), ('L5', 0.5, 0.375)]

        csv = ''.join(csv_lines(data, sep=' ', columns=['a', 'b', 'c']))
        assert csv == '''
a b c
L1 1.0 0.0
L2 0.0 0.25
L3 0.75 0.5
L4 0.25 0.375
L5 0.5 0.375
'''[1:]

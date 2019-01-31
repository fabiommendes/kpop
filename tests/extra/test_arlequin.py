import io

import pytest

from kpop.io.arlequin import export_arlequin, create_sample_data


class TestArlequin:
    @pytest.fixture
    def pop(self, popA, popB):
        return popA + popB

    def test_convert_to_arlequin(self, pop):
        file = io.StringIO()
        export_arlequin(pop, file)
        assert file.getvalue() == POPAB_ARL

    def test_create_sample_data(self, popB):
        file = io.StringIO()
        create_sample_data(popB, file)
        assert file.getvalue() == POPB_DATA

    def test_create_single_sample(self, popB):
        file = io.StringIO()
        create_sample_data(popB[:1], file)
        assert file.getvalue() == '    B1 1 2 2 2 1 2\n         2 2 2 2 2\n'


#
# Examples
#
POPB_DATA = """
    B1 1 2 2 2 1 2
         2 2 2 2 2
    B2 1 2 1 1 2 2
         2 2 1 2 1
    B3 1 2 1 1 1 2
         2 2 1 2 1
    B4 1 2 2 2 2 1
         2 2 2 1 2
"""[1:]

POPAB_ARL = """
[Profile]
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
SampleName="A"
SampleSize=8
SampleData={
    A1 1 1 2 1 1 1
         1 2 2 2 2
    A2 1 1 2 1 2 2
         1 2 1 2 1
    A3 1 1 2 1 2 2
         1 2 1 2 1
    A4 1 1 2 1 2 1
         1 2 1 1 2
    A5 1 1 2 2 2 2
         1 2 1 1 1
    A6 1 1 2 2 2 2
         1 2 1 2 1
    A7 1 1 2 1 1 1
         1 2 1 2 2
    A8 1 1 2 2 2 1
         1 2 1 2 2
}

SampleName="B"
SampleSize=4
SampleData={
    B1 1 2 2 2 1 2
         2 2 2 2 2
    B2 1 2 1 1 2 2
         2 2 1 2 1
    B3 1 2 1 1 1 2
         2 2 1 2 1
    B4 1 2 2 2 2 1
         2 2 2 1 2
}

"""[1:]

import io
import pytest

from kpop import Population
from kpop.io import export_arlequin


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
LocusSeparator=SPACE
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

def test_export_population_to_arlequin(popA, popB):
    F = io.StringIO()
    export_arlequin(popA + popB, F)
    data = F.getvalue()
    assert data == arldata_AB
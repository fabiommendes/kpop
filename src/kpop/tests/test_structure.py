from pprint import pprint
import numpy as np
import pytest
from kpop.structure import run_structure


@pytest.mark.slow
def test_structure_can_detect_easy_parental_populations(popA, popB):
    res = run_structure(popA + popB, 2)
    ancestryA = [[0, 1] for x in range(10)]
    ancestryB = [[1, 0] for x in range(10)]
    ancestry = [list(x) for x in res.ancestry]
    assert ancestry == (ancestryA + ancestryB) or \
           ancestry == (ancestryB + ancestryA)
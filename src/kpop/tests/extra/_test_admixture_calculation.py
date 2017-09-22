import numpy as np
import pytest

from kpop import Population
from kpop.admixture.mixture_algorithm import AdmixtureCalculation


@pytest.fixture
def N():
    return 5


@pytest.fixture
def J():
    return 50


@pytest.fixture
def q_right(N):
    q_right = [[1, 0, 0]] * N
    q_right.extend([[0, 1, 0]] * N)
    q_right.extend([[0, 0, 1]] * N)
    return np.array(q_right)


@pytest.fixture
def pop(N, J):
    pA = Population.random(N, J, id='A', seed=1)
    pB = Population.random(N, J, id='B', seed=2)
    pC = Population.random(N, J, id='C', seed=3)
    return pA + pB + pC


@pytest.fixture
def alg(pop):
    return AdmixtureCalculation.from_pop(pop, 3)


def test_func_j_decomposition(alg):
    func = alg.func()
    func_ = sum(alg.func_j(j) for j in range(alg.num_loci))
    np.testing.assert_almost_equal(func, func_)


def test_func_i_decomposition(alg):
    func = alg.func()
    func_ = sum(alg.func_i(i) for i in range(alg.num_ind))
    np.testing.assert_almost_equal(func, func_)


def test_cython_implementations_matches_reference_qi(alg):
    for j in range(alg.num_loci):
        # Function fj
        cy = alg._func_j(j, alg.freqs[j])
        py = alg._func_j_python(j, alg.freqs[j])
        np.testing.assert_almost_equal(cy, py)

        # Jacobian fj
        cy = alg._jac_fj(j, alg.freqs[j])
        py = alg._jac_fj_python(j, alg.freqs[j])
        np.testing.assert_almost_equal(cy, py)

        # Hessian fj
        cy = alg._hess_fj(j, alg.freqs[j])
        py = alg._hess_fj_python(j, alg.freqs[j])
        np.testing.assert_almost_equal(cy, py)


def test_cython_implementations_matches_reference_fj(alg):
    for i in range(alg.num_ind):
        # Function qi
        cy = alg._func_i(i, alg.q[i])
        py = alg._func_i_python(i, alg.q[i])
        np.testing.assert_almost_equal(cy, py)

        # Jacobian qi
        cy = alg._jac_qi(i, alg.q[i])
        py = alg._jac_qi_python(i, alg.q[i])
        np.testing.assert_almost_equal(cy, py)

        # Hessian qi
        cy = alg._hess_qi(i, alg.q[i])
        py = alg._hess_qi_python(i, alg.q[i])
        np.testing.assert_almost_equal(cy, py)


def test_optimize_em(alg):
    alg2 = AdmixtureCalculation(alg.data, alg.k, alg.freqs, alg.q)

    for _ in range(5):
        alg._optimize_em()
        alg2._optimize_em_python()
        np.testing.assert_almost_equal(alg.q, alg2.q)
        np.testing.assert_almost_equal(alg.freqs, alg2.freqs)

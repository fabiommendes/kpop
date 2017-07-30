import os

import pytest

from kpop.external.structure import run_structure
from kpop.external.structure.parser import parse_file
from kpop.tests import is_ci

pytestmark = pytest.mark.skipif(is_ci,
                                reason='structure is not installed in travis')


@pytest.mark.slow
class TestStructureProgram:
    def test_structure_can_detect_easy_parental_populations(self, popA, popB):
        res = run_structure(popA + popB, 2, disp=2)
        ancestryA = [[0, 1] for x in range(10)]
        ancestryB = [[1, 0] for x in range(10)]
        # ancestry = [list(x) for x in res.ancestry]
        # assert ancestry == (ancestryA + ancestryB) or \
        #    ancestry == (ancestryB + ancestryA)


class TestStructureParser:
    @pytest.fixture(scope='class')
    def example(self):
        folder = os.path.dirname(__file__)
        with open(os.path.join(folder, 'structure.example')) as F:
            data = parse_file(F)
        return data

    def test_structure_parser_extracts_basic_info(self, example):
        assert example.k == 3
        assert example.individuals == 542
        assert example.loci == 46
        assert example.burnin == 10000
        assert example.RANDOMIZE == False
        assert 'SelP' in example.infile
        assert example.reps == 100000


    def test_structure_parser_extracts_extended_info(self, example):
        assert {'-e', '-m'}.issubset(example.argv)
        assert sorted(example.options.keys()) == [
            'ADMBURNIN', 'ALPHA', 'ALPHAMAX', 'ALPHAPRIORA', 'ALPHAPRIORB',
            'ALPHAPROPSD', 'ANCESTDIST', 'ANCESTPINT', 'BURNIN', 'COMPUTEPROB',
            'DATAFILE', 'EXTRACOLS', 'FPRIORMEAN', 'FPRIORSD', 'FREQSCORR',
            'GENSBACK', 'INFERALPHA', 'INFERLAMBDA', 'INTERMEDSAVE', 'LABEL',
            'LAMBDA', 'LINKAGE', 'LOCDATA', 'LOCISPOP', 'LOCPRIOR',
            'LOCPRIORINIT', 'LOCPRIORSTEP', 'MARKOVPHASE', 'MAXLOCPRIOR',
            'MAXPOPS', 'METROFREQ', 'MIGRPRIOR', 'MISSING', 'NOADMIX',
            'NUMBOXES', 'NUMINDS', 'NUMLOCI', 'NUMREPS', 'ONEFST', 'OUTFILE',
            'PHASED', 'PHASEINFO', 'PHENOTYPE', 'PLOIDY', 'POPALPHAS',
            'POPDATA', 'POPFLAG', 'POPSPECIFICLAMBDA', 'PRINTKLD',
            'PRINTLAMBDA', 'PRINTLIKES', 'PRINTNET', 'PRINTQHAT', 'PRINTQSUM',
            'RANDOMIZE', 'REPORTHITRATE', 'SEED', 'STARTATPOPINFO',
            'UNIFPRIORALPHA', 'UPDATEFREQ', 'USEPOPINFO',
        ]
        assert sorted(example.strat_options.keys()) == [
            'EMERROR', 'LOCUSxONLY', 'MISSINGPHENO', 'NUMSIMSTATS',
            'PHENOTYPECOL', 'POOLFREQ',
        ]

    def test_structure_parser_extracts_frequency_matrix(self, example):
        F = example.F_matrix
        assert list(F.columns) == \
            ['name', 'n_alleles', 'missing', 'data']

        locus1 = F.iloc[0]
        assert locus1['name'] == 'Locus 1'
        assert locus1['n_alleles'] == 2
        assert locus1['missing'] == 0

        data1 = locus1['data']
        assert list(data1.columns) == \
            ['total', 'cluster_0', 'cluster_1', 'cluster_2']
        assert (data1.as_matrix() == [
            [0.605,  0.813,  0.840,  0.434],
            [0.395,  0.187,  0.160,  0.566]
        ]).all()

    def test_structure_parser_extracts_fst(self, example):
        assert (example.mean_fst == [0.2843, 0.3614, 0.2333]).all()
        assert (example.avg_heterozygosity == [0.3061, 0.2843, 0.3244]).all()
        assert (example.allele_freq_divergence == [
            [0.,  0.1911,  0.119],
            [0.1911,  0.,  0.1857],
            [0.119,  0.1857,  0.]
        ]).all()

    def test_structure_parser_extracts_likelihood_data(self, example):
        assert example.mean_alpha == 0.0322
        assert example.data_log_like == -23816.3
        assert example.mean_log_like == -23695.0
        assert example.var_log_like == 242.6

    def test_structure_parser_extracts_mixture_data(self, example):
        Q = example.q_matrix
        assert list(Q.columns) == [
            'label', 'missing', 'population', 'cluster_0', 'cluster_1', 
            'cluster_2'
        ]
        ind0 = Q.iloc[0]
        assert ind0['label'] == 'HGDP00452'
        assert ind0['missing'] == 0
        assert ind0['population'] == 1
        assert ind0['cluster_0'] == .003
        assert ind0['cluster_1'] == .995
        assert ind0['cluster_2'] == .002

        # There are small rounding errors for mixture data. Maybe we should 
        # consider fixing them on the parser.
        totals = Q[['cluster_0', 'cluster_1', 'cluster_2']].sum(1)
        assert (abs(totals - 1) < 2e-3).all()

    def test_structure_parser_extracts_memberships(self, example):
        assert (example.pre_defined_membership.as_matrix() == [
            [0.012, 0.978, 0.010, 105],
            [0.011, 0.008, 0.980, 154],
            [0.944, 0.008, 0.048,  59],
            [0.981, 0.006, 0.012, 224],
        ]).all()
        
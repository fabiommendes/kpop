from ..libs import np
from sidekick import Record, field, record_to_dict

from ..io.utils import file_or_path

__all__ = ['parse', 'parse_out', 'parse_qfile', 'parse_pfile', 'Admixture']


class AdmixtureOut(Record):
    """
    A record that stores information from collected from an ADMIXTURE run.
    """

    method = field(str)
    acceleration = field(str)
    loglike = field(np.ndarray)
    loglike_em = field(np.ndarray)
    fst_matrix = field(np.ndarray)
    seed = field(int)
    duration = field(float)


class Admixture(AdmixtureOut):
    """
    Full result of an ADMIXTURE run.
    """

    pmatrix = field(np.ndarray)
    qmatrix = field(np.ndarray)


class AdmixtureOutParser:
    """
    ADMIXTURE output file parser.

    This is only an implementation detail. Please use the functional interface.
    """

    def __init__(self, file):
        self.file = file

    def parse(self) -> AdmixtureOut:
        """
        Execute parsing
        """
        raise NotImplementedError


@file_or_path
def parse_pfile(file):
    """
    Reads a .P file and return the corresponding numpy frequency matrix.
    """

    data = []
    for line in file:
        data.append(list(map(float, line.strip().split())))
    return np.array(data).T


@file_or_path
def parse_qfile(file):
    """
    Reads a .Q file and return the corresponding numpy admixture matrix.
    """

    data = []
    for line in file:
        data.append(list(map(float, line.strip().split())))
    return np.array(data)


@file_or_path
def parse_out(file) -> AdmixtureOut:
    """
    Parse a string the the output of a run of the ADMIXTURE program.

    It returns a record with the following fields:

    method ('em' or 'block'):
        Method used to compute results.
    acceleration ('sqs' or 'qn'):
        Acceleration method.
    loglike_em (array of floats):
        An array with the log likelihood for each run of the EM step.
    loglike (array of floats):
        An array with the log likelihood for each run of the main block
        relaxation step.
    fst_matrix (array):
        The resulting Fst divergence matrix for the estimated populations.
    seed (int):
        Random seed that initializes the algorithm.
    duration (float):
        Duration of the run.
    """

    parser = AdmixtureOutParser(file)
    return parser.parse()


def parse(outfile, pfile, qfile) -> Admixture:
    """
    Parse the result of an admixture run, including its outfile, .P file and
    .Q file and return the parsed result.
    """

    out = parse_out(outfile)
    pmatrix = parse_pfile(pfile)
    qmatrix = parse_qfile(qfile)
    kwargs = record_to_dict(out)
    kwargs.update(pmatrix=pmatrix, qmatrix=qmatrix)
    return Admixture(**kwargs)

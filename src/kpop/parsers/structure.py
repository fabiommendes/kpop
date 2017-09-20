import re
from collections import deque
from functools import lru_cache
from types import SimpleNamespace

import numpy as np
import pandas as pd


@lru_cache(maxsize=256)
def matcher(regex):
    "A cached version of re.match"

    return re.compile(regex).match


def match(re, st):
    "Match string with given regex."

    return matcher(re)(st)


class StructureParser:
    """
    Parse structure output files.
    """

    skip_whitespace_re = r'^\s*'

    def __init__(self, lines):
        self.lines = deque(lines)
        self.result = SimpleNamespace()

    def parse(self):
        self.parse_header()
        self.parse_config()
        self.parse_membership()
        self.parse_divergence()
        self.parse_avg_heterozygosity()
        self.parse_likelihood()
        self.parse_fst()
        self.parse_q()
        self.parse_f()
        self.parse_options()
        self.parse_strat()
        return self.result

    def pop(self):
        "Pops the first line."

        return self.lines.popleft()

    def push(self, value):
        "Append string to the begining of file"

        self.lines.appendleft(value)

    def parse_re(self, re):
        """
        Match regex with the first line and return the corresponding match
        object.
        """

        line = self.lines[0]
        m = match(re, line)
        if m is None:
            raise ValueError('unexpected line: %r' % line)

        new_line = self.lines[0][m.span()[1]:]
        if new_line:
            self.lines[0] = new_line
        else:
            self.lines.popleft()
        return m

    def parse_data(self, re, func=str):
        match = self.parse_re(re)
        return func(match.groups(0)[0])

    def skip_fragment(self, st):
        if not self.lines:
            raise ValueError('EOF')
        elif not self.lines[0].startswith(st):
            raise ValueError(
                'expect\n  %r, got\n  %r' % (st, self.lines[0].rstrip())
            )
        new_line = self.lines[0][len(st):]
        if new_line:
            self.lines[0] = new_line
        else:
            self.lines.popleft()

    def skip_blank(self):
        lines = self.lines

        while lines and (lines[0].isspace() or not lines[0]):
            lines.popleft()

    def skip_dashed_line(self):
        line = self.pop()
        if set(line) != {'-'}:
            raise ValueError('unexpected line: %r' % line)

    def skip_lines(self, n):
        for _ in range(n):
            self.lines.popleft()

    #
    # STRUCTURE FILE SECTIONS
    #
    def parse_header(self):
        self.skip_blank()
        self.pop()
        while not self.pop().startswith('-'):
            pass
        self.skip_blank()

    def parse_config(self):
        result = self.result

        # Cli args
        self.skip_fragment('Command line arguments:')
        result.argv = self.lines.popleft().strip().split()

        # Infile
        self.skip_fragment('Input File:')
        result.infile = self.lines.popleft().strip()
        self.skip_blank()

        # Run parameters
        self.skip_fragment('Run parameters:')
        result.individuals = self.parse_data(
            r'^\s*(\d+)\s*individuals$', int
        )
        result.loci = self.parse_data(
            r'^\s*(\d+)\s*loci$', int
        )
        result.k = self.parse_data(
            r'^\s*(\d+)\s*populations assumed$', int
        )
        result.burnin = self.parse_data(
            r'^\s*(\d+)\s*Burn-in period$', int
        )
        result.reps = self.parse_data(
            r'^\s*(\d+)\s*Reps$', int
        )

        # Randomize flag
        # TODO: accept other options
        # e.g.: STARTATPOPINFO, LOCPRIOR
        randomize = self.parse_data(
            r'^RANDOMIZE turned (\w+)$'
        )
        result.RANDOMIZE = {'off': False, 'on': True}[randomize]

        self.skip_blank()

    def parse_membership(self):
        # Skip table header
        self.skip_dashed_line()
        self.skip_lines(2)
        self.skip_blank()
        self.skip_lines(2)
        self.skip_blank()

        # Read lines until it reaches the dashed line
        lines = []
        index = []
        while True:
            line = self.pop()
            if line.startswith('-'):
                break
            else:
                cells = line.strip().split()
                line = [float(x) for x in cells[1:-1]]
                line.append(int(cells[-1]))
                lines.append(line)
                index.append(cells[0].rstrip(': '))

        # Convert lines to a Pandas dataframe
        self.result.pre_defined_membership = df = pd.DataFrame(index=index)
        for col in range(len(lines[0]) - 1):
            df['cluster_%s' % col] = [line[col] for line in lines]
        df['size'] = [line[-1] for line in lines]

        self.skip_blank()

    def parse_divergence(self):
        self.skip_lines(2)
        self.skip_blank()
        self.skip_lines(1)

        # Read lines
        lines = []
        for _ in range(self.result.k):
            lines.append(
                [float(x) if x != '-' else 0.0
                 for x in self.pop().strip().split()[1:]]
            )

        # Convert to matrix
        self.result.allele_freq_divergence = np.array(lines)
        self.skip_blank()

    def parse_avg_heterozygosity(self):
        pattern = r'^cluster\s*\d+\s*:\s*(\d+[.]\d+)\s*$'
        self.skip_lines(1)
        self.result.avg_heterozygosity = np.array(
            [self.parse_data(pattern, float) for _ in range(self.result.k)]
        )
        self.skip_blank()

    def parse_likelihood(self):
        result = self.result
        self.skip_dashed_line()
        result.data_log_like = self.parse_data(
            r'Estimated Ln Prob of Data\s*= (-?\d+[.]\d+)', float
        )
        result.mean_log_like = self.parse_data(
            r'Mean value of ln likelihood\s*= (-?\d+[.]\d+)', float
        )
        result.var_log_like = self.parse_data(
            r'Variance of ln likelihood\s*= (-?\d+[.]\d+)', float
        )
        result.mean_alpha = self.parse_data(
            r'Mean value of alpha\s*= (-?\d+[.]\d+)', float
        )
        self.skip_blank()

    def parse_fst(self):
        pattern = r'Mean value of Fst_\d+\s*= (-?\d+[.]\d+)'
        self.result.mean_fst = np.array(
            [self.parse_data(pattern, float) for _ in range(self.result.k)]
        )
        self.skip_blank()

    def parse_q(self):
        self.skip_lines(2)

        # Prepare line matcher
        regex = r'^\s*\d+\s*(\w+)\s*[(](\d+[.]?\d*)[)]\s*(\d+)\s*:\s*'
        regex += r'(\d+[.]\d+)\s*' * self.result.k
        regex += r'$'
        regex = re.compile(regex)

        def groups(x):
            return regex.match(x).groups()

        types = [str, float, int] + [float] * self.result.k

        def values(x):
            return (f(x) for f, x in zip(types, groups(x)))

        # Read mixture data
        data = []
        for _ in range(self.result.individuals):
            label, miss, pop, *q_probs = values(self.pop())
            data.append([label, miss, pop] + q_probs)

        # Save on data frame
        columns = ['label', 'missing', 'population'] + \
                  ['cluster_%s' % i for i in range(self.result.k)]
        self.result.q_matrix = pd.DataFrame(data, columns=columns)
        self.skip_blank()

    def parse_f(self):
        self.skip_lines(2)
        self.skip_blank()

        # Match each line of allele data
        regex = r'^\s*(\w+)\s*[(](\d+[.]?\d*)[)]\s*'
        regex += r'(\d+[.]?\d*)\s*' * self.result.k
        regex += r'$'
        regex = re.compile(regex)
        types = [str, float] + [float] * self.result.k

        def values(line):
            m = regex.match(line)
            if m is None:
                raise ValueError('unexpected line: %r' % line)
            return [f(x) for f, x in zip(types, m.groups())]

        # Extract data
        k = self.result.k
        data = []
        for _ in range(self.result.loci):
            name = self.pop().strip(': ')
            n_alleles = self.parse_data(r'^(\d+) alleles$', int)
            missing = self.parse_data(r'^(\d+[.]\d+)% missing data$', float)

            columns = ['name', 'total'] + ['cluster_%s' % i for i in range(k)]
            alleles = [values(self.pop()) for x in range(n_alleles)]
            alleles = pd.DataFrame(alleles, columns=columns)
            alleles.index = alleles.pop('name')

            data.append([name, n_alleles, missing, alleles])
            self.skip_blank()

        # Convert to dataframe
        columns = ['name', 'n_alleles', 'missing', 'data']
        self.result.F_matrix = pd.DataFrame(data, columns=columns)

    def parse_options(self):
        self.skip_lines(1)

        def flag(x):
            return bool(int(x))

        converters = dict(
            NUMINDS=int,
            NUMLOCI=int,
            MISSING=flag,
            id=flag,
            POPDATA=flag,
            POPFLAG=flag,
            PHENOTYPE=flag,
            EXTRACOLS=flag,
            MAXPOPS=int,
            BURNIN=int,
            NUMREPS=int,
            USEPOPINFO=flag,
            INFERALPHA=flag,
            INFERLAMBDA=flag,
            POPSPECIFICLAMBDA=flag,
            POPALPHAS=flag,
            COMPUTEPROB=flag,
            NOADMIX=flag,
            ADMBURNIN=int,
            UPDATEFREQ=flag,
            PRINTLIKES=flag,
            INTERMEDSAVE=flag,
            PRINTKLD=flag,
            PRINTNET=flag,
            PRINTLAMBDA=flag,
            ANCESTDIST=flag,
            NUMBOXES=int,
            ANCESTPINT=float,
            GENSBACK=int,
            MIGRPRIOR=float,
            PRINTQHAT=flag,
            PRINTQSUM=flag,
            ALPHA=float,
            FREQSCORR=flag,
            FPRIORMEAN=float,
            FPRIORSD=float,
            ONEFST=flag,
            LAMBDA=float,
            UNIFPRIORALPHA=flag,
            ALPHAMAX=float,
            ALPHAPRIORA=float,
            ALPHAPRIORB=float,
            ALPHAPROPSD=float,
            STARTATPOPINFO=flag,
            RANDOMIZE=flag,
            LINKAGE=flag,
            METROFREQ=int,
            REPORTHITRATE=flag,
            MARKOVPHASE=int,
            PHASED=flag,
            PLOIDY=int,
            PHASEINFO=flag,
            LOCPRIOR=flag,
            LOCPRIORINIT=float,
            LOCDATA=flag,
            LOCISPOP=flag,
            LOCPRIORSTEP=float,
            MAXLOCPRIOR=float,
            SEED=int,
        )

        self.result.options = self.parse_dict(self.pop(), converters)

    def parse_strat(self):
        def flag(x):
            return bool(int(x))

        line = self.pop().partition(':')[-1].strip()
        converters = dict(
            NUMSIMSTATS=int,
            PHENOTYPECOL=int,
            POOLFREQ=int,
            LOCUSxONLY=flag,
            EMERROR=float,
            MISSINGPHENO=int,
        )
        self.result.strat_options = self.parse_dict(line, converters)

    def parse_dict(self, line, converters={}):
        line = line.replace(',', '')
        options = line.split('\t')
        options = map(lambda x: x.strip().partition('='), options)
        options = filter(lambda x: x[0], options)
        options = dict(map((lambda x: (x[0], x[2])), options))
        return {k: converters.get(k, str)(v) for k, v in options.items()}


def parse_lines(lines, parser=StructureParser):
    """
    Use the given parser to parse the list of lines.
    """
    parser = parser(lines)
    return parser.parse()


def parse_file(file, parser=StructureParser):
    """
    Parse the given file.
    """
    return parse_lines(file.read().splitlines(), parser)

import os
from collections import deque, Counter

import numpy as np

from kpop.structure_params import MAINPARAMS_DEFAULTS, EXTRAPARAMS_DEFAULTS


def run_structure(pop, k=2, *, method='parental', job_dir=None,
                  outfile='out.structure', disp=1):
    """
    Runs the 'Structure' program for the given population.

    Args:
        pop (ndarray):
            A population with SNPs data.
        k (int):
            Number of sub-populations.
        job_dir (str):
            Location to save the mainparams, extraparams, sample and outfile
            files. Chooses a random directory if not given.
        outfile (str):
            Name of the output file.
        method:
            {'parental', 'admixture'}
        disp:
            Controls the verbosity level. (0 = no output, 1 = normal, and
            2 = verbose).
        ...

    Returns:
        A :cls:`StructureResult` object
    """

    # Prepare input files
    kwargs = {}
    if method == 'parental':
        kwargs.update(noadmix=1)
    elif method in ['admixture', 'admix']:
        kwargs.update(noadmix=0)
    else:
        raise ValueError('unknown method')
    main, extra, popfile = prepare_setup_files(pop,
                                               save_dir=job_dir,
                                               outfile=outfile,
                                               maxpops=k, **kwargs)

    # Display input files
    if disp > 2:
        print('MAINPARAMS\n' + '*' * 80, '\n' + main, '\n\n')
        print('EXTRAPARAMS\n' + '*' * 80, '\n' + extra, '\n\n')
        print('POPULATION FILE\n' + '*' * 80, '\n' + popfile, '\n\n')

    # Prepare structure command
    cmd = 'structure'
    if disp < 1:
        cmd += ' &> /dev/null'
    if disp >= 1:
        print('Running structure command\n\n    $', cmd)

    # Execute
    error = os.system(cmd)
    if error != 0:
        raise RuntimeError('structure returned with error code: %s' % error)

    # Read output file
    # (I don't know why Structure prepends an "_f" to the end of file)
    with open(outfile + '_f', encoding='utf8') as F:
        data = F.read()

    if disp >= 1:
        print(data)

    return StructureResult(data)


def structure_population(pop, *, label='ind', onerowperind=False,
                         popdata=None, popflag=None, locdata=None,
                         phenotype=None, sep=' ', linesep='\n'):
    """
    Return a string representation of a population as a structure input file.

    Args:
        pop:
            An array of individuals.
        label:
            An optional list of row. If a string is given, treats as the
            prefix to all individual names.
        onerowperind:
            Enable the ONEROWPERIND option of structure configuration, i.e.,
            each individual is represented in a single row.
        popdata, popflag, locdata, phenotype:
            Optional Structure row that are inserted in the beginning of
            each row to label individual with some characteristic.

    Example:
        >>> pop = [[[0, 1], [0, 0], [1, 1]],
        ...        [[1, 0], [0, 0], [0, 1]],
        ...        [[0, 1], [0, 0], [1, 1]],
        ...        [[1, 0], [0, 0], [0, 1]]]
        >>> beatles = ['John', 'Paul', 'George', 'Ringo']
        >>> out = structure_population(pop, label=beatles)
        >>> print(out)
        John 0 0 1
        John 1 0 1
        Paul 1 0 0
        Paul 0 0 1
        George 0 0 1
        George 1 0 1
        Ringo 1 0 0
        Ringo 0 0 1
        >>> out = structure_population(pop, label=beatles,onerowperind=True)
        >>> print(out)
        John 0 1 0 0 1 1
        Paul 1 0 0 0 0 1
        George 0 1 0 0 1 1
        Ringo 1 0 0 0 0 1
    """
    table = []
    for i, ind in enumerate(pop):
        if isinstance(label, (list, tuple)):
            row = [label[i]]
        else:
            row = [getattr(ind, 'label', label + str(i + 1))]
        table.append(row)

        # Add extra columns, if necessary
        for c in [popdata, popflag, locdata, phenotype]:
            if c is not None:
                row.append(str(c[i]))

        # One or two rows per individual?
        if onerowperind:
            row.extend(str(x) for x in ind.flatten())
        else:
            row_second = row[:]
            table.append(row_second)
            row.extend(str(x) for x in ind[:, 0])
            row_second.extend(str(x) for x in ind[:, 1])

    return linesep.join([sep.join(L) for L in table])


def params_file(data):
    """
    Creates a minimal, but functional, params file from a dictionary of
    arguments.
    """

    lst = []
    for k, v in sorted(data.items()):
        lst.append('#define %s %s' % (k.upper(), v))
    return '\n'.join(sorted(lst))


class StructureResult:
    """
    Result of computation of a Structure run.

    Attributes:
        All Structure options are saved as attributes (e.g., result.NUMINDS,
        result.MAXPOPS, etc). It also defines a few attributes.

        k:
            Number of parental/cluster populations (alias to result.NUMINDS).
        num_inds:
            Number of analyzed individuals (alias to result.MAXPOPS).
        memberships:
            An array with the membership ratios for each cluster.
        allele_divergence:
            ???
        log_prob_data:
            Logarithm of probability of yielding given data.
        log_like_mean, log_like_var:
            Mean and variance of logarithmic likelihood.
        mean_fst:
            An array with the mean FST values for each cluster compared to a
            complete parental population.
        ancestry:
            Mixture coefficients or membership ratios for each individual.
        missing_allele_freqs:
            An array with the ratio of missing alleles for each locus.
        total_allele_freqs:
            A list of Counter mappings. Each item corresponds to a specific
            locus and maps allele values to their relative frequency on the
            overall population.
        allele_freqs:
            Similar to above, but creates one list per sub-population/cluster.
    """

    NUMINDS = 0
    MAXPOPS = 0
    log_prob_data = None
    log_like_mean = None
    log_like_var = None
    alpha_mean = None
    alpha_var = None

    @property
    def k(self):
        return self.MAXPOPS

    @property
    def num_inds(self):
        return self.num_inds

    def __init__(self, body):
        self.raw = body
        body, _, tail = body.partition('Values of parameters used in '
                                       'structure:\n')
        params, strat = tail.splitlines()
        for k, v in self._parse_params(params).items():
            setattr(self, k, v)
        for k, v in self._parse_strat(strat).items():
            setattr(self, k, v)
        self._parse_body(body)

    def _parse_params(self, params):
        def convert(x):
            if x.isdigit():
                return int(x)
            try:
                return float(x)
            except ValueError:
                return x

        params = params.split(',')
        params = (x.strip().partition('=') for x in params if '=' in x)
        params = {k: convert(v) for k, _, v in params}
        return params

    def _parse_strat(self, data):
        _, _, data = data.partition(':')
        return self._parse_params(data.strip())

    def _parse_body(self, data):
        sections = {}
        lines = deque(data.splitlines())

        def skip_blank():
            while lines[0].isspace() or not lines[0]:
                lines.popleft()

        def parse_section():
            section_lines = []
            skip_blank()
            while lines:
                line = lines.popleft()
                if line.startswith('----'):
                    break
                else:
                    section_lines.append(line)
            return '\n'.join(section_lines).strip()

        # Remove initial banner (skipped)
        __empty = parse_section()
        __banner = parse_section()
        assert __banner.startswith('STRUCTURE'), __banner

        # Section describing the parameters of the simulation (skipped)
        __params = parse_section()
        assert __params.startswith('Input File:'), __params

        # Parse data sections
        self.memberships = self._parse_membership(parse_section())
        self.allele_divergence = self._parse_allele_freq_divergence(
            parse_section())
        self._parse_main_results(parse_section())

        return sections

    def _parse_membership(self, data):
        assert data.startswith('Overall proportion of membership'), data
        memberships = data.splitlines()[-1].strip().split()
        return [float(x) for x in memberships]

    def _parse_allele_freq_divergence(self, data):
        assert data.startswith('Allele-freq. divergence'), data
        return NotImplemented

    def _parse_main_results(self, data):
        lines = deque(data.splitlines())

        # Probability metrics
        names = {
            'estimated ln prob of data': 'log_prob_data',
            'mean value of ln likelihood': 'log_like_mean',
            'variance of ln likelihood': 'log_like_var',
            'mean value of alpha': 'alpha_mean',
        }
        while lines[0]:
            line = lines.popleft()
            name, _, value = line.partition('=')
            name = name.strip().lower()
            setattr(self, names[name], float(value))

        # Fst metrics
        while lines[0].isspace() or not lines[0]:
            lines.popleft()
        fst = []
        for i in range(self.MAXPOPS):
            fst.append(float(lines.popleft().partition('=')[-1]))
        self.mean_fst = np.array(fst)

        # Inferred ancestries
        while not lines[0].startswith('Inferred ancestry of individuals'):
            lines.popleft()
        lines.popleft()
        lines.popleft()
        ancestry = []
        for i in range(self.NUMINDS):
            values = lines.popleft().partition(':')[-1].split()
            values = list(map(float, values))
            ancestry.append(values)
        self.ancestry = np.array(ancestry)

        # Estimated allele frequencies for each locus
        while not lines[0].startswith('Estimated Allele Frequencies'):
            lines.popleft()
        lines.popleft()
        lines.popleft()
        data = '\n'.join(line for line in lines)
        loci = data.strip().split('\n\n')
        total_allele_freqs = []
        missing_allele_freqs = []
        allele_freqs = [[] for _ in range(self.MAXPOPS)]
        for locus in loci:
            lines = deque(locus.splitlines()[1:])
            n_alleles = int(lines.popleft().split()[0])
            missing_fraction = float(lines.popleft().split('%')[0]) / 100
            missing_allele_freqs.append(missing_fraction)
            if n_alleles != len(lines):
                raise ValueError('unexpected number of alleles: %s' % n_alleles)

            parental_freqs = Counter()
            total_allele_freqs.append(parental_freqs)
            cluster_freqs = [Counter() for _ in range(self.MAXPOPS)]
            for line in lines:
                allele_id, *tail = line.split()
                allele_id = int(allele_id)
                parental_freq, *tail = tail
                parental_freq = float(parental_freq[1:-1])
                parental_freqs[allele_id] = parental_freq
                tail = list(map(float, tail))
                for freq, cluster in zip(tail, cluster_freqs):
                    cluster[allele_id] = freq
            for cluster, locus in zip(allele_freqs, cluster_freqs):
                cluster.append(locus)
        self.missing_allele_freqs = np.array(missing_allele_freqs)
        self.total_allele_freqs = total_allele_freqs
        self.allele_freqs = allele_freqs


def prepare_setup_files(pop, save_files=True, save_dir=None,
                        label='ind', population_file='population.dat',
                        popdata=None, popflag=None, locdata=None,
                        phenotype=None,
                        **kwargs):
    """
    Creates and return the 'mainparams', 'extraparams', and 'sample' files
    for a structure job.
    """

    if not pop:
        raise ValueError('cannot proceed with empty population.')

    B = lambda x: int(bool(x))
    kwargs.setdefault('infile', population_file)
    kwargs.setdefault('numinds', len(pop))
    kwargs.setdefault('numloci', len(pop[0]))
    kwargs.setdefault('ploidy', len(pop[0][0]))
    kwargs.update(
        popdata=B(popdata),
        popflag=B(popflag),
        locdata=B(locdata),
        phenotype=B(phenotype))

    # Compute the mainparams and extraparams files
    main_params = MAINPARAMS_DEFAULTS.copy()
    extra_params = EXTRAPARAMS_DEFAULTS.copy()
    main_params.update({k: v for (k, v) in kwargs.items() if k in main_params})
    extra_params.update(
        {k: v for (k, v) in kwargs.items() if k in extra_params})
    mainparams = params_file(main_params)
    extraparams = params_file(extra_params)

    # Compute sample file
    popfile = structure_population(pop, label=label,
                                   onerowperind=main_params['onerowperind'])

    # Save files in the specified directory
    if save_files:
        if save_dir:
            save_dir = os.path.abspath(save_dir)
        else:
            save_dir = os.getcwd()
        path = lambda x: os.path.join(save_dir, x)
        with open(path('mainparams'), 'w', encoding='utf8') as F:
            F.write(mainparams)
        with open(path('extraparams'), 'w', encoding='utf8') as F:
            F.write(extraparams)
        with open(path(population_file), 'w', encoding='utf8') as F:
            F.write(popfile)

    return mainparams, extraparams, popfile

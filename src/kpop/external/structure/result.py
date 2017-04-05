from collections import deque, Counter

import numpy as np


class StructureResult:
    """
    Result of computation of a Structure run.

    Attributes:
        All Structure options are saved as attributes (e.g., result.NUMINDS,
        result.MAXPOPS, etc). It also defines a few attributes.

        k:
            Number of parental/admix populations (alias to result.NUMINDS).
        num_inds:
            Number of analyzed individuals (alias to result.MAXPOPS).
        memberships:
            An array with the membership ratios for each admix.
        allele_divergence:
            ???
        log_prob_data:
            Logarithm of probability of yielding given data.
        log_like_mean, log_like_var:
            Mean and variance of logarithmic like.
        mean_fst:
            An array with the mean FST values for each admix compared to a
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
            Similar to above, but creates one list per sub-population/admix.
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

    def __init__(self, body, job_dir=None):
        self.raw = body
        self.job_dir = job_dir
        body, _, tail = body.partition('Values of parameters used in '
                                       'structure:\n')
        params, strat = tail.splitlines()
        for k, v in self._parse_params(params).items():
            setattr(self, k, v)
        for k, v in self._parse_strat(strat).items():
            setattr(self, k, v)
        # self._parse_body(body)

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
        parse_section()
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
            'mean value of ln like': 'log_like_mean',
            'variance of ln like': 'log_like_var',
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
                raise ValueError(
                    'unexpected number of alleles: {0!s}'.format(n_alleles))

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

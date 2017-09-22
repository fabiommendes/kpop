import numpy as np

from ..prob import Prob


#
# Generic functions
#
def discard_attrs(obj, attrs):
    """
    Remove all attributes from object, if they exist.
    """

    for attr in attrs:
        if attr in obj.__dict__:
            obj.__delattr__(attr)


def get_freqs(pop):
    """
    Return a list of Prob instances representing the frequencies in each locus.
    """

    if pop._freqs is None:
        pop._freqs = pop.stats.empirical_freqs()
    return pop._freqs


def set_freqs(pop, value):
    """
    Sets the freqs attribute of a population.
    """
    if not len(value) == pop.num_loci or not \
            all(isinstance(x, Prob) for x in value):
        raise ValueError('not a valid loci frequency vector')

    pop._freqs = list(value)
    discard_attrs(pop, ['freqs_vector', 'freqs_matrix'])


def hfreqs_vector(pop):
    """
    Compute the heterozygote frequency on population.

    This is the probability of a heterozygote genotype for each locus.
    """

    if pop.ploidy == 2:
        data = np.array(pop)
        heterozygote = data[:, :, 0] != data[:, :, 1]
        return heterozygote.mean(axis=0)
    else:
        raise NotImplementedError('require diploid populations')


def parse_population_data(data):
    """
    Return population (genotype data, ids) from a string.
    """

    part_labels = [line.rpartition(':') for line in data.splitlines()]
    labels = [label or None for label, _, _ in part_labels]
    if set(labels) == {None}:
        labels = None
    data_raw = [line.split() for _, _, line in part_labels]
    data_raw = np.array([[list(e) for e in line] for line in data_raw])

    # Create a list of locus maps. Each element is a dictionary mapping
    # chars in raw data to integers
    locus_map = []
    for idx in range(data_raw.shape[1]):
        locus = data_raw[:, idx, :]
        alleles = set(locus.flat)
        alleles.discard('-')

        if all(tk.isdigit() for tk in alleles):
            map = {tk: int(tk) for tk in alleles}
        else:
            map = dict(enumerate(sorted(alleles), 1))
            map = {v: k for k, v in map.items()}
        map['-'] = 0
        locus_map.append(map)

    # Convert string data to numeric data
    data = np.zeros(data_raw.shape, dtype='uint8')
    ii, jj, kk = data.shape
    for j in range(jj):
        map = locus_map[j]
        for i in range(ii):
            for k in range(kk):
                data[i, j, k] = map[data_raw[i, j, k]]

    return data, labels


def id_from_parents(l1, l2):
    """
    Creates a new id label from parents id labels.
    """

    if l1 == l2:
        return l1
    elif l1 is None:
        return l2 + '_'
    elif l2 is None:
        return l1 + '_'

    common = []
    for c1, c2 in zip(l1, l2):
        if c1 == c2:
            common.append(c1)
            continue
        break
    common = ''.join(common)
    l1 = l1[len(common):] or '_'
    l2 = l2[len(common):] or '_'
    return '%s-%s,%s' % (common, l1, l2)


def random_individual_data(freqs, ploidy=2, seed=None):
    """
    Creates a random biallelic individual data with given ploidy.

    Freqs can be a list of (p, 1 - p) pairs or a list of p's.
    """

    if seed is not None:
        np.random.seed(seed)

    freq_array = np.asarray(freqs)
    num_loci = len(freq_array)
    values = np.random.uniform(size=num_loci * ploidy)
    values = values.reshape((num_loci, ploidy))

    if freq_array.ndim == 2:
        if freq_array.shape[1] != 2:
            raise NotImplementedError('must be biallelic!')
        return np.asarray(values >= freq_array[:, 0, None], dtype=np.uint8) + 1
    else:
        return np.asarray(values >= freq_array[:, None], dtype=np.uint8) + 1

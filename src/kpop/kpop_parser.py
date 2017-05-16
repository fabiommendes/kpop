import collections

import numpy as np

from kpop.exceptions import MissingDataError


def check_valid_names(names):
    """
    Raises ValueError if names mapping contains invalid values.
    """
    if 0 in names or '0' in names:
        raise MissingDataError('missing data should be represented by a dash')


def tokenize_locus(st):
    """
    Split string with genotypes into its components.
    """
    return st.split(',') if ',' in st else list(st)


def max_non_numeric_value(names_map):
    """
    Return the maximum value associated with a non-numeric key.
    """
    return max(v for k, v in names_map.items() if not k.isdigit())


def update_names_map(names_map, tokens):
    """
    Update names mapping with given tokens.
    """

    names_map['-'] = 0
    if '0' in tokens:
        raise MissingDataError(
            'missing data must be represented by a dash "-". Do not use zero '
            'directly.')

    # Update mapping with new tokens.
    if not all(tk in names_map for tk in tokens):
        numeric = (tk for tk in tokens if tk.isdigit())
        names_map.update({tk: int(tk) for tk in numeric})
        non_numeric = sorted(tk for tk in tokens
                             if not tk.isdigit() and tk not in names_map)
        value = max_non_numeric_value(names_map) + 1
        for tk in non_numeric:
            names_map[tk] = value
            value += 1


def str_to_data(data, names=None):
    """
    Convert string of the form "aA aa AA" to a tuple of
    (genotype_data, allele_names).

    Args:
        data:
            Input string data.
        names:
            Optional list of mappings between allele name to its numerical
            value. If a single dictionary is given, it assumes that it uses a
            common mapping to all alleles.
    """
    genotypes = data.split()

    # Normalize names
    if names is None:
        names = [{} for _ in genotypes]
    elif isinstance(names, collections.Mapping):
        mapping = dict(names)
        names = [mapping for _ in genotypes]
    else:
        names = list(names)
        if len(names) != len(genotypes):
            raise IndexError('list of names mapping do not match data size')
    for x in names:
        check_valid_names(x)

    # Save data
    result = []
    for names_map, genotype in zip(names, genotypes):
        tokens = tokenize_locus(genotype)
        update_names_map(names_map, tokens)
        locus = [names_map[tk] for tk in tokens]
        result.append(locus)
    return np.array(result), [{v: k for k, v in x.items()} for x in names]

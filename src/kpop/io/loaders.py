import csv
import operator as op
import sys
from collections import OrderedDict
from typing import List, Union, Dict

import numpy as np

from ..individual import Individual
from ..population import Population, MultiPopulation

NOT_GIVEN = object()


def load_csv(file=sys.stdout, *,
             missing='-',
             label_col: Union[str, int] = NOT_GIVEN,
             pop_col: Union[str, int] = NOT_GIVEN,
             ignore_cols: List[Union[str, int]] = None,
             ignore_rows: List[int] = None,
             meta: Dict[str, Union[str, int]] = None,
             label: str = None,
             ploidy=2):
    """
    Load population data from CSV file.

    Args:
        file:
            Name of the file or an open file-like object. If no file is given,
            uses sys.stdout and prints data on the screen.
        pop_col (optional):
            Index or name of the population column. This column identify each
            individual with its respective population. Default value for k-pop
            generated CSV files is "population".
        label_col (optional):
            Index or name for the label column. This column identifies each
            individual with a unique label within its population. Default value
            for k-pop generated CSV files is "label".
        missing:
            Character used to specify missing data. The default value is a
            single dash character per loci. Empty cells are also treated as
            missing.
        ignore_cols:
            An optional list of column indexes or labels for columns that
            should be ignored.
        ignore_rows:
            A list of rows indexes to ignore.
        meta:
            A dictionary mapping meta information labels to their corresponding
            columns.
        label:
            The label for the new resulting population.

    Returns:
        If pop_col is given, return a new MultiPopulation object. Otherwise,
        return a single Population.
    """
    if isinstance(file, str):
        with open(file) as file:
            return load_csv(**locals())

    # Read csv file and extract data
    rows, row_map = extract_rows(csv.reader(file), ignore_rows)
    header, *body = rows
    body = np.array(body, dtype=str)
    col_normalize = col_normalizer(header)
    n_ind = len(body)

    # Extract pop_values, ind_labels, etc
    pop_label = label
    pop_col = col_normalize(pop_col, 'population')
    pop_values = [None] * n_ind if pop_col is None else body[:, pop_col]
    label_col = col_normalize(label_col, 'label')
    label_values = [None] * n_ind if label_col is None else body[:, label_col]
    meta = {k: col_normalize(v) for k, v in (meta or {}).items()}

    # Compute ignore cols and remove them from the body
    ignore_cols = {col_normalize(x) for x in (ignore_cols or ())}
    ignore_cols.update(meta.values())
    ignore_cols.update([x for x in [pop_col, label_col] if x is not None])
    body, col_map = extract_cols(body, ignore_cols)

    # Normalize missing data
    body[body == ''] = '--'
    shape = body.shape
    norm_missing = op.methodcaller('replace', missing, '-')
    body = np.array(list(map(norm_missing, body.flat)))
    body = body.reshape(shape)

    # Assert all cells have the same length
    try:
        ploidy, = set(map(len, body.flat))
    except ValueError:
        raise inconsistent_length_error(body, row_map, col_map)

    # Convert body to a 3D array
    body = np.array([list(x) for x in body.flat]).reshape(shape + (2,))

    # Get all values for each column and prepare numerical conversions
    new_body = np.zeros(shape + (2,), dtype='uint8')
    for i in range(body.shape[1]):
        col_data = body[:, i].flatten()
        genotypes = set(np.unique(col_data))
        genotypes.discard('-')
        gmap = dict((y, x) for x, y in enumerate(sorted(genotypes), 1))
        gmap['-'] = 0
        col_data = np.array(list(map(gmap.get, col_data)), dtype='uint8')
        new_body[:, i] = col_data.reshape((n_ind, 2))

    # Recreate populations from body data
    body = new_body
    populations = OrderedDict((p, []) for p in pop_values)
    for i, (label, pop, row) in enumerate(zip(label_values, pop_values, body)):
        ind = Individual(row, label=label)
        populations[pop].append(ind)

    # Create population or multipopulation return value
    if pop_col is None:
        result = Population(populations[None], label=pop_label)
    else:
        def factory(x): return Population(x[1], label=x[0])
        populations = map(factory, populations.items())
        result = MultiPopulation(populations, label=pop_label)
    result.loci_names = header
    return result


def load_pickle(file):
    """
    Load population from file.

    This method is useful to load populations saved with the pop.save() method.
    Pickle is an internal binary protocol used in Python. It is very  efficient
    method, but it is not portable to other languages.
    """

    return Population.load(file, 'pickle')


def load(file):
    """
    Loads population from file.

    It automatically recognizes the format of the input file using its
    name/extension.
    """

    name = file if isinstance(file, str) else file.name
    if name.endswith('.csv'):
        return load_csv(file)
    elif name.endswith('.pickle'):
        return load_pickle(file)
    else:
        raise ValueError('could not determine file type for %r' % file)


#
# Auxiliary functions
#
def extract_rows(csv, ignore):
    """
    Extract rows from csv file skipping ones that should be ignored and return
    a tuple of (data, row_map) where row_map is a dictionary of new indexes
    to the old indexes.
    """

    ignore = set(ignore or ())
    rows = []
    row_map = {}
    diff = 0
    for idx, row in enumerate(csv):
        if (idx in ignore) or len(row) == 0:
            diff += 1
        else:
            rows.append(row)
            row_map[idx - diff] = idx
    return rows, row_map


def extract_cols(data, ignore):
    """
    Extract cols from numpy array skipping ones that should be ignored and
    return a tuple of (data, col_map) where col_map is a dictionary of new
    indexes to the old indexes.
    """

    ignore = set(ignore or ())
    result = []
    col_map = {}
    diff = 0
    for idx, col in enumerate(data.T):
        if idx in ignore:
            diff += 1
        else:
            result.append(col)
            col_map[idx - diff] = idx
    return np.array(result).T, col_map


def col_normalizer(head):
    """
    Return a function that normalizes values to indexes of the head object.
    """

    def func(x, default=None, optional=True):
        if x is NOT_GIVEN:
            x = default
        if x is None:
            return None
        elif isinstance(x, int):
            return x
        else:
            try:
                return head.index(x)
            except ValueError:
                if optional:
                    return None
                raise
    return func


def print_csv(csv):
    for row in csv:
        print(', '.join(map(str, row)))


def inconsistent_length_error(data, row_map, col_map):
    n = len(data[0, 0])
    for i, row in enumerate(data):
        for j, cell in enumerate(row):
            if len(cell) != n:
                i, j = row_map[i], col_map[j]
                return ValueError('bad data at (%s, %s): %r' % (i, j, cell))

    return RuntimeError

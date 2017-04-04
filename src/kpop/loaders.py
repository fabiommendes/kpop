import csv
from collections import OrderedDict

from kpop import Individual, Population, MultiPopulation

NOT_GIVEN = object()


def load_csv(file, label_col=NOT_GIVEN, pop_col=NOT_GIVEN, missing='-',
             ignore=None, meta=None, label=None):
    """
    Load population data from CSV file.

    Args:
        file:
            Name of the file or an open file-like object.
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
        ignore:
            An optional list of column indexes or labels for columns that
            should be ignored.
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
    rows = iter(csv.reader(file))

    header = next(rows)
    index_map = {name: idx for idx, name in enumerate(header)}

    # Normalize pop_label index
    if pop_col is NOT_GIVEN and 'population' in index_map:
        pop_col = index_map['population']
    elif pop_col is not NOT_GIVEN:
        pop_col = pop_col if isinstance(pop_col, int) else index_map[pop_col]
    else:
        pop_col = None

    # Normalize label index
    if label_col is NOT_GIVEN and 'label' in index_map:
        label_col = index_map['label']
    elif label_col is not NOT_GIVEN:
        label_col = label_col if isinstance(label_col, int) else index_map[
            label_col]
    else:
        label_col = None

    # Normalize skip and meta indexes
    if meta is not None:
        meta = {k: (v if isinstance(v, int) else index_map[v])
                for k, v in meta.items()}
    ignore = {x if isinstance(x, int) else index_map[x] for x in ignore or ()}
    if pop_col is not None:
        ignore.add(pop_col)
    if label_col is not None:
        ignore.add(label_col)
    if meta:
        ignore.update(meta.values())

    # Process lines
    populations = OrderedDict()
    for row in rows:
        # Ignore blank lines
        if not row:
            continue

        pop = row[pop_col] if pop_col is not None else None
        ind_label = row[label_col] if label_col is not None else None
        data = [cell.replace(missing, '-') or '--'
                for i, cell in enumerate(row) if i not in ignore]
        ind = Individual(' '.join(data), label=ind_label)
        try:
            populations[pop].append(ind)
        except KeyError:
            populations[pop] = [ind]

    # Loci names
    names = [x for i, x in enumerate(header) if i not in ignore]

    if pop_col is None:
        result = Population(populations[None], label=label)
    else:
        populations = [Population(data, label=label) for label, data in
                       populations.items()]
        result = MultiPopulation(populations, label=label)
    result.loci_names = names
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

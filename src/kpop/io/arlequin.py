import itertools
import sys
from typing import List

import numpy as np

from ..population import Population


def export_arlequin(pop: Population,
                    file=sys.stdout, *,
                    title: str = None,
                    names: List[str] = None,
                    missing: str = '-'):
    """
    Export a population to an Arlequin file.

    Args:
        pop:
            A MultiPopulation object with several sub-populations.
        file:
            A file object or file path. If not given, uses sys.stdout and writes
            output to the console.
        names:
            A list of sub-population names to override the default population
            labels.
        title:
            A name to override the default multipopulation label.
        missing:
            The missing data label.
    """
    # Handle string files
    if isinstance(file, str):
        with open(file, 'w') as file:
            return export_arlequin(**locals())

    create_header(pop, file, title, missing)
    file.write('\n\n')
    create_data(pop, file, names, missing)


#
# Auxiliary functions.
#
def create_header(pop, file=sys.stdout, title=None, missing='-'):
    title = title or getattr(pop, 'label', None) or 'Sample Arlequin data'
    file.write('[Profile]\n')
    file.write('Title="%s"\n' % title)
    file.write('NbSamples=%s\n' % len(pop.populations))
    file.write('DataType=STANDARD\n')
    file.write('GenotypicData=1\n')
    file.write('LocusSeparator=WHITESPACE\n')
    file.write('GameticPhase=0\n')
    file.write('RecessiveData=0\n')
    file.write("MissingData='%s'\n" % missing)


def create_sample_data(subpop, file, missing):
    labels = (ind.id or 'ind%s' % i for i, ind in enumerate(subpop))
    labels = ['    %s ' % label for label in labels]
    label_size = len(max(labels, key=len))
    indent = ' ' * (label_size + 2)

    for label, ind in zip(labels, subpop):
        # First line
        file.write(label.ljust(label_size))
        file.write('1 ')
        file.write(' '.join(render(ind[:, 0], missing=missing)))
        file.write('\n')

        # Second line
        file.write(indent)
        file.write(' '.join(render(ind[:, 0], missing=missing)))
        file.write('\n')


def create_data(pop, file=sys.stdout, names=None, missing='-'):
    names = itertools.repeat(None) if names is None else names
    names = iter(names)

    file.write('[Data]\n')
    file.write('[[Samples]]\n')

    for idx, subpop in enumerate(pop.populations):
        name = next(names) or subpop.id or 'pop%s' % idx
        file.write('SampleName="%s"\n' % name.replace('"', ''))
        file.write('SampleSize=%s\n' % len(subpop))
        file.write('SampleData={\n')
        create_sample_data(subpop, file, missing)
        file.write('}\n')
        if idx != len(pop) - 1:
            file.write('\n')


def render(data, missing='-'):
    """
    Convert array of numbers to an array of strings and converts all zeros
    to the given missing data char.
    """

    data = np.asarray(data)
    missing = missing or '-'
    missing_idxs = data == 0
    new = data.astype(str)
    new[missing_idxs] = missing
    return new

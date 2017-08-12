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
def create_header(pop, F=sys.stdout, title=None, missing='-'):
    title = title or getattr(pop, 'label', None) or 'Sample Arlequin data'
    F.write('[Profile]\n')
    F.write('Title="%s"\n' % title)
    F.write('NbSamples=%s\n' % len(pop.populations))
    F.write('DataType=STANDARD\n')
    F.write('GenotypicData=1\n')
    F.write('LocusSeparator=WHITESPACE\n')
    F.write('GameticPhase=0\n')
    F.write('RecessiveData=0\n')
    F.write("MissingData='%s'\n" % missing)


def create_sample_data(subpop, F, missing):
    labels = (ind.label or 'ind%s' % i for i, ind in enumerate(subpop))
    labels = ['    %s ' % label for label in labels]
    label_size = len(max(labels, key=len))
    indent = ' ' * (label_size + 2)

    for label, ind in zip(labels, subpop):
        # First line
        F.write(label.ljust(label_size))
        F.write('1 ')
        F.write(' '.join(render(ind[:, 0], missing=missing)))
        F.write('\n')

        # Second line
        F.write(indent)
        F.write(' '.join(render(ind[:, 0], missing=missing)))
        F.write('\n')


def create_data(pop, F=sys.stdout, names=None, missing='-'):
    names = itertools.repeat(None) if names is None else names
    names = iter(names)

    F.write('[Data]\n')
    F.write('[[Samples]]\n')

    for idx, subpop in enumerate(pop.populations):
        name = next(names) or subpop.label or 'pop%s' % idx
        F.write('SampleName="%s"\n' % name.replace('"', ''))
        F.write('SampleSize=%s\n' % len(subpop))
        F.write('SampleData={\n')
        create_sample_data(subpop, F, missing)
        F.write('}\n')
        if idx != len(pop) - 1:
            F.write('\n')


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

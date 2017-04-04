import os

import click

from .util import printif, verbose, file_to_pop


@click.command()
@click.option('--format', '-f', prompt='structure, admixture, or arlequin?',
              help='output format (structure, admixture, or arlequin)')
@verbose()
@file_to_pop()
def export(pop, format=None, output=None, verbose=False):
    """
    Export population to external formats.
    """

    # Normalize
    if output is None:
        output = os.path.splitext(pop.file)[0]

    if format == 'structure':
        return export_structure(pop, output, verbose)
    elif format == 'admixture':
        return export_admixture(pop, output)
    elif format == 'arlequin':
        return export_arlequin(pop, output)
    else:
        raise ValueError('invalid format: {0!r}'.format(format))


def export_structure(pop, output, verbose):
    from ..external.structure import structure_population, mainparams, \
        extraparams

    # Load folder
    folder = os.path.splitext(output)[0]
    printif(verbose, 'Exported data will be saved at ', nl=False)
    printif(verbose, folder + os.sep, bold=True)

    if not os.path.exists(folder):
        os.mkdir(folder)
    if not os.path.isdir(folder):
        raise SystemExit('{0!s}{1!s} is not a folder'.format(folder, os.sep))

    # Save params files
    printif(verbose, 'Creating param files')
    with open(os.path.join(folder, 'mainparams'), 'w') as F:
        F.write(mainparams(pop, outfile='result'))
    with open(os.path.join(folder, 'extraparams'), 'w') as F:
        F.write(extraparams(pop))

    # Create database
    printif(verbose, 'Converting database')
    data = structure_population(pop, onerowperind=True)
    with open(os.path.join(folder, 'data.txt'), 'w') as F:
        F.write(data)
        F.write('\n')


def export_arlequin(pop, output):
    pass


def export_admixture(pop, output):
    pass

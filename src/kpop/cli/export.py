import os

import click

from .util import printif, verbose, file_to_pop


@click.command()
@click.option('--format', '-f', prompt='structure, admixture, or arlequin?',
              help='output format (structure, admixture, or arlequin)')
@click.option('--output', '-o', default=None,
              help='output file')
@click.option('--id', '-i', default=None,
              help='a string identifier for the population')
@verbose()
@file_to_pop()
def export(pop, format=None, output=None, verbose=False, id=None):
    """
    Export population to external formats.
    """

    # Normalize
    if output is None:
        output = os.path.splitext(pop.file)[0]

    if format == 'structure':
        return export_structure(pop, output, verbose)
    elif format == 'admixture':
        return export_admixture(pop, output, id)
    elif format == 'arlequin':
        return export_arlequin(pop, output, id)
    else:
        raise SystemExit('invalid format: %r' % format)


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
        raise SystemExit('%s%s is not a folder' % (folder, os.sep))

    # Save params files
    printif(verbose, 'Creating param files')
    with open(os.path.join(folder, 'mainparams'), 'w') as F:
        F.write(mainparams(pop, outfile='result'))
    with open(os.path.join(folder, 'extraparams'), 'w') as F:
        F.write(extraparams())

    # Create database
    printif(verbose, 'Converting database')
    data = structure_population(pop, onerowperind=True)
    with open(os.path.join(folder, 'data.txt'), 'w') as F:
        F.write(data)
        F.write('\n')


def export_arlequin(pop, output, title):
    from kpop.io import export_arlequin

    export_arlequin(pop, output + '.arp', title=title)


def export_admixture(pop, output):
    pass

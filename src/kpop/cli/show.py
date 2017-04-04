import click

from .util import verbose, file_to_pop


@click.command()
@click.option('--method', '-m', prompt='method', help='visualization method')
@verbose()
@file_to_pop()
def show(pop, method, verbose=False):
    """
    Dimensionality reduction and visualization.
    """

    if method == 'pca':
        pop.plot.pca()
    else:
        raise SystemExit('invalid method: %r' % method)

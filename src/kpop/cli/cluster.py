import click

from .util import verbose, file_to_pop


@click.command()
@click.option('--method', '-m', prompt='method', help='clusterization method')
@verbose()
@file_to_pop()
def cluster(pop, method, verbose=False):
    """
    Clusterization methods.
    """

    raise NotImplementedError('sorry, not implemented yet ;)')
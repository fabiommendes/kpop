import click

from .util import verbose, file_to_pop


@click.command()
@click.option('--method', '-m', prompt='method', help='clusterization method')
@verbose()
@file_to_pop()
def admix(pop, method, verbose=False):
    """
    Admixture coefficients and soft clustering.

    Also implements some fuzzy clustering methods.
    """

    raise NotImplementedError('sorry, not implemented yet ;)')

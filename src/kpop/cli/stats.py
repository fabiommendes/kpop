import click

from .util import verbose, file_to_pop


@click.command()
@click.option('--json', '-j', is_flag=True, help='print results JSON instead of YAML')
@verbose()
@file_to_pop()
def stats(pop, verbose=False, json=False):
    """
    Statistics about population.
    """

    click.echo('summary:')
    click.echo('    size: {0!s}'.format(pop.size))
    click.echo('    num_loci: {0!s}'.format(pop.num_loci))
    click.echo('    ploidy: {0!s}'.format(pop.ploidy))
    click.echo('    num_alleles: {0!s}'.format(pop.num_alleles))
    click.echo('    missing_ratio: {0:.3e}'.format(pop.missing_ratio))

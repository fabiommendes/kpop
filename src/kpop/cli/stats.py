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
    click.echo('    size: %s' % pop.size)
    click.echo('    num_loci: %s' % pop.num_loci)
    click.echo('    ploidy: %s' % pop.ploidy)
    click.echo('    num_alleles: %s' % pop.num_alleles)
    click.echo('    missing_ratio: %.3e' % pop.missing_ratio)

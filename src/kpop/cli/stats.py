import click

from .util import verbose, file_to_pop


@click.command()
@click.option('--json', '-j', is_flag=True,
              help='print results JSON instead of YAML')
@click.option('--sub-populations', '-s', is_flag=True,
              help='detailed stats for each sub-population.')
@verbose()
@file_to_pop()
def stats(pop, verbose=False, json=False, sub_populations=False):
    """
    Statistics about population.
    """

    click.echo(click.style('summary:', bold=True))
    pop_summary(pop)
    if sub_populations:
        click.echo(click.style('populations:', bold=True))
        for i, subpop in enumerate(pop.populations, 1):
            entry('- %s' % (subpop.label or 'pop%s' % i))
            pop_summary(subpop, indent=3, label=False)


def pop_summary(pop, indent=1, label=True):
    if label and pop.label:
        entry('label', pop.label, indent)
    entry('size', pop.size, indent)
    entry('num_loci', pop.num_loci, indent)
    entry('ploidy', pop.ploidy, indent)
    entry('num_alleles', pop.num_alleles, indent)
    entry('missing_ratio', pop.missing_ratio, indent)
    entry('num_populations', pop.num_populations, indent)


def entry(k, v='', indent=1):
    click.echo('  ' * indent + click.style(k + ': ', bold=True) + str(v))

import click

from kpop.population.population_base import PopulationBase
from .util import verbose, file_to_pop


@click.command()
@click.option('--json', '-j', is_flag=True,
              help='print results JSON instead of YAML')
@click.option('--sub-populations', '-s', is_flag=True,
              help='detailed stats for each sub-population.')
@click.option('--fst', '-f', is_flag=True,
              help='Print the pairwise F_st statistics for sub-populations.')
@click.option('--freqs-biallelic', '-fb', is_flag=True,
              help='Prints a table with frequencies of the first allele.')
@click.option('--sep', default='  ',
              help='String separator for csv-like frequency data.')
@click.option('--csv', is_flag=True,
              help='renders table as CSV data.')
@verbose()
@file_to_pop()
def stats(pop: PopulationBase, verbose=False, json=False,
          sub_populations=False, fst=False,
          freqs_biallelic=False, sep='  ', csv=False):
    """
    Statistics about population.
    """

    if fst:
        click.echo(pop.stats.render_fst())
        return
    elif freqs_biallelic:
        kwargs = {'sep': sep}
        if csv:
            kwargs['align'] = None
            kwargs['sep'] = sep.strip() or ','
        click.echo(pop.stats.render_biallelic_freqs(**kwargs))
        return
    else:
        click.echo(click.style('summary:', bold=True))
        pop_summary(pop)
        if sub_populations:
            click.echo(click.style('populations:', bold=True))
            for i, subpop in enumerate(pop.populations, 1):
                entry('- %s' % (subpop.id or 'pop%s' % i))
                pop_summary(subpop, indent=3, id=False)


def pop_summary(pop, indent=1, id=True):
    if id and pop.id:
        entry('label', pop.id, indent)
    entry('size', pop.size, indent)
    entry('num_loci', pop.num_loci, indent)
    entry('ploidy', pop.ploidy, indent)
    entry('num_alleles', pop.num_alleles, indent)
    entry('missing_ratio', pop.missing_data_ratio, indent)
    entry('num_populations', pop.num_populations, indent)


def entry(k, v='', indent=1):
    click.echo('  ' * indent + click.style(k + ': ', bold=True) + str(v))

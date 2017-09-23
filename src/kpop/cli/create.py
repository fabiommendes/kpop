import string

import click

from .. import Population, MultiPopulation

SYMBOL_TABLE = string.ascii_uppercase


@click.command()
@click.argument('file')
@click.option('--size', '-s', type=int, prompt='size',
              help='population size')
@click.option('--num-loci', '-l', type=int, prompt='num loci',
              help='number of loci')
@click.option('--sub-populations', '-k', type=int, default=1,
              help='number of sub populations')
@click.option('--id', '-i', default='random',
              help='a string identifier for the population')
def create(file, size=100, num_loci=100, sub_populations=1, id='random'):
    """
    Random synthetic populations.
    """

    pops = []
    for i in range(sub_populations):
        symbol = get_pop_symbol(i, sub_populations)
        pop_id = id if sub_populations == 1 else '%s-%s' % (id, symbol)
        pop = Population.random(size, num_loci, id=pop_id)
        pops.append(pop)

    if sub_populations == 1:
        pop_total = pops[0]
    else:
        pop_total = MultiPopulation(pops)
    pop_total.io.save(file, format='auto')


def get_pop_symbol(i, sub_populations):
    if sub_populations > 26:
        return 'pop%s' % (i + 1)
    else:
        return 'pop' + SYMBOL_TABLE[i]

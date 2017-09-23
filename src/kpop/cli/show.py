import click

from kpop.libs import plt
from .util import verbose, file_to_pop


@click.command()
@click.option(
    '--method', '-m', prompt='method',
    help='visualization method'
)
@click.option(
    '--alpha', '-a', default=1.0,
    help='opacity in the [0, 1] range'
)
@click.option(
    '--no-legend', default=False, is_flag=True,
    help='disable legends',
)
@click.option(
    '--output', '-o',
    help='saves to output file',
)
@verbose()
@file_to_pop()
def show(pop, method, verbose=False, alpha=1, no_legend=False, output=None):
    """
    Dimensionality reduction and visualization.
    """

    from kpop.population.projection import Projection

    method_ = method.replace('-', '_').replace(' ', '_').lower()
    if method_ in Projection._methods:
        getattr(pop.plot, method_)(alpha=alpha, legend=not no_legend)
        if output:
            plt.savefig(output)
        else:
            plt.show()
    else:
        raise SystemExit('invalid method: %r' % method)

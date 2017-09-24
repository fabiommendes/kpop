import click
from colorama import init

from .create import create
from .export import export
from .import_ import import_
from .shell import shell
from .show import show
from .stats import stats


# from .cluster import cluster
# from .admix import admix


@click.group()
def cli():
    """
    Kpop software for population genetics.
    """
    pass


def main():
    """
    Main entry point for kpop console interface.
    """

    init()
    # cli.add_command(admix)
    # cli.add_command(cluster)
    cli.add_command(create)
    cli.add_command(export)
    cli.add_command(shell)
    cli.add_command(show)
    cli.add_command(stats)
    cli.add_command(import_)
    cli()


if __name__ == '__main__':
    main()

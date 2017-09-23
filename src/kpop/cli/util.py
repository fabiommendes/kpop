import functools

import click

from .. import load


def verbose():
    """
    Decorator that marks the '--verbose' option to function.
    """

    return click.option('--verbose', '-v', is_flag=True,
                        help='print detailed information')


def file_to_pop():
    """
    Adds a file argument to function and convert it to a Population object.
    """

    def decorator(func):
        @click.argument('file')
        @functools.wraps(func)
        def decorated(file, *args, **kwargs):
            verbose = kwargs.get('verbose', False)
            pop = load(file)
            pop.file = file
            printif(verbose, 'Population loaded successfully!')
            return func(pop, *args, **kwargs)

        return decorated

    return decorator


def printif(cond, *args, **kwargs):
    """
    Prints if first argument is True.
    """

    if cond:
        click.secho(*args, **kwargs)

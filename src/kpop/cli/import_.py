import os

import click

from kpop.io import load_csv
from kpop.io import load_pickle


@click.command('import')
@click.argument('file')
@click.option('--output', '-o', default=None,
              help='output file (.csv or .pickle)')
@click.option('--pop-col', '-p', default=None,
              help='index or name of the population column')
@click.option('--ind-col', '-i', default=None,
              help='index or name of the individual labels column')
@click.option('--missing-data', '-m', default='-',
              help='character used to represent missing data')
@click.option('--skip-cols', '-c', default=None,
              help='a comma separated list of columns to ignore')
@click.option('--skip-rows', '-r', default=None,
              help='a comma separated list of rows to ignore')
@click.option('--label', '-l', default=None,
              help='a string label that describes population')
@click.option('--meta', default=None,
              help='a comma separated list of meta:col pairs mapping meta '
                   'information with their respective columns')
@click.option('--debug/--no-debug', '-d', default=False,
              help='if true, shows debug information')
def import_(file, output=None, id=None,
            pop_col=None, ind_col=None, missing_data='-',
            skip_rows=None, skip_cols=None,
            meta=None, debug=False):
    """
    Import file to kpop's .csv or .pickle formats.
    """
    name, ext = os.path.splitext(file)

    # Load data from csv
    if ext == '.csv':
        kwargs = dict(
            pop_col=number_or_text(pop_col, -1),
            label_col=number_or_text(ind_col, -1),
            missing=missing_data,
            ignore_cols=normalize_list(skip_cols),
            ignore_rows=normalize_list(skip_rows),
            meta=meta,
            id=label
        )
        kwargs = {k: v for k, v in kwargs.items() if v is not None}
        try:
            pop = load_csv(file, **kwargs)
        except ValueError as ex:
            if debug:
                raise

            error(ex)
            hint('Call this command with --help for more options.')
            raise SystemExit(1)

    # Load data from pickle
    elif ext == '.pickle':
        pop = load_pickle(file)

    # Invalid input data
    elif not ext:
        raise SystemExit('unknown input file format')
    else:
        raise SystemExit('unknown file format: %r' % ext)

    # Saves loaded population
    if output:
        pop.save(output)
    else:
        click.echo(click.style('Population data', bold=True))
        click.echo('\n' + repr(pop) + '\n')
        click.echo(click.style('Hint: ', bold=True), nl=False)
        click.echo('append `-o FILE` option to determine the output file.')
        click.echo('    You may want to use the .pickle extension.')


def error(msg):
    error_prompt = click.style('Error: ', fg='red', bold=True, blink=True)
    click.echo(error_prompt, nl=False)
    click.echo(msg)


def hint(msg):
    click.echo(click.style('Hint: ', bold=True), nl=False)
    click.echo(msg)


def number_or_text(x, inc=0):
    if x is None:
        return None
    try:
        return int(x) + inc
    except ValueError:
        pass
    try:
        return float(x) + inc
    except ValueError:
        pass
    return str(x)


def normalize_list(seq):
    if seq is None:
        return None
    return [number_or_text(x, -1) for x in seq.split(',')]

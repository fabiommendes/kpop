import os
import re
import sys
from collections import defaultdict

import click

from kpop import __version__, load
from .util import verbose, printif

# Module constants
NAME_RE = re.compile(r'^[a-zA-Z_]\w*$')
BANNER = r'''
Kpop is powered by IPython interactive shell.
    Kpop %s.
    Python %s.
 _ __
| / /___  ___ ___
|  \| . \/ . \ . \
|_\_\  _/\___/  _/
    |_|      |_|

Population genetics and machine learning. Please refer to https://readthedocs.org/kpop
for the complete documentation.
'''.lstrip('\n') % (__version__, sys.version.replace('\n', ' '))


@click.command()
@click.option('--notebook', '-n', is_flag=True, help='starts a notebook')
@verbose()
def shell(verbose=False, notebook=False):
    """
    Start an IPython shell.
    """

    if notebook:
        start_notebook(verbose)
    else:
        start_ipython_shell(verbose)


def start_ipython_shell(verbose):
    from IPython.terminal.embed import InteractiveShellEmbed as _Shell

    class Shell(_Shell):
        banner = BANNER

    namespace = start_shell_namespace(os.getcwd(), verbose)
    exec('from kpop.all import *', namespace)
    shell = Shell.instance()
    shell(local_ns=namespace)


def start_notebook(verbose, filename='kpop.ipynb'):
    from notebook.notebookapp import NotebookApp
    from jupyter_client import KernelManager

    file_path = os.path.join(os.getcwd(), filename)
    if not os.path.exists(file_path):
        create_default_notebook_file(file_path)

    sys.argv[:] = ['notebook', 'kpop.ipynb']

    NotebookApp.launch_instance(
        kernel_manager=KernelManager,
        kernel_name='python',
    )


def create_default_notebook_file(file_path):
    """
    Creates an empty kpop-powered notebook file.
    """
    import nbformat

    dirname = os.path.dirname(file_path)
    data = {
        'metadata': {
            'kernel_info': {
                'name': 'python'
            },
        },
        'nbformat': 4,
        'nbformat_minor': 0,
        'cells': [
            {
                'cell_type': 'markdown',
                'metadata': {},
                'source': (
                    '# Kpop notebook\n'
                    'Execute next cell to initialize kpop and all populations.'),
            },
            {
                'cell_type': 'code',
                'execution_count': 0,
                'metadata': {},
                'source': get_first_cell_source_for_notebook(dirname),
                'outputs': []
            },
        ],
    }
    with open(file_path, 'w') as F:
        nb = nbformat.from_dict(data)
        nbformat.write(nb, F, version=4)


def get_first_cell_source_for_notebook(path):
    """
    Create first cell of code that should be displayed on new notebooks.
    """
    names_to_exts = get_name_to_exts_map(path)
    lines = ['from kpop.all import *', '']

    if not names_to_exts:
        lines.append('# no populations were found!')

    # Fill lines with population values
    for file, exts in names_to_exts.items():
        ext = load_best_ext(exts)
        name = python_name(file)
        file = os.path.join(path, file + ext)
        lines.append('%s = load(%r)' % (name, file))

    return '\n'.join(lines)


def start_shell_namespace(path, verbose=False):
    """
    Load all population files in the given path.

    Returns:
        A dictionary mapping population file (without extension) to the
        corresponding population.
    """

    names_to_exts = get_name_to_exts_map(path)
    result = {}

    printif(verbose, click.style('Loading populations:', bold=True, fg='green'))
    for file, exts in names_to_exts.items():
        ext = load_best_ext(exts)
        name = python_name(file)
        fpath = os.path.join(path, file + ext)

        try:
            result[name] = load(fpath)
        except:
            click.echo(
                click.style('Invalid input: ', fg='red', bold=True) +
                click.style(repr(file + ext), fg='green')
            )
        else:
            click.echo(click.style(name, bold=True, fg='green'), nl=False)
            click.echo(' from %r' % (file + ext))

    click.echo('')
    return result


def get_name_to_exts_map(path):
    """
    Map file names name to available extensions.
    """
    valid_extensions = {'.csv', '.pickle'}

    # Map file names to possible extensions
    files = os.listdir(path)
    names_to_exts = defaultdict(set)
    for fname in files:
        fname, ext = os.path.splitext(fname)
        if ext in valid_extensions:
            names_to_exts[fname].add(ext)

    return names_to_exts


def load_best_ext(exts):
    """
    Return the best extension from a sequence of extensions.
    """

    if not exts:
        raise ValueError('empty sequence')
    if len(exts) == 1:
        return next(iter(exts))
    if '.pickle' in exts:
        return '.pickle'
    elif '.csv' in exts:
        return '.csv'
    else:
        return next(iter(exts))


def python_name(name):
    """
    Normalize string to a valid Python name.
    """

    old_name = name
    transforms = [
        lambda x: x,
        lambda x: x.replace('-', '')
    ]
    for f in transforms:
        name = f(name)
        if NAME_RE.match(name):
            return name
    else:
        raise ValueError('canot convert to valid python name: %r' % old_name)

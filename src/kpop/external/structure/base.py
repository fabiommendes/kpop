import os
import tempfile

import subprocess

from kpop.population.population_base import PopulationBase
from .result import StructureResult
from .params import MAINPARAMS_DEFAULTS, EXTRAPARAMS_DEFAULTS


def run_structure(pop, k=2, *, method='parental', job_dir=None,
                  outfile='out.structure', disp=1, keep_dir=False):
    """
    Runs the 'Structure' program for the given population.

    Args:
        pop (ndarray):
            A population with SNPs data.
        k (int):
            Number of sub-populations.
        job_dir (str):
            Location to save the mainparams, extraparams, sample and outfile
            files. Chooses a random directory if not given.
        outfile (str):
            Name of the output file.
        method:
            {'parental', 'admixture'}
        disp:
            Controls the verbosity level. (0 = no output, 1 = normal, and
            2 = verbose).
        keep_dir:
            If True, do not remove the temporary build directory. If job_dir
            is explicitly given, it is never removed.

    Returns:
        A :class:`StructureResult` object
    """

    # Create temporary job directory
    if job_dir is None and keep_dir:
        job_dir = tempfile.mkdtemp(prefix='structure-')
    elif job_dir is None:
        tempdir = tempfile.TemporaryDirectory(prefix='structure-')
        job_dir = tempdir.name

    # Prepare input files
    kwargs = {}
    if method == 'parental':
        kwargs.update(noadmix=1)
    elif method in ['admixture', 'admix']:
        kwargs.update(noadmix=0)
    else:
        raise ValueError('unknown method')
    main, extra, popfile = prepare_setup_files(pop,
                                               save_dir=job_dir,
                                               outfile=outfile,
                                               maxpops=k, **kwargs)

    # Display input files
    if disp >= 1:
        print('Started structure analysis at %r' % job_dir)
        for file in os.listdir(job_dir):
            size = os.path.getsize(os.path.join(job_dir, file))
            print('    %s (%.2f kb)' % (file, size / 1024))
    if disp >= 2:
        print()
        print('MAINPARAMS\n' + '-' * 80, '\n' + main, '\n\n')
        print('EXTRAPARAMS\n' + '-' * 80, '\n' + extra, '\n\n')
        print('POPULATION FILE\n' + '-' * 80, '\n' + popfile, '\n\n')

    # Prepare structure command
    cmd = 'structure'
    if disp < 1:
        cmd += ' &> /dev/null'
    if disp >= 1:
        print('Running structure command\n\n    $', cmd)

    # Execute
    error = subprocess.check_call(cmd, shell=True, cwd=job_dir)
    if error != 0:
        raise RuntimeError('structure returned with error code: %s' % error)

    # Read output file
    # (I don't know why Structure prepends an "_f" to the end of file)
    outfile = os.path.join(job_dir, outfile + '_f')
    with open(outfile, encoding='utf8') as F:
        data = F.read()

    if disp >= 2:
        print('\nOUTFILE')
        print('-' * 80)
        print(data)

    return StructureResult(data, job_dir=job_dir)


def structure_population(pop, *, id='ind', onerowperind=False,
                         popdata=None, popflag=None, locdata=None,
                         phenotype=None, sep=' ', linesep='\n'):
    """
    Return a string representation of a population as a structure input file.

    Args:
        pop:
            An array of individuals.
        label:
            An optional list of row. If a string is given, treats as the
            prefix to all individual names.
        onerowperind:
            Enable the ONEROWPERIND option of structure configuration, i.e.,
            each individual is represented in a single row.
        popdata, popflag, locdata, phenotype:
            Optional Structure row that are inserted in the beginning of
            each row to label individual with some characteristic.

    Example:
        >>> pop = [[[0, 1], [0, 0], [1, 1]],
        ...        [[1, 0], [0, 0], [0, 1]],
        ...        [[0, 1], [0, 0], [1, 1]],
        ...        [[1, 0], [0, 0], [0, 1]]]
        >>> beatles = ['John', 'Paul', 'George', 'Ringo']
        >>> out = structure_population(pop, id=beatles)
        >>> print(out)
        John 0 0 1
        John 1 0 1
        Paul 1 0 0
        Paul 0 0 1
        George 0 0 1
        George 1 0 1
        Ringo 1 0 0
        Ringo 0 0 1
        >>> out = structure_population(pop, id=beatles,onerowperind=True)
        >>> print(out)
        John 0 1 0 0 1 1
        Paul 1 0 0 0 0 1
        George 0 1 0 0 1 1
        Ringo 1 0 0 0 0 1
    """
    table = []
    if isinstance(label, (list, tuple)):
        minsize = max(len(x) for x in label)
    elif isinstance(pop, PopulationBase):
        minsize = max(len(x.label) for x in pop)
    else:
        minsize = len(label) + len(str(len(pop)))

    for i, ind in enumerate(pop):
        if isinstance(label, (list, tuple)):
            row = [label[i]]
        else:
            row = [getattr(ind, 'label', label + str(i + 1))]
        row[0] = row[0].ljust(minsize + 1)
        table.append(row)

        # Add extra columns, if necessary
        for c in [popdata, popflag, locdata, phenotype]:
            if c is not None:
                row.append(str(c[i]))

        # One or two rows per individual?
        if onerowperind:
            row.extend(str(x) for x in ind.flatten())
        else:
            row_second = row[:]
            table.append(row_second)
            row.extend(str(x) for x in ind[:, 0])
            row_second.extend(str(x) for x in ind[:, 1])

    return linesep.join([sep.join(L) for L in table])


def prepare_setup_files(pop, save_files=True, save_dir=None,
                        id='ind', population_file='population.dat',
                        outfile='out',
                        popdata=None, popflag=None, locdata=None,
                        phenotype=None, onerowperind=True,
                        **kwargs):
    """
    Creates and return the 'mainparams', 'extraparams', and 'sample' files
    for a structure job.
    """

    if not pop:
        raise ValueError('cannot proceed with empty population.')

    def B(x): return int(bool(x))                                      # noqa: E731
    mainparams_data = mainparams(
        pop,
        infile=population_file,
        outfile=outfile,
        numinds=pop.shape[0],
        numloci=pop.shape[1],
        ploidy=pop.shape[2],
        popdata=B(popdata),
        popflag=B(popflag),
        locdata=B(locdata),
        phenotype=B(phenotype),
    )
    extraparams_data = extraparams(EXTRAPARAMS_DEFAULTS.copy())

    # Compute sample file
    popfile = structure_population(pop, id=label,
                                   onerowperind=onerowperind)

    # Save files in the specified directory
    if save_files:
        save_dir = os.path.abspath(save_dir) if save_dir else os.getcwd()

        def path(x): return os.path.join(save_dir, x)                  # noqa: E731

        with open(path('mainparams'), 'w', encoding='utf8') as F:
            F.write(mainparams_data)
        with open(path('extraparams'), 'w', encoding='utf8') as F:
            F.write(extraparams_data)
        with open(path(population_file), 'w', encoding='utf8') as F:
            F.write(popfile)

    return mainparams_data, extraparams_data, popfile


def mainparams(pop, infile='data.txt', **kwargs):
    """
    Creates data for a mainparams file.
    """

    data = dict(
        MAINPARAMS_DEFAULTS,
        numinds=pop.shape[0],
        numloci=pop.shape[1],
        ploidy=pop.shape[2],
        infile=infile,
    )
    data.update(**kwargs)
    return params_file(data)


def extraparams(pop, **kwargs):
    """
    Creates data for a mainparams file.
    """

    data = dict(EXTRAPARAMS_DEFAULTS)
    data.update(**kwargs)
    return params_file(data)


def params_file(data):
    """
    Creates a minimal, but functional, params file from a dictionary of
    arguments.
    """

    def normalize(x):
        if isinstance(x, bool):
            return int(x)
        return x

    lst = []
    for k, v in sorted(data.items()):
        lst.append('#define %s %s' % (k.upper(), normalize(v)))
    return '\n'.join(sorted(lst))

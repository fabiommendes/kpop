import io
import os
import subprocess
import tempfile
from subprocess import PIPE

from ..parsers.admixture import parse


def run_admixture(pop, k, job_dir=None, disp=1, supervised=False):
    """
    Runs the ADMIXTURE program.

    Args:
        pop (Population):
            Input population.
        k (int):
            Number of clusters/parental populations.
        job_dir:
            The directory used to store the result files. If not given,
            creates a temporary disposable directory.
        disp:
            The verbosity level in the 0-2 range.
        supervised (sequence):
            If given, must be a list of populations labels for each individual.
            Label should be None for individuals that have unknown ancestry.

    Returns:
        An :cls:`kpop.parsers.admixture.Admixture` object.
    """

    if job_dir is None:
        tempdir = tempfile.TemporaryDirectory()
        job_dir = tempdir.name
    else:
        job_dir = os.path.abspath(job_dir)

    # Create .ped and .map files
    with open(os.path.join(job_dir, 'job.ped'), 'w') as F:
        F.write(pop.io.plink_ped())
    with open(os.path.join(job_dir, 'job.map'), 'w') as F:
        F.write(pop.io.plink_map())

    # For supervised learning, create a .pop file
    if supervised:
        create_pop_file(os.path.join(job_dir, 'job.pop'), supervised)

    # Convert everything to a binary format since it seems that ADMIXTURE does
    # not support the regular .ped files (?)
    convert_to_bed('job', cwd=job_dir)

    # Set command line options to the admixture program
    cmd = ['admixture', 'job.bed', str(k)]
    if supervised:
        cmd.append('--supervised')
    if disp > 0:
        print('Running ADMIXTURE with command flags:\n'
              '  $ %s' % ' '.join(cmd))

    # Run the admixture program in a subprocess
    try:
        result = subprocess.run(cmd, stdout=PIPE, check=True, cwd=job_dir)
    except subprocess.CalledProcessError as ex:
        msg = 'ADMIXTURE ended with runtime code %s.\n\n' % ex.returncode
        msg += ex.stdout.decode('utf8')
        raise RuntimeError(msg)
    admix_out = io.StringIO(result.stdout.decode('utf8'))

    # Read the output of .P and .Q files with allele frequencies
    with open(os.path.join(job_dir, 'job.%s.P' % k)) as pfile, \
            open(os.path.join(job_dir, 'job.%s.Q' % k)) as qfile:
        result = parse(admix_out, qfile=qfile, pfile=pfile)

    return result


def create_pop_file(path, labels):
    """
    Creates a .pop file from the given list of labels.
    """
    with open(path, 'w') as F:
        for label in labels:
            if label is None:
                F.write('-\n')
            else:
                F.write('%s\n' % label)


def convert_to_bed(file, cwd):
    """
    Convert plink's .ped files to a .bed file.
    """

    cmd = ['plink', '--file', file, '--out', file, '--make-bed', '--noweb']
    try:
        subprocess.run(cmd, stdout=PIPE, check=True, cwd=cwd)
    except subprocess.CalledProcessError as ex:
        msg = 'plink ended with runtime code %s.\n\n' % ex.returncode
        msg += ex.stdout.decode('utf8')
        raise RuntimeError(msg)

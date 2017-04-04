import os
import subprocess
import tempfile
from subprocess import PIPE

import numpy as np

from .result import AdmixtureResult


def run_admixture(pop, k, job_dir=None, disp=1, supervised=False):
    """
    Runs the ADMIXTURE program.

    Args:
        pop:
            Input population.
        k:
            Number of clusters/parental populations.
        job_dir:
            The directory used to store the temporary files. If not given,
            creates a temporary disposable directory.
        disp:
            The verbosity level in the 0-2 range.
        supervised:
            If given, must be a list of populations labels for each individual.
            Label should be None for individuals that have unknown ancestry.

    Returns:
        A AdmixtureResult object.

    """
    if job_dir is None:
        tempdir = tempfile.TemporaryDirectory()
        job_dir = tempdir.name
    else:
        job_dir = os.path.abspath(job_dir)

    # Create .ped and .map files
    with open(os.path.join(job_dir, 'job.ped'), 'w') as F:
        F.write(pop.render_ped())

    # Creates .map file
    with open(os.path.join(job_dir, 'job.map'), 'w') as F:
        F.write(pop.render_map())

    # For supervised learning, create a .pop file
    if supervised:
        with open(os.path.join(job_dir, 'job.pop'), 'w') as F:
            for label in supervised:
                if label is None:
                    F.write('-\n')
                else:
                    F.write('%s\n' % label)

    # Convert everything to a binary format since it seems that ADMIXTURE does
    # not support the regular .ped files (?)
    cmd = ['plink', '--file', 'job', '--out', 'job', '--make-bed', '--noweb']
    try:
        subprocess.run(cmd, stdout=PIPE, check=True, cwd=job_dir)
    except subprocess.CalledProcessError as ex:
        msg = 'plink ended with runtime code %s.\n\n' % ex.returncode
        msg += ex.stdout.decode('utf8')
        raise RuntimeError(msg)

    # Set command line options to the admixture program
    cmd = ['admixture', 'job.bed', str(k)]
    if supervised:
        cmd.append('--supervised')
    if disp > 0:
        print('Running ADMIXTURE with command flags:\n'
              '  $ %s' % ' '.join(cmd))
    try:
        result = subprocess.run(cmd, stdout=PIPE, check=True, cwd=job_dir)
    except subprocess.CalledProcessError as ex:
        msg = 'ADMIXTURE ended with runtime code %s.\n\n' % ex.returncode
        msg += ex.stdout.decode('utf8')
        raise RuntimeError(msg)
    admix_out = result.stdout.decode('utf8')

    # Read the output .P file with allele frequencies
    with open(os.path.join(job_dir, 'job.%s.P' % k)) as F:
        data = []
        for line in F:
            data.append(list(map(float, line.strip().split())))
    freqs = np.array(data).T

    # Read the output .Q file with admixture coefficients
    with open(os.path.join(job_dir, 'job.%s.Q' % k)) as F:
        data = []
        for line in F:
            data.append(list(map(float, line.strip().split())))
    proportions = np.array(data)

    return AdmixtureResult(admix_out, freqs, proportions, ploidy=pop.ploidy)
